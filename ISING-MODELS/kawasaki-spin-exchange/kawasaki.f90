program kawasaki
    implicit none

    ! Parameters
    integer, parameter :: L = 32
    integer, parameter :: steps = 1000000
    integer, parameter :: n = 10
    integer, parameter :: equil_steps = 100000 !this can be varied
    integer, parameter :: sample_interval = 100 !interval length
    real(8), parameter :: kb = 1.0
    real(8), parameter :: Jval = 1.0

    ! Temperature and results
    real(8), dimension(n) :: temp = (/0.0, 0.5, 1.0, 1.5, 2.0, 2.26, 2.5, 3.0, 3.5, 4.0/)
    real(8), dimension(n) :: mag, energy_avg, cv
    real(8) :: T

    ! Variables
    integer :: i, j, k, step, num_samples
    integer :: x1, y1, x2, y2, dir
    integer :: lattice(L, L)
    integer :: up, down, left, right, neighbor_sum
    real(8) :: energy, deltaE, r
    real(8) :: e_before, e_after
    real(8) :: energy_sum, energy2_sum, mag_sum, magnetization
    character(len=40) :: filename
    integer :: unit_step_file

    call random_seed()

    open(unit=3, file="results_kawasaki.dat", status="unknown")
    write(3,*) "T   Mag_per_site   E_per_site   Cv_per_site"

    do k = 1, n
        T = temp(k)

        ! Random initial lattice
        do i = 1, L
            do j = 1, L
                call random_number(r)
                !Condensed form used instead of if-else statements:
                lattice(i,j) = merge(1, -1, r < 0.5) !condensed form 
            end do
        end do

        ! Compute initial energy
        energy = 0.0
        do i = 1, L
            do j = 1, L
            !Periodic boundary conditions:
                up    = lattice(mod(i - 2 + L, L) + 1, j)
                down  = lattice(mod(i, L) + 1, j)
                left  = lattice(i, mod(j - 2 + L, L) + 1)
                right = lattice(i, mod(j, L) + 1)
                neighbor_sum = up + down + left + right
                energy = energy - 0.5 * Jval * lattice(i, j) * neighbor_sum
            end do
        end do

        ! Initializing:
        mag_sum = 0.0
        energy_sum = 0.0
        energy2_sum = 0.0

        ! Opening file to store magnetization vs steps
        write(filename, '("mag_vs_steps_T", F4.2, ".dat")') T
        filename = adjustl(filename)
        open(newunit=unit_step_file, file=trim(filename), status="unknown")
        write(unit_step_file, *) "# Step    Magnetization_per_site"

        ! Kawasaki dynamics
        do step = 1, steps
            call random_number(r); x1 = int(r * L) + 1
            call random_number(r); y1 = int(r * L) + 1
            call random_number(r); dir = int(r * 4) + 1
            !the following is used instead of if-else statements
            select case(dir)
            case(1); x2 = mod(x1, L) + 1;       y2 = y1       ! down
            case(2); x2 = mod(x1 - 2 + L, L)+1; y2 = y1       ! up
            case(3); x2 = x1;                   y2 = mod(y1 - 2 + L, L)+1 ! left
            case(4); x2 = x1;                   y2 = mod(y1, L)+1         ! right
            end select

            if (lattice(x1,y1) /= lattice(x2,y2)) then
                e_before = -Jval * ( &
                    lattice(x1, y1) * (lattice(mod(x1,L)+1,y1) + lattice(mod(x1-2+L,L)+1,y1) + &
                                       lattice(x1,mod(y1,L)+1) + lattice(x1,mod(y1-2+L,L)+1)) + &
                    lattice(x2, y2) * (lattice(mod(x2,L)+1,y2) + lattice(mod(x2-2+L,L)+1,y2) + &
                                       lattice(x2,mod(y2,L)+1) + lattice(x2,mod(y2-2+L,L)+1)) )

                call swap(lattice(x1,y1), lattice(x2,y2))

                e_after = -Jval * ( &
                    lattice(x1, y1) * (lattice(mod(x1,L)+1,y1) + lattice(mod(x1-2+L,L)+1,y1) + &
                                       lattice(x1,mod(y1,L)+1) + lattice(x1,mod(y1-2+L,L)+1)) + &
                    lattice(x2, y2) * (lattice(mod(x2,L)+1,y2) + lattice(mod(x2-2+L,L)+1,y2) + &
                                       lattice(x2,mod(y2,L)+1) + lattice(x2,mod(y2-2+L,L)+1)) )

                deltaE = e_after - e_before

                if (T == 0.0) then
                    if (deltaE > 0.0) call swap(lattice(x1,y1), lattice(x2,y2))  ! Reject
                    if (deltaE <= 0.0) energy = energy + deltaE
                else
                    if (deltaE <= 0.0) then
                        energy = energy + deltaE
                    else
                        call random_number(r)
                        if (r < exp(-deltaE / (kb * T))) then
                            energy = energy + deltaE
                        else
                            call swap(lattice(x1,y1), lattice(x2,y2))  ! Reject
                        end if
                    end if
                end if
            end if

            ! Sampling
            if (step > equil_steps .and. mod(step, sample_interval) == 0) then
                magnetization = sum(lattice)
                mag_sum = mag_sum + abs(magnetization / (L*L))
                energy_sum = energy_sum + energy
                energy2_sum = energy2_sum + energy**2

                ! Storing magnetization vs step
                write(unit_step_file, '(I10, F14.8)') step, magnetization / (L*L)
            end if
        end do

        num_samples = (steps - equil_steps) / sample_interval
        mag(k) = mag_sum / num_samples
        energy_avg(k) = energy_sum / num_samples / (L*L)
        cv(k) = (energy2_sum / num_samples - (energy_sum / num_samples)**2) / (kb * T**2 * L * L)

        write(3,'(F5.2,3X,F10.6,3X,F10.6,3X,F10.6)') T, mag(k), energy_avg(k), cv(k)
        print*, "T =", T, "Mag =", mag(k), "E/site =", energy_avg(k), "Cv/site =", cv(k)

        close(unit_step_file)
    end do

    close(3)
end program kawasaki


subroutine swap(a, b)
        integer, intent(inout) :: a, b
        integer :: temp
        temp = a
        a = b
        b = temp
end subroutine swap
