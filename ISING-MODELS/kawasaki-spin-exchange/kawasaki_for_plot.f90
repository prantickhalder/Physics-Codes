program ising_kawasaki
    implicit none

    ! Parameters
    integer, parameter :: L = 32
    integer, parameter :: steps = 300000
    real(8), parameter :: T = 2.5
    real(8), parameter :: kb = 1.0
    real(8), parameter :: Jval = 1.0

    ! Variables
    integer :: i, j, step
    integer :: lattice(L, L)
    integer :: x1, y1, x2, y2, dir
    integer :: up, down, left, right
    integer :: neighbor_sum
    real(8) :: energy, magnetization
    real(8) :: deltaE, r
    real(8) :: e_before, e_after

    ! Random seed
    call random_seed()

    ! Initialize lattice with random spins
    do i = 1, L
        do j = 1, L
            call random_number(r)
            if (r < 0.5) then
                lattice(i, j) = 1
            else
                lattice(i, j) = -1
            end if
        end do
    end do

    ! Initial energy and magnetization
    energy = 0.0
    magnetization = 0.0
    do i = 1, L
        do j = 1, L
            up    = lattice(mod(i - 2 + L, L) + 1, j)
            down  = lattice(mod(i, L) + 1, j)
            left  = lattice(i, mod(j - 2 + L, L) + 1)
            right = lattice(i, mod(j, L) + 1)

            neighbor_sum = up + down + left + right
            energy = energy - 0.5 * Jval * lattice(i, j) * neighbor_sum
            magnetization = magnetization + lattice(i, j)
        end do
    end do

    ! Save initial configuration
    open(unit = 1, file = 'ini_ising_re.dat', status = 'unknown')
        write(1,*) L * L
        write(1,*)
        do i = 1, L
            do j = 1, L
                write(1,*) i, j, lattice(i, j)
            end do
        end do
    close(1)

    ! Open output file for configuration updates
    open(unit = 2, file = 'mont_ising_re.dat', status = 'unknown')

    ! Kawasaki Monte Carlo loop
    do step = 1, steps
        call random_number(r)
        x1 = int(r * L) + 1
        call random_number(r)
        y1 = int(r * L) + 1

        call random_number(r)
        dir = int(r * 4) + 1

        select case(dir)
        case(1)  ! down
            x2 = mod(x1, L) + 1
            y2 = y1
        case(2)  ! up
            x2 = mod(x1 - 2 + L, L) + 1
            y2 = y1
        case(3)  ! left
            x2 = x1
            y2 = mod(y1 - 2 + L, L) + 1
        case(4)  ! right
            x2 = x1
            y2 = mod(y1, L) + 1
        end select

        ! Swap only if spins differ
        if (lattice(x1, y1) /= lattice(x2, y2)) then
            e_before = -Jval * ( &
                lattice(x1, y1) * ( &
                    lattice(mod(x1, L) + 1, y1) + lattice(mod(x1 - 2 + L, L) + 1, y1) + &
                    lattice(x1, mod(y1, L) + 1) + lattice(x1, mod(y1 - 2 + L, L) + 1) ) + &
                lattice(x2, y2) * ( &
                    lattice(mod(x2, L) + 1, y2) + lattice(mod(x2 - 2 + L, L) + 1, y2) + &
                    lattice(x2, mod(y2, L) + 1) + lattice(x2, mod(y2 - 2 + L, L) + 1) ) )

            call swap(lattice(x1, y1), lattice(x2, y2))

            e_after = -Jval * ( &
                lattice(x1, y1) * ( &
                    lattice(mod(x1, L) + 1, y1) + lattice(mod(x1 - 2 + L, L) + 1, y1) + &
                    lattice(x1, mod(y1, L) + 1) + lattice(x1, mod(y1 - 2 + L, L) + 1) ) + &
                lattice(x2, y2) * ( &
                    lattice(mod(x2, L) + 1, y2) + lattice(mod(x2 - 2 + L, L) + 1, y2) + &
                    lattice(x2, mod(y2, L) + 1) + lattice(x2, mod(y2 - 2 + L, L) + 1) ) )

            deltaE = e_after - e_before

            if (T == 0.0) then
                if (deltaE > 0.0) then
                    call swap(lattice(x1, y1), lattice(x2, y2)) ! Reject
                else
                    energy = energy + deltaE
                end if
            else
                if (deltaE <= 0.0) then
                    energy = energy + deltaE
                else
                    call random_number(r)
                    if (r < exp(-deltaE / (kb * T))) then
                        energy = energy + deltaE
                    else
                        call swap(lattice(x1, y1), lattice(x2, y2)) ! Reject
                    end if
                end if
            end if
        end if

        ! Optional: write only every 1000 steps
        if (mod(step, 1) == 0) then
            write(2,*) L * L
            write(2,*)
            do i = 1, L
                do j = 1, L
                    write(2,*) i, j, lattice(i, j)
                end do
            end do
        end if
    end do

    close(2)

    ! Final magnetization
    magnetization = sum(lattice)

    ! Output results
    print*, "At temperature =", T
    print*, "Final Energy per site: ", energy / (L * L)
    print*, "Final Magnetization per site: ", magnetization / (L * L)

contains
    subroutine swap(a, b)
        integer, intent(inout) :: a, b
        integer :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap

end program ising_kawasaki

