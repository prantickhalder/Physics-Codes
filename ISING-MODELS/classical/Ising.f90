program ising2d
    implicit none

    ! Parameters for the Ising model
    integer, parameter :: L = 32        ! Lattice size (L x L)
    integer, parameter :: steps = 1000000 ! Total Monte Carlo steps
    integer, parameter :: n = 10 ! for array length
    
    
    real(8) :: T
    
    real(8), dimension(n) :: mag ! to store magnetisation per site

    ! Variables
    integer :: i, j, step, k               ! Loop counters
    integer :: lattice(L, L)           ! 2D lattice for spins
    integer :: x, y                    ! Randomly selected spin coordinates
    integer :: neighbor_sum            ! Sum of neighbor spins
    real(8) :: energy, magnetization   ! Energy and magnetization
    real(8) :: deltaE                  ! Energy change for a spin flip
    real(8) :: prob                    ! Probability for accepting a flip
    real(8) :: r                       ! Random number
    real(8), parameter :: kb = 1.0     ! Boltzmann constant (set to 1)
    real(8), parameter :: Jval = 1.0      ! Interaction strength (set to 1)
    integer(8) :: up, down, left , right
    
    real(8), dimension(n) :: temp = (/0.0, 0.5, 1.0, 1.5, 2.0, 2.26, 2.5, 3.0, 3.5, 4.0/) !(in units of kT/Jval)
    
    
    
    ! Random number initialization
    call random_seed()
    
    
    ! Open file to store magnetization vs. temperature
    open(unit=3, file="mag_vs_temp.dat", status="unknown")
    write(3,*) "Temperature Magnetization_per_site"

    ! Loop over temperature values
    do k = 1, n
        T = temp(k)

    ! Initialize the lattice with random spins (+1 or -1)
    do i = 1, L
        do j = 1, L
            call random_number(r)
            if (r < 0.5) then ! DEBUG --> DONE
                lattice(i, j) = 1
            else
                lattice(i, j) = -1
            end if
        end do
    end do

    ! Calculate initial energy and magnetization
    energy = 0.0
    magnetization = 0.0
    do i = 1, L
        do j = 1, L
        
        
        if(i == 1) then
        	up = lattice(L, j)
        else
        	up = lattice(i - 1, j)
        end if
        
        
        if(i == L) then
        	down = lattice(1, j)
        else
        	down = lattice(i + 1, j)
        end if
        
        
        if(j == 1) then
        	left = lattice(i, L)
        else
        	left = lattice(i, j - 1)
        end if
        
        
        if(j == L) then
        	right = lattice(i, 1)
        else
        	right = lattice(i, j + 1)
        end if
        
        
            neighbor_sum = (up + down + left + right)      ! DEBUG --> DONE
            energy = energy - 0.5 * Jval * real(lattice(i, j)) * real(neighbor_sum) !QUESTION - Why 0.5 factor --> since every lattice interaction was considered twice
            magnetization = magnetization + real(lattice(i, j))
        end do
    end do
    
    
    ! Monte Carlo simulation loop
       
    do step = 1, steps
        ! Randomly select a spin
        call random_number(r)
        x = int(r * L) + 1
        call random_number(r)
        y = int(r * L) + 1
        
        
	if(x == 1) then
        	up = lattice(L, y)
        else
        	up = lattice(x - 1, y)
        end if
        
        
        if(x == L) then
        	down = lattice(1, y)
        else
        	down = lattice(x + 1, y)
        end if
        
        
        if(y == 1) then
        	left = lattice(x, L)
        else
        	left = lattice(x, y - 1)
        end if
        
        
        if(y == L) then
        	right = lattice(x, 1)
        else
        	right = lattice(x, y + 1)
        end if

        ! Calculate the change in energy if the spin is flipped
        neighbor_sum = (up + down + left + right) ! DEBUG --> DONE
        deltaE = 2.0 * Jval * lattice(x, y) * neighbor_sum !IS IT CORRECT? --> Yes --> 2.0 due to flipping since 1-(-1)) = 2

        ! Decide whether to flip the spin
        if (deltaE <= 0.0) then
            magnetization = magnetization - 2.0 * real(lattice(x, y))
            lattice(x, y) = -lattice(x, y)  ! DEBUG --> DONE
            energy = energy + deltaE        ! Update energy
            
        else
            call random_number(r)
            prob = exp(-deltaE / (kb * T))
            if (r < prob) then
            	magnetization = magnetization - 2.0 * real(lattice(x, y))
                lattice(x, y) = -lattice(x, y)  ! DEBUG --> DONE
                energy = energy + deltaE        ! Update energy
               
            end if
        end if
             
    end do
    
    
    mag(k) = magnetization/(L * L)
    
    write(3,*) T, abs(mag(k))
    

    ! Print final results
    print*, "At temperature =", T, "Final Magnetization per site: ", abs(mag(k))

    end do
    
    close(3)
        
end program ising2d
