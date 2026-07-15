program ising2d
    implicit none

    ! Parameters for the Ising model
    integer, parameter :: L = 32        ! Lattice size (L x L)
    integer, parameter :: steps = 100000 ! Total Monte Carlo steps
    real(8), parameter :: T = 0.0       ! Temperature (in units of kT/Jval)

    ! Variables
    integer :: i, j, step               ! Loop counters
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
    
    ! Random number initialization
    call random_seed()

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
    
    
	OPEN(UNIT = 1, FILE = 'ini_ising_re.dat', STATUS = 'UNKNOWN')
        	WRITE(1,*) (L*L)
        	WRITE(1,*)
        	DO i = 1, L
                	DO j =1, L
                		WRITE(1,*)i, j, lattice(i, j)
                	END DO
        	END DO 
    	CLOSE(1)
    	

    ! Monte Carlo simulation loop
    
    OPEN(UNIT = 2, FILE = 'mont_ising_re.dat', STATUS = 'UNKNOWN')
    
    
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
            lattice(x, y) = -lattice(x, y)  ! DEBUG --> DONE
            energy = energy + deltaE        ! Update energy
            magnetization = magnetization + 2.0 * real(lattice(x, y))
        else
            call random_number(r)
            prob = exp(-deltaE / (kb * T))
            if (r < prob) then
                lattice(x, y) = lattice(x, y)  ! DEBUG --> DONE
                energy = energy + deltaE        ! Update energy
                magnetization = magnetization + 2.0 * real(lattice(x, y))
            end if
        end if
    !end do
    
    WRITE(2,*)(L*L)
        WRITE(2,*)
        DO i = 1, L
                DO j = 1, L
                	WRITE(2,*)i ,j, lattice(i, j)
                END DO
        END DO
      
    end do
    
    
    CLOSE(2)
    
    

    ! Print final results
    print*, "At tempratrure =", T, ":"
    print *, "Final Energy per site: ", energy / (L * L)
    print *, "Final Magnetization per site: ", magnetization/(L * L)
    
end program ising2d
