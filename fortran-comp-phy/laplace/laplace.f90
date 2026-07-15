program laplace
implicit none

real(8), dimension(:, :), allocatable :: phi, old_phi
real(8) :: eps, tol
integer :: i, j, n

n = 50 ! grid of dimension n by n

allocate(phi(n, n), old_phi(n, n))

phi = 0.0d0

! Boundary condition: only line y=1 i.e., phi(i, n) = 1, rest(y=0, x=1, and x=0) are 0
phi(:, n) = 1.0d0

tol = 1.0e-5
eps = 1.0d0

do while(eps > tol)

	old_phi = phi
	
	do i = 2, n-1
		do j =2, n-1
			phi(i, j) = 0.25d0 * (phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + phi(i, j-1))
						
		end do
	end do


	eps = maxval(abs(phi - old_phi))
end do

open(unit = 1, file = 'laplace.dat', status = 'unknown')
do i = 1, n
	write(1, *) (phi(i, j), j = 1, n)
end do
close(1)

end program laplace
