program ode_rk2
implicit none

real(8) :: m, k
real(8), dimension(2) :: x
real(8), dimension(2) :: k1, k2
real(8) :: en
real(8) :: t, dt
integer :: i, nt

m = 1.0d0
k = 1.0d0
dt = 0.01d0	!time step
nt = 1000
x(1) = 1.0d0	!initial position	
x(2) = 1.0d0	!initial velocity
t = 0.0d0	!initial time point
!en = 0.5

open(unit = 1, file = 'harmonic_osc_ode_rk2_results.dat', status = 'unknown')
do i = 1, nt
	en = 0.5d0 * k * x(1) ** 2 + 0.5d0 * m * x(2) ** 2
	write(1, *) t, x, en
	call rhs(m, k, x, k1)
	call rhs(m, k, x + 0.5d0 * k1 * dt, k2)
	x = x + dt * k2
	t = t + dt
	
end do
close(1)

end program ode_rk2


!---------------------------------------------------------

subroutine rhs(m, k, x, f)
implicit none

real(8), intent(in) :: m, k
real(8), dimension(2), intent(in) :: x
real(8), dimension(2), intent(out) :: f

f(1) = x(2)
f(2) = - (k/m) * x(1)
end subroutine rhs
