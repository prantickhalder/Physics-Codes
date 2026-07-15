program ode_euler
implicit none

real(8) :: m, k
real(8), dimension(2) :: x
real(8), dimension(2) :: f
real(8) :: en
real(8) :: t, dt
integer :: i, nt

m = 1.0d0
k = 1.0d0
dt = 0.1d0
nt = 1000
x(1) = 1.0d0	!initial position	
x(2) = 1.0d0	!initial velocity
t = 0.0d0
en = 0.5

open(unit = 1, file = 'harmonic_osc_ode_euler_x_vs_t.dat', status = 'unknown')
do i = 1, nt
	write(1, *) t, x, en
	call rhs(m, k, x, f)
	x = x + dt*f
	t = t + dt
	en = 0.5d0 * k * x(1) ** 2 + 0.5 * m * x(2) ** 2
end do
close(1)

end program ode_euler


!---------------------------------------------------------

subroutine rhs(m, k, x, f)
implicit none

real(8), intent(in) :: m, k
real(8), dimension(2), intent(in) :: x
real(8), dimension(2), intent(out) :: f

f(1) = x(2)
f(2) = - (k/m) * x(1)
end subroutine rhs
