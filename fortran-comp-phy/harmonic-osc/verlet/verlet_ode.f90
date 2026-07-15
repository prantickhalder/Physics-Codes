program ode_verlet
implicit none

real(8) :: m, k, mi, a, a_new
real(8), dimension(1000) :: x, v ! nt = 1000
real(8) :: en
real(8) :: t, dt
integer :: i, nt

m = 1.0d0
k = 1.0d0
dt = 0.01d0
nt = 1000	!n*t
x(1) = 1.0d0	!initial position	
v(1) = 1.0d0	!initial velocity
t = 0.0d0	!initial time
!en = 0.5d0	!energy at initial conditions

mi = 1.0d0/m	!mass inverse

open(unit = 1, file = 'harmonic_osc_ode_VERLET_x_vs_t.dat', status = 'unknown')
do i = 1, (nt - 1)
	a = - k * x(i) * mi
	x(i + 1) = x(i) + v(i) * dt + 0.5d0 * a * dt * dt	!F = -kx
	a_new = - k * x(i + 1) * mi
	v(i + 1) = v(i) + 0.5d0 * dt * (a_new + a)
	t = t + dt
	en = 0.5d0 * k * x(i + 1) * x(i + 1) + 0.5d0 * m * v(i + 1) * v(i + 1)	
	write(1, *) t, x(i + 1), v(i + 1), en
end do
close(1)

end program ode_verlet


!---------------------------------------------------------

