!Find out the root of the following equation by using Newton Raphson method:

!m = tanh(m/T)

!Consider the values of T from 0.1 to 1.1 in the steps of 0.1. Plot root Vs T.


program newt_raph
implicit none

real(8) :: step, a, b, eps, f, deri, step_len, ini
real(8), dimension(10) :: T, m
integer :: i, step_len_int, iter, max

m = 0.0

step = 0.1

!Limits:
a = 0.1 
b = 1.1 

eps = 1e-8

!Initial guess:
ini = 1.0

!Maximum iterations to prevent infinite loop:
max = 100

step_len = (b - a)/step + 1
step_len_int = int(step_len)

do i = 1, step_len_int
	T(i) = a + (i-1) * step
end do

!Newton-Raphson Steps:
do i = 1, step_len_int
	m(i) = ini
	do iter = 1, max
		call func(m(i), T(i), f)
		call deriv(m(i), T(i), deri)
		
		if(abs(f) < eps .or. abs(deri) < eps) exit !preventing division by very small or zero derivative
			m(i) = m(i) - f/deri
	end do
	
end do

open(unit = 1, file = 'new_raph_root_vs_T.txt', status = 'unknown')
do i = 1, step_len_int
	write(1, *) T(i), m(i)
end do	

close(1)

end program newt_raph





subroutine func(m, T, f)
implicit none

real(8), intent(in) :: m, T
real(8), intent(out) :: f

f = tanh(m/T) - m

end subroutine func



subroutine deriv(m, T, deri)

real(8), intent(in) :: m, T
real(8), intent(out) :: deri

deri = (1 / (T * cosh(m/T)**2)) - 1

end subroutine deriv

