program taylor_exp
	implicit none

	real :: x, f_true, f
	real, dimension(5) :: remi
	integer :: i
	real :: h

	x = 0
	h = 1
	
	call funct(1.0,f_true)
	
	print*, "Actual f(0) =", f_true
	
	call tay(x, h, f)

	print*,"By Taylor expansion:", f
	
	call reminder(x, h, f, remi)
	
	print*,"Reminders:"
	
	do i = 1, 5
		print*, "n =", i,"-->", remi(i)
	end do

end program taylor_exp




subroutine funct(x,f)
	implicit none
	real, intent(in) :: x
	real, intent(out) :: f
	
	f = -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x +1.2
	
end subroutine funct




subroutine tay(x,h,f)
	implicit none
	real, intent(in) :: x, h
	real, intent(out) :: f
	real :: f0, f1, f2, f3, f4
	
	f1 = -0.4*x**3 - 0.45*x*x - x - 0.25
	f2 = -1.2*x*x - 0.9*x - 1
	f3 = -2.4*x - 0.9
	f4 = -2.4
	
	call funct(0.0, f0)
	
	f = f0 + f1*h + f2*0.5*h**2 + f3*h**3/6 + f4*h**4/24 
	
end subroutine tay




subroutine reminder(x, h, f, rem)
	implicit none
	real, intent(in) :: x, h
	real, intent(in) :: f
	real, dimension(5), intent(out) :: rem
	integer :: i
	real :: f0, f1, f2, f3, f4
	
	f2 = -1.2*x*x - 0.9*x - 1
	f3 = -2.4*x - 0.9
	f4 = -2.4
	
	
	rem(1) = f2*h**2/2
	rem(2) = f3*h**3/6
	rem(3) = f4*h**4/24
	
	do i = 4, 5
		rem(i) =  0
	end do
	
end subroutine reminder
