PROGRAM num_deri
	IMPLICIT NONE

	INTEGER :: i, n
	REAL :: h
	REAL, DIMENSION(5) :: x, f
	REAL, DIMENSION(5) :: deri_for, deri_back, deri_cen
	n  = 5
	h = 0.5
	x = (/0.0, 0.5, 1.0, 1.5, 2.0/)
	
	DO i = 1, n
		f(i) = -0.1*x(i)**4 - 0.15*x(i)**3 - 0.5*x(i)*x(i) - 0.25*x(i) + 1.2
	END DO
	
	CALL forward(n, x, h, f, deri_for)
	CALL backward(n, x, h, f, deri_back)
	CALL central(n, x, h, f, deri_cen)
	
	PRINT*, "Forward:"
	DO i = 1, n
		PRINT*, x(i), deri_for(i)
	END DO
	
	PRINT*, "Backward:"
	DO i = 1, n
		PRINT*, x(i), deri_back(i)
	END DO
	
	PRINT*, "Central:"
	DO i = 1, n
		PRINT*, x(i), deri_cen(i)
	END DO
	! actual derivative at x = 0.5 is -0.9125	
	PRINT*, "Differencr from actual derivative at x = 0.5:"
	PRINT*, "Forward:", - 0.9123 - deri_for(2), "Backward:", - 0.9123 - deri_back(2),"Central:", - 0.9123 - deri_cen(2)
	
	
END PROGRAM num_deri

SUBROUTINE forward(n, x, h, f, deri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	REAL, INTENT(IN) :: h
	REAL, DIMENSION(n), INTENT(IN) :: x, f
	REAL, DIMENSION(n), INTENT(OUT) :: deri
	INTEGER :: i
	
	DO  i = 1,  4
		deri(i) = (f(i+1) - f(i))/h
	END DO
	
	deri(n) = deri(n-1)
	
END SUBROUTINE

SUBROUTINE backward(n, x, h, f, deri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	REAL, INTENT(IN) :: h
	REAL, DIMENSION(n), INTENT(IN) :: x, f
	REAL, DIMENSION(n), INTENT(OUT) :: deri
	INTEGER :: i
	
	DO i = n, 2, -1
		deri(i) = (f(i) - f(i - 1))/h
	END DO
	
	deri(1) = deri(2)
END SUBROUTINE backward

SUBROUTINE central(n, x, h, f, deri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	REAL, INTENT(IN) :: h
	REAL, DIMENSION(n), INTENT(IN) :: x, f
	REAL, DIMENSION(n), INTENT(OUT) :: deri
	INTEGER :: i
	
	DO i = 2, n-1
		deri(i) = (f(i + 1) - f(i - 1))/(2*h)
	END DO
	
	deri(n) = deri(n-1)
END SUBROUTINE central
