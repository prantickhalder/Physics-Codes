PROGRAM quad
        IMPLICIT NONE
        REAL :: a, b , c, r1, r2, d, re, im1, im2
        PRINT*, "Enter coefficients a, b and c in ax^2 + bx + c = 0:"
        READ*, a, b, c

        d = b*b - 4.0*a*c

        IF(d > 0.0) THEN 
                r1= 0.5*(-b + SQRT(d))/a
                r2 = 0.5*(-b - SQRT(d))/a
                PRINT*, "The roots are real and distinct and they are:", r1, r2
        ELSE IF(d == 0.0) THEN
                r1 = 0.5*(-b + SQRT(d))/a
                PRINT*, "The roots are real and equal and they are:", r1, r1
        ELSE
                re = 0.5 * (-b)/a
                im1 = 0.5 * ((-1)*SQRT(-d))/a
                im2 = 0.5 * (-(-1)*SQRT(-d))/a

                PRINT*,"The roots are imaginary"
                PRINT*, "The roots are:",re,"+",im1,"i","   and",re,"+",im2,"i"
                
        END IF

END PROGRAM quad
