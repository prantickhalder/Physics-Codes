PROGRAM facto
        IMPLICIT NONE
        INTEGER(kind=16) :: n, f, i
        

        PRINT*,"Enter a number:"
        READ*, n
        f=1

        DO i = n,1,-1
                f= f*i
        ENDDO

        PRINT*,"The facatorial is:", f

END PROGRAM facto
