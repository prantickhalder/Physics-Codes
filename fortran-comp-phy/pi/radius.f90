PROGRAM RADIUS
        IMPLICIT NONE
        
        REAL :: p(4), M1, M2, x(500000), y(500000) 
        INTEGER :: i, j,  N(4), fu
        
        CHARACTER (len = *), PARAMETER :: OUT = 'pi_data.txt'
        CHARACTER (len = *), PARAMETER :: IN = 'pi_plot.plt'

        !PRINT*, "Enter  N:"
        !READ*, N
        
        N = (/10000, 50000, 100000, 500000/)
        M1 = 0.0
        M2 = 0.0

        DO j =1, 4 

        DO i = 1, N(j)
              
                CALL RANDOM_NUMBER(x(i))
                CALL RANDOM_NUMBER(y(i))
                x(i) = x(i) - 0.5 ! scaling to [-0.5, 0.5] since it is not centered at (0,0) 
                y(i) = y(i) - 0.5

                IF(x(i)*x(i) + y(i)*y(i) <= 0.25) THEN
                M1 = M1 + 1.0
                ENDIF
                

                M2 = M2 + 1.0

                ENDDO
        
        p(j) = 4.0 * (M1/M2)
        
        ENDDO

        OPEN (ACTION = 'write', FILE = OUT, NEWUNIT = fu, STATUS = 'replace')

        DO j = 1, 4
                WRITE(fu, *) N(j), p(j)
        ENDDO

        CLOSE(fu)


        PRINT*, "PI = ", p

        
        
END PROGRAM RADIUS





