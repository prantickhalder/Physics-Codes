PROGRAM temp
        IMPLICIT NONE
        REAL :: c, f
        REAL :: t1, t2
        PRINT*, "Enter the temperature in Fahrenheit:"
        READ*, f

        CALL CPU_TIME(t1)
        c = 5.0 * (f - 32.0)/9.0
        CALL CPU_TIME(t2)

        PRINT*, "The temperature in Celcius is:", c, "Run time:", t2-t1

        END PROGRAM temp
