PROGRAM compare
        IMPLICIT NONE
        REAL  :: a, b, c
        READ*, a, b, c

        IF(a>b) THEN
                IF(a>c) THEN
                        PRINT*,"The largest number is:", a
                ELSE
                        PRINT*,"The largest number is:",c
                ENDIF
        ELSE
                IF(b>c) THEN
                        PRINT*, "The largest number is:", b
                ELSE
                        PRINT*, "The largest number is :", c
                ENDIF
        ENDIF
END PROGRAM compare
