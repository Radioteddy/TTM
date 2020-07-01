        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 01 14:16:30 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TEST__genmod
          INTERFACE 
            PURE ELEMENTAL FUNCTION TEST(X,Y) RESULT(Z)
              REAL(KIND=8), INTENT(IN) :: X
              REAL(KIND=8), INTENT(IN) :: Y
              REAL(KIND=8) :: Z
            END FUNCTION TEST
          END INTERFACE 
        END MODULE TEST__genmod
