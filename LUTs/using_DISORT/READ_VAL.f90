      SUBROUTINE READ_VAL( FNAME,COT,MU0,PMOM)
! c To read input values from a text file
      CHARACTER*(*)      FNAME
      REAL            COT, MU0
      INTEGER K
      REAL, DIMENSION (900) :: PMOM
      OPEN ( 10, FILE=FNAME )
      READ ( 10, *) COT
      READ ( 10, *) MU0
      DO  60  K = 1, 900
         READ (10,*) PMOM(K)
60    CONTINUE
      CLOSE(10)

      END
