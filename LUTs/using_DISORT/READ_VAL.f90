      SUBROUTINE READ_VAL( FNAME,COT,SCA,MU0,PMOM)
! c To read input values from a text file
      CHARACTER*(*)      FNAME
      REAL            COT, MU0, SCA
      INTEGER K
      REAL, DIMENSION (500) :: PMOM
      OPEN ( 10, FILE=FNAME )
      READ ( 10, *) COT
      READ ( 10, *) SCA 
      READ ( 10, *) MU0
      DO  60  K = 1, 500
         READ (10,*) PMOM(K)
60    CONTINUE
      CLOSE(10)

      END
