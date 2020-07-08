      SUBROUTINE READ_VAL( FNAME,COT,MU0)
! c To read input values from a text file
      CHARACTER*(*)      FNAME
      REAL            COT, MU0
      OPEN ( 10, FILE=FNAME )
      READ ( 10, *) COT
      READ ( 10, *) MU0
      CLOSE (10)
      END
