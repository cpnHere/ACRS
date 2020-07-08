      SUBROUTINE READ_VAL( FNAME,RVAL)
      CHARACTER*(*)      FNAME
      REAL            RVAL
      OPEN ( 10, FILE=FNAME )
      READ ( 10, *) RVAL
      END
