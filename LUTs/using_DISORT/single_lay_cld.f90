module disort_vars
implicit none 

LOGICAL DOPROB( 17 )
DATA DOPROB / 17*.TRUE. /
LOGICAL  USRANG, USRTAU, ONLYFL, PRNT(5), &
         PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE
INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
INTEGER  NUMU_O
LOGICAL  DEBUG
real(kind=4) :: ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, &
                PHI0, TEMIS, TTEMP, WVNMLO, WVNMHI, UMU0 
real(kind=4),parameter :: EARTH_RADIUS = 6371.0     
DATA PRNT  / .TRUE., 3*.FALSE., .TRUE. /

INTEGER       BRDF_TYPE
REAL          BRDF_ARG(4)
LOGICAL       DO_SHADOW
REAL          WIND_SPD, REFRAC_INDEX
REAL          B0, HH, W
REAL          K_VOL, K_ISO, K_GEO
REAL          RHO_0, KAPPA, G, H0
REAL          FLUX_UP, DFDTAU
INTEGER       NMUG
REAL          BDREF
EXTERNAL      BDREF

real(kind=4),dimension(:),allocatable     :: DTAUC, PHI, SSALB, TEMPER, UMU, UTAU                             
real(kind=4),dimension(:,:),allocatable   :: PMOM          
real(kind=4),dimension(:,:,:),allocatable :: RHOQ, RHOU 
real(kind=4),dimension(:),allocatable     :: EMUST, BEMST   
real(kind=4),dimension(:,:),allocatable   :: RHO_ACCURATE                             
real(kind=4),dimension(:),allocatable     :: RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED
real(kind=4),dimension(:,:,:),allocatable :: UU
real(kind=4),dimension(:),allocatable     :: H_LYR

contains 

subroutine allocate_disort_allocatable_arrays(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)
implicit none
integer,intent(in) :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU 

allocate( DTAUC( NLYR ), SSALB( NLYR ), PMOM( 0:NMOM, NLYR ), &
          TEMPER( 0:NLYR ), UTAU( NTAU ), UMU( NUMU ), PHI( NPHI ), H_LYR( 0:NLYR ) )  
allocate( RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
          EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI) )                
allocate( RFLDIR( NTAU ), RFLDN( NTAU ), FLUP( NTAU ), DFDT( NTAU ), UAVG( NTAU ),&
          ALBMED( NUMU ), TRNMED( NUMU ), UU( NUMU, NTAU, NPHI ) )   
DTAUC = 0.0; SSALB = 0.0; PMOM = 0.0; TEMPER = 0.0; UTAU = 0.0; UMU = 0.0; PHI = 0.0;
H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0;
RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; DFDT = 0.0; UAVG = 0.0; UU = 0.0;
ALBMED = 0.0; TRNMED = 0.0; 

end subroutine allocate_disort_allocatable_arrays

subroutine deallocate_disort_allocatable_arrays()

deallocate( DTAUC, SSALB, PMOM, TEMPER, UTAU, UMU, PHI, H_LYR )  
deallocate( RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE )                
deallocate( RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED, UU )  

end subroutine deallocate_disort_allocatable_arrays

end module disort_vars


program SINGLE_LAY_CLD
use disort_vars
implicit none

!** Correct answers
INTEGER  MXPROB, MXCASE, MXTAU, MXMU, MXPHI
PARAMETER  ( MXPROB = 16, MXCASE = 8, MXTAU = 5,  &
             MXMU = 90, MXPHI = 5 )
REAL  TSTFIR( MXTAU, MXCASE, MXPROB ), &
      TSTFDN( MXTAU, MXCASE, MXPROB ), &
      TSTFUP( MXTAU, MXCASE, MXPROB ), &
      TSTDFD( MXTAU, MXCASE, MXPROB ), &
      TSTUU ( MXTAU, MXMU, MXPHI, MXCASE, MXPROB )
COMMON / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU
REAL,DIMENSION(:),ALLOCATABLE     :: CMPFIR, CMPFDN, CMPFUP, CMPDFD
REAL,DIMENSION(:,:,:),ALLOCATABLE :: CMPUU 
REAL COT,MU0

real(kind=4),parameter :: PI = 2.*ASIN(1.0)
INTEGER :: NTEST, NPASS
INTEGER :: ICAS, IOD, IU, I, J, K, LC, LENTIT, LU, NPROB
CHARACTER  HEADER*127
CHARACTER  ABC(18)*1, TITLE*100, BLANKS*3
DATA  ABC / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', &
            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r' /, &
      BLANKS / '   ' / 

ACCUR = 0.0
NTEST = 0
NPASS = 0


IF( DOPROB(17) ) THEN

! c  *****************************************************
! c  ****   Test Problem 17: New delta-M-Plus
! c  *****************************************************

  USRTAU    = .TRUE.
  USRANG    = .TRUE.
  LAMBER    = .TRUE.
  PLANK     = .FALSE.
  ONLYFL    = .FALSE.
  DO_PSEUDO_SPHERE = .FALSE.
  DELTAMPLUS = .TRUE.

  NSTR = 32; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
  NLYR = 1; 
  NMOM = 32
  NTAU      = 2; IF(.NOT.USRTAU) NTAU = NLYR + 1
  NUMU      = 90; IF(.NOT.USRANG) NUMU = NSTR
  NPHI      = 3; 
  IBCND     = 0

  PRNT( 5 ) = .FALSE.

  ICAS = 2

       
  NMOM = 900
   
  call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

  CALL READ_VAL("inputFile.dat",COT,MU0)
  UMU0      = MU0;  PHI0      = 0.0
  J = 1
  DO I = 89, 0, -1
    UMU( J ) = COS(I*PI/180.0)
    J = J + 1
  ENDDO
  PHI( 1 )  = 0.0; PHI( 2 )   = 90.0; PHI( 3 )   = 180.0
  ALBEDO    = 0.0
  CALL  GETMOM( 7, 0.0, NMOM, PMOM )
  DTAUC( 1 ) = COT; SSALB( 1 ) = 1.0
  UTAU( 1 ) = 0.0; UTAU( 2 )  = DTAUC(1)
  FBEAM      = 1.0; FISOT      = 0.0
  WRITE( TITLE, '(2A)' )'Test Case No. 17', ABC(ICAS)
  LENTIT = INDEX( TITLE, BLANKS )
  HEADER  = TITLE(1:LENTIT)//                    &
                ': One Layer Kokhanovsky Cloud '


  CALL DISORT( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           &
               USRANG, USRTAU, IBCND, ONLYFL, PRNT,          &
               PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,  &          
               DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,   & 
               UTAU, UMU0, PHI0, UMU, PHI, FBEAM,            &                        
               FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,           &
               EARTH_RADIUS, H_LYR,                          &
               RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       &
               ACCUR,  HEADER,                               &
               RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          &
               ALBMED, TRNMED )          
                 
  call PRTFIN2(NTAU, NUMU, NUMU, NTAU, UTAU, UMU, UU, NPHI, PHI)                      
  
  call deallocate_disort_allocatable_arrays()

  
ENDIF
end program SINGLE_LAY_CLD

