! *** important *** !
! To use this program properly, the input refrative index
! must be in the form mr + i*mi where mr >1 and mi >0
       program mie_single_size
       use netcdf
       implicit none

       integer,parameter :: MX_Nband  = 1000
       integer,parameter :: MX_Nsize  = 50000
       integer,parameter :: MX_Numang = 18001

!----------variables for MIEV0-----------------!
       real    :: XX
       complex :: CREFIN
       logical :: PERFCT
       real    :: MIMCUT
       logical :: ANYANG
!       real    :: XMU( MAXANG )
       integer :: NMOM
       integer :: IPOLZN
       integer,parameter :: MOMDIM = 200
       logical :: PRNT(2)
       real    :: QEXT, QSCA, GQSC
       real    :: PMOM( 0:MOMDIM, 4 )
       complex :: SFORW,SBACK!, S1( MAXang ), S2( MAXang )
       complex :: TFORW, TBACK
       real    :: SPIKE
       character*4 :: schar,bchar
!-----------------------------------------------!

!-----------local variables---------------------!
       real    :: PI
       integer :: Nband, Nsize, Numang
       real,allocatable  :: wl(:), mr(:), mi(:)
       real,allocatable  :: Dmax(:)
       real,allocatable  :: ANGLE(:), xmu(:)
       real,allocatable  :: P11_all(:,:,:),P12_all(:,:,:),P33_all(:,:,:),P34_all(:,:,:)
       real,allocatable  :: P11(:),P12(:),P33(:),P34(:)
       complex,allocatable::S1(:),S2(:)
       real,allocatable  :: Qe_all(:,:), al_all(:,:), gf_all(:,:)
       real  :: Qe, al, gf

       real    :: i1,i2
       real    :: sum_P11
       integer :: i,j,k
       character*127 :: header
       complex :: img
       real    :: S11,S12,S33,S34
       character*200 :: output_filename
       character*3 :: bandsub_char
       character*5 :: size_char
!--------- variables for netcdf -----------------!
       integer :: ncStatus, ncFileId, sizeDimId, bandDimId,angDimId,&
                  dirDimId, ncVarId, &
                  numX, numY, numZ, numIntensityDirs

!------------------------------------------------!
!******** code begin********!
       PI = 2.0*asin(1.0)
       !-------read the wavelength and refractive index------!
       open(unit=1,file='wl_ref.dat',status='old',action='read')
       read(1,*)Nband
       allocate(wl(nband),mr(nband),mi(nband))
       if (Nband .gt. MX_Nband) then
           print*,'too many wavelength, please increase MX_Nband'
           stop
       endif
       do i=1, Nband
           read(1,*)wl(i),mr(i),mi(i)
       end do
       close(1)
       print*,'finish reading wavelength information'
       !------read the particle sizes--------!
       open(unit=2,file='size.dat',status='old',action='read')
       read(2,*)Nsize
       allocate(Dmax(Nsize))
       if (Nsize .gt. MX_Nsize) then
           print*,'too many sizes, please increase MX_Nsize'
           stop
       endif
       do i=1, Nsize
           read(2,*)Dmax(i)
       end do
       close(2)
       !print*,'finish reading size information'
       !------- read in scattering angles-----!
       open(unit=3,file='ang.dat',status='old',action='read')
       read(3,*)Numang
       allocate(angle(Numang),xmu(Numang))
       if (Numang .gt. MX_Numang) then
           print*,'too many angle, please increase MX_Numang'
           stop
       endif
       do i=1, Numang
           read(3,*)angle(i)
           xmu(i)=cos(angle(i)*PI/180.0)
       end do
       close(3)
       print*,'finish reading angle information'
!------ allocate variables---------!
       allocate(Qe_all(Nband,Nsize),al_all(Nband,Nsize),gf_all(Nband,Nsize))
       allocate(P11_all(Numang,Nband,Nsize),P12_all(Numang,Nband,Nsize))
       allocate(P33_all(Numang,Nband,Nsize),P34_all(Numang,Nband,Nsize))
       allocate(P11(Numang),P12(Numang),P33(Numang),P34(Numang),S1(Numang),S2(Numang))

       !---- read the output file name ---
       open(unit=3,file='filename.dat',status='old',action='read')
       read(3,'(A200)')output_filename
	  close(3)
!------- begin to compute the Mie scattering properties ------------!
       !------initialize some MIEV0 variables----------!
       PERFCT = .False.
       MIMCUT = 1.E-12
       ANYANG = .True.
       NMOM   = 0
       IPOLZN = 1234
       PRNT( 1 ) = .False.
       PRNT( 2 ) = .False.

       write(header,'(A)')'Mie_mono'

size:     do j = 1, Nsize
             print*,'processing size',j
band:      do i = 1, Nband
             CREFIN = CMPLX(mr(i),  mi(i))
             img=(0,1)
             XX = PI * Dmax(j) / wl(i) ! size parameter
             P11 = 0.0
             P12 = 0.0
             P34 = 0.0
             P33 = 0.0
             S1 = 0.0
             S2 = 0.0
             !print*,XX, CREFIN, PERFCT, MIMCUT
             !print*,ANYANG, NUMANG, NMOM, IPOLZN,MOMDIM
              CALL MIEV0( XX, CREFIN, PERFCT, MIMCUT, &
                   ANYANG, NUMANG, XMU, NMOM, &
                   IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC, PMOM, SFORW,&
                   SBACK, S1, S2, TFORW, TBACK, SPIKE )
              !print*,'MIE computation finished'

              Qe = QEXT
              al = QSCA / QEXT
              gf = GQSC / QSCA
angle1:       do k=1,NUMANG
                !the equations below follow those in van de Hulst's book Page 44
                ! note that the definition of (I Q U V) of this mie code is
                ! I=E_lxE_l* + E_rxE_r*
                ! Q=E_lxE_l* - E_rxE_r*
                ! U=E_lxE_r* + E_rxE_l*
                ! V=i(E_lxE_r* - E_rxE_l*)
                ! i.e right handness

                S11=0.5*(ABS( S1(k) )**2 + ABS( S2(k) )**2)
                S12=0.5*(ABS( S2(k) )**2 - ABS( S1(k) )**2)
                S33=0.5*(CONJG(S2(k))*S1(k) + S2(k)*CONJG(S1(k)))
                S34=0.5*img*(S1(k)*CONJG(S2(k)) - S2(k)*CONJG(S1(k)))
                !-----------------------------------------------------!
                P11(k)=S11
                P12(k)=S12
                P33(k)=S33
                P34(k)=S34
              end do angle1
              P12(1:NUMANG) =  P12(1:NUMANG)/P11(1:NUMANG)
              P33(1:NUMANG) =  P33(1:NUMANG)/P11(1:NUMANG)
              P34(1:NUMANG) =  P34(1:NUMANG)/P11(1:NUMANG)
              sum_P11 = sum( 0.5*( P11(1:NUMANG-1) + P11(2:NUMANG) )*&
                           ( XMU(1:NUMANG-1)- XMU(2:NUMANG) ) )/2.0
              P11 = P11/sum_P11
              Qe_all(i,j)=Qe
              al_all(i,j)=al
              gf_all(i,j)=gf
              P11_all(:,i,j)=P11
              P12_all(:,i,j)=P12*P11
              P33_all(:,i,j)=P33*P11
              P34_all(:,i,j)=P34*P11
              end do band
         end do size
       !---- open a netcdf file to store the output-----!
       ncStatus = nf90_create(trim(output_FileName), nf90_clobber, ncFileId)
       ncStatus = nf90_put_att(ncFileId, NF90_Global, "description", &
                        "Output from Mie scattering computation")
       ncStatus = nf90_def_dim(ncFileId, "Diameter", Nsize, SizeDimId)
       ncStatus = nf90_def_dim(ncFileId, "Wavelength", Nband, bandDimId)
       ncStatus = nf90_def_dim(ncFileId, "PhaseFunctionAngle",Numang,angDimId)
       ncStatus = nf90_def_var(ncFileId, "Diameter", &
                                nf90_float, sizeDimId, ncVarId)
       ncStatus = nf90_def_var(ncFileId, "Wavelength", &
                                nf90_float, bandDimId, ncVarId)
       ncStatus = nf90_def_var(ncFileId, "Refr_real", &
                                nf90_float, bandDimId, ncVarId)
       ncStatus = nf90_def_var(ncFileId, "Refr_img", &
                                nf90_float, bandDimId, ncVarId)
       ncStatus = nf90_def_var(ncFileId, "PhaseFunctionAngle", &
                                nf90_float, angDimId, ncVarId)
       ncStatus = nf90_def_var(ncFileId, "SingleScatteringAlbedo", &
                                nf90_float, (/ bandDimId, sizeDimId /), ncVarId)
       ncStatus = nf90_def_var(ncFileId, "ExtinctionEfficiency", &
                                nf90_float, (/ bandDimId, sizeDimId /), ncVarId)
       ncStatus = nf90_def_var(ncFileId, "AsymmetryFactor", &
                                nf90_float, (/ bandDimId, sizeDimId /), ncVarId)
        ncStatus = nf90_def_var(ncFileId, "P11", &
                                nf90_float, (/ angDimId,bandDimId, sizeDimId /), ncVarId)
        ncStatus = nf90_def_var(ncFileId, "P12", &
                                nf90_float, (/ angDimId,bandDimId, sizeDimId /), ncVarId)
        ncStatus = nf90_def_var(ncFileId, "P33", &
                                nf90_float, (/ angDimId,bandDimId, sizeDimId /), ncVarId)
        ncStatus = nf90_def_var(ncFileId, "P34", &
                                nf90_float, (/ angDimId,bandDimId, sizeDimId /), ncVarId)
        ncStatus = nf90_EndDef(ncFileId)
        !---- file set up now write each variables-------!
         ncStatus = nf90_inq_varid(ncFileId, "Diameter", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, Dmax(1:Nsize))
         ncStatus = nf90_put_att(ncFileId, ncVarId, 'Units','Micron')

         ncStatus = nf90_inq_varid(ncFileId, "Wavelength", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, wl(:))
         ncStatus = nf90_put_att(ncFileId, ncVarId, 'Units','Micron')

         ncStatus = nf90_inq_varid(ncFileId, "Refr_real", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, mr(:))
         ncStatus = nf90_put_att(ncFileId, ncVarId, 'long_name','Real part of refractive index')

         ncStatus = nf90_inq_varid(ncFileId, "Refr_img", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, mi(:))
         ncStatus = nf90_put_att(ncFileId, ncVarId, 'long_name','Imaginary  part of refractive index')

         ncStatus = nf90_inq_varid(ncFileId, "PhaseFunctionAngle", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, angle(:))
         ncStatus = nf90_put_att(ncFileId, ncVarId, 'units','Degree')

          ncStatus = nf90_inq_varid(ncFileId, "SingleScatteringAlbedo", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, al_all(:,:))
         ncStatus = nf90_inq_varid(ncFileId, "ExtinctionEfficiency", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, Qe_all(:,:))
         ncStatus = nf90_inq_varid(ncFileId, "AsymmetryFactor", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, gf_all(:,:))
         ncStatus = nf90_inq_varid(ncFileId, "P11", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, P11_all(:,:,:))
         ncStatus = nf90_put_att(ncFileId, ncVarId, 'long_name', 'P11 is normailzed')
         ncStatus = nf90_inq_varid(ncFileId, "P12", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, P12_all(:,:,:))
         ncStatus = nf90_inq_varid(ncFileId, "P33", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, P33_all(:,:,:))
         ncStatus = nf90_inq_varid(ncFileId, "P34", ncVarId)
         ncStatus = nf90_put_var(ncFileId, ncVarId, P34_all(:,:,:))
         ncStatus = nf90_close(ncFileId)
         deallocate(Dmax,wl,mr,mi,Qe_all,al_all,gf_All,angle,xmu)
         deallocate(P11_all,P12_all,P33_all,P34_All)
         deallocate(P11,P12,P33,P34)
         end program mie_single_size
