	program driver
	implicit none
	
	real :: band(2) 
    integer :: Nwlpbd
	real :: Deff_st,Deff_ed,veff
	complex :: refr
	complex,external :: REFWAT
	integer,parameter :: Nbin = 10000!
	integer :: Nsize    !=10000
    real:: rad(Nbin),fn(Nbin)
	character*3::bin_type
	integer :: i,j,k
	real,parameter :: T=290.0     
	real::wl,tmp1,tmp2
	real::size_st,size_ed,Deff
	character*100::output_filename
      integer::nang
      real,allocatable :: ang(:)
	open(unit=1,file='driver_input.dat',status='old',action='read')
	read(1,*)band(1:2),Nwlpbd
	read(1,*)Deff_st,Deff_ed,veff,Nsize
	read(1,*)bin_type
	read(1,*)nang
      read(1,*)output_filename
    close(1)
    open(unit=1,file='wl_ref.dat',status='unknown',action='write')
    write(1,*)Nwlpbd
    do i=1,Nwlpbd
		wl = band(1)+(band(2)-band(1))/real(Nwlpbd) * real(i)
		refr = refwat(wl,T)
		write(1,*)wl,real(refr),imag(refr)
	enddo
	close(1)
    if (Deff_st .lt. Deff_ed) then ! gamma or log-normal distributioni
     print*,Deff_st,Deff_ed
    !call gamma_ditribution(Nbin,Deff_st,Veff,rad,fn)
!	size_st = rad(1)
 !   call gamma_ditribution(Nbin,Deff_ed,Veff,rad,fn)
!	size_ed = rad(Nbin)
    open(unit=1,file='size.dat',status='unknown',action='write')
      write(1,*) Nsize
    if (bin_type .eq. 'LOG') then
            tmp1=log10(Deff_st)
            tmp2=log10(Deff_ed)
            do i=1,Nsize
                  Deff = 10.0**(tmp1+(tmp2-tmp1)/real(Nsize) * real(i))*2.0
                  write(1,*)Deff
            enddo
    elseif(bin_type .eq. 'BLG') then
            tmp1=log10(Deff_st)
            tmp2=log10(Deff_ed)
            do i=1,Nsize
                  Deff = 10.0**(tmp1+(tmp2-tmp1)/real(Nsize) * real(i))
                  write(1,*)Deff
            enddo

    elseif (bin_type .eq. 'SIG') then
      write(1,*)Deff_st
    elseif (bin_type .eq. 'BIN') then
        tmp1=(Deff_st)
        tmp2=(Deff_ed)
        do i=1,Nsize
            Deff = (tmp1+(tmp2-tmp1)/real(Nsize) * real(i))
            write(1,*)Deff
        enddo
     else
        tmp1=(size_st)
            tmp2=(size_ed)
            do i=1,Nsize
                  Deff = (tmp1+(tmp2-tmp1)/real(Nsize) * real(i))*2.0
                  write(1,*)Deff
            enddo
    endif
      close(1)

!	elseif (Deff_st .eq. Deff_ed) then ! single size computation
!	 size_st = Deff_st
!	 size_ed = Deff_ed
!	 open(unit=1,file='size.dat',status='unknown',action='write')
!	 write(1,*)'1'
!	 write(1,*)Deff_st
!	 close(1)
	else
	 print*,'input Deff_st > Deff_ed...stopping'
	 stop
	endif
!    open(unit=1,file='size.dat',status='unknown',action='write')
!	write(1,*) Nsize
!    if (bin_type .eq. 'LOG') then
!		tmp1=log10(size_st)
!		tmp2=log10(size_ed)
!		do i=1,Nsize
!			Deff = 10.0**(tmp1+(tmp2-tmp1)/real(Nsize) * real(i))*2.0
!			write(1,*)Deff
!		enddo
!    elseif (bin_type .eq. 'SIG') then
!	write(1,*)Deff_st
!	else
!        tmp1=(size_st)
!		tmp2=(size_ed)
!		do i=1,Nsize
!			Deff = (tmp1+(tmp2-tmp1)/real(Nsize) * real(i))*2.0
!			write(1,*)Deff
!		enddo
!   endif
!	close(1)

	open(unit=1,file='filename.dat',status='unknown',action='write')
	write(1,'(A)')trim(output_filename)
	close(1)
      open(unit=1,file='ang.dat',status='unknown',action='write')
      write(1,*)nang
     do i=1,nang
      write(1,'(f9.3)') real( I - 1 )*180. / real( Nang - 1 )
       end do 

	end program driver	
