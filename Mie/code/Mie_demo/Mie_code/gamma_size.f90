       subroutine gamma_ditribution(Nbin,Deff,Veff,rad,fn)
          implicit none
          !------I/O--------!
          integer :: Nbin
          real :: Deff,Veff
          real :: rad(Nbin),fn(Nbin)
          
          !----local variables-----!
          real :: Reff
          real :: rad1,rad2,delrad,rlow,rhigh
          integer :: i,j,k,ilow,ihigh
          real,parameter  :: eps =1e-5,large_num=1e30 
		  logical,parameter::dump_fn = .true.
          real,parameter  :: r_max =1000.0
          ! Because the numerical computation of a 
          ! size distribution must be trucated some where
          ! there is always some area of the size distribution
          ! that will be ignored
          ! the eps here means the accuracy of the trucation
          ! the area ignored by the trucation is smaller that eps

          reff = deff/2.0        
          rad1 = 0.0;
          !rad2 = min(600.0, exp(log(large_num)/( (1.0-3.0*veff)/veff) ))
          rad2=r_max
           delrad = (rad2-rad1) / float(nbin)
          
          do k = 1, nbin
             rad(k) = rad1 + (float(k) -0.5) * delrad
             fn(k) = exp(-rad(k) /(reff * veff) )    &
            *(rad(k) ** ((1.0-3.0*veff)/veff))             
          enddo
           
          ! compute the cumulative distribution
          do k = 2, nbin
             fn(k) = fn(k) + fn(k-1) 
          enddo
         
          ! normalized the cumulative distribution
          do k = 1, Nbin
             fn(k) = fn(k) / fn(Nbin)
          enddo
            
          ! search for the lower end of the size
          ! note that the cumulative distribution 
          ! for size < rlow is < 0.001
          ilow = 1
          ihigh = Nbin
R1:     do while(ihigh - ilow .gt. 1) 
              k = (ilow+ihigh) / 2
              if(fn(k) < eps) then
                ilow = k
              else
                ihigh = k
              end if
          end do R1   
 
          rlow = (rad(ilow) + rad(ihigh) ) / 2.0
      
          ilow = 1
          ihigh = Nbin 
r2:    do while(ihigh - ilow .gt. 1) 
             k = (ilow+ihigh) / 2
             if(fn(k) < 1-eps) then
                ilow = k
             else
                ihigh = k
             end if
          end do r2
          
          rhigh = (rad(ilow) + rad(ihigh) ) / 2.0
          delrad = (rhigh-rlow) / float(Nbin)
            
          do k = 1, Nbin               
               rad(k) = rlow + (float(k) -0.5) * delrad
               fn(k) = exp(-rad(k) /(reff * veff) ) &
              *(rad(k) ** ((1.0-3.0*veff)/veff))
          enddo   
          
          fn = fn/maxval(fn)
          if (dump_fn) then
			open(unit=1,file='fndump.dat',status='unknown',action='write')
			do i=1,Nbin
			 write(1,*)rad(i),fn(i)
			enddo
			close(1)
          endif
          return
          end subroutine gamma_ditribution
