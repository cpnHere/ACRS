       subroutine write_data( file_head,&
                              Dmax,wl,mr,mi,&
                              Qe, al, gf,&
                              NUMANG,ANGLE,&
                              P11,P12,P33,P34)
       implicit none
       character*127 :: file_head
       real    :: Dmax ,wl, mr, mi,Qe, al, gf       
       integer :: NUMANG
       real    :: ANGLE(NUMANG)
       real    :: P11(NUMANG),P12(NUMANG),P33(NUMANG),P34(NUMANG)
       
       integer :: i,j,k
       integer :: LA1,LA2,LB1,LB2
       character*30 :: F_name_size,F_name_wl,Filename
       print*,'input Dmax',Dmax 
       WRITE(F_NAME_WL,"(F7.3)") wl
       WRITE(F_NAME_SIZE,"(F10.3)") Dmax

       DO I = 1, 30
       LA1 = I
       IF(F_NAME_WL(I:I) /= " ") EXIT
       END DO

       DO I = 30, 1, -1
       LA2 = I
       IF(F_NAME_WL(I:I) /= " ") EXIT
       END DO

       DO I = LA1, LA2
       IF(F_NAME_WL(I:I) == ".") F_NAME_WL(I:I)="p"
       END DO


       DO I = 1, 30
       LB1 = I
       IF(F_NAME_SIZE(I:I) /= " ") EXIT
       END DO


       DO I = 30, 1, -1
       LB2 = I
       IF(F_NAME_SIZE(I:I) /= " ") EXIT
       END DO

       DO I = LB1, LB2
       IF(F_NAME_SIZE(I:I) == ".") F_NAME_SIZE(I:I)="p"
       END DO

       OPEN(UNIT=1,FILE=trim(file_head)//"_WL_"//F_NAME_WL(LA1:LA2) &
       //"_SIZE_"//F_NAME_SIZE(LB1:LB2),STATUS="REPLACE")

       WRITE(1,*)'Diameter (micron) of the sphere  D=  '
       write(1,'(1x,E13.7)')Dmax
  
       write(1,*)'WL,mr,mi :'
       WRITE(1, '((3(1x,E12.4)))') wl, MR, MI
       
       WRITE(1,7103 )Qe, al, gf
 7103  FORMAT(1X,'  Qe,    albedo,     g-factor:',/1X, 4(E13.6, 3X))


      WRITE(1,7104)
 7104 FORMAT(10X/1X, ' scattering angle , phase function  '/1X )
      write(1,7106)
7106  FORMAT(1X, ' ang      P11     P12/P11   P33/P11  P34/P11   '/1X )
      DO k = 1, Numang
         WRITE(1,7105)angle(k),p11(k),p12(k), p33(k),p34(k)
      END DO

7105  format(1x,f6.2,1x,6(E13.7,2x))
      close(1)
      end subroutine write_data
       
