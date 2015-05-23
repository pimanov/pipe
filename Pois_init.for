*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Imax=2049)
      character*12 fnbcp
      character*256 comment
      dimension
     >buf(Imax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
     >/cf/cf
*
      open(5,file='Pois_init.car')
      read(5,*) cf
      read(5,*) Xmax, lx
      read(5,*) Jm, epsr
      read(5,*) nsym, lt
      read(5,*) fnbcp
      read(5,*, err=500) comment
500   write(*,*) 'base_init.car: cf=',cf,' fncp=',fncp,' fnbcp=',fnbcp
      write(*,*) 'comment:',comment
*
      Npm=1
      call com
      write(*,*)' ***************************************************'
      write(*,*)' *        Number of processors =',Npm,'          *'
      write(*,200) cf,Xmax,epsr,Imm,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    cf=',1pf15.10,'    Xmax=',e9.2,/,
     >'    epsr=',e9.2,' Im=',i4,' Jm=',i4,' Km=',i4,' Nsym=',i3)
*
      open(9,file=fnbcp,form='unformatted')
      write(9)cf,Xmax,epsr,lx,Jm,lt,nsym
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            buf(i)=1.0-yt(j)**2
          end do
          write(9)(buf(i),i=1,Imm)
          do i=1,Im
            buf(i)=0.0
          end do
          write(9)(buf(i),i=1,Imm)
          write(9)(buf(i),i=1,Imm)
        end do
      end do
      write(9) comment
      close(9)
*
      stop
      end
