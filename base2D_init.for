*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Imax=1025)
      character*12 fncp,fnbcp
      character*256 comment
      dimension buf(2*Imax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
     >/cf/cf
*
      open(5,file='base2D_init.car')
      read(5,*) cf
      read(5,*) ics
      read(5,*) Xmax_new, lx_new
      read(5,*) fncp
      read(5,*) fnbcp
      read(5,*, err=500) comment
500   close(5)
      write(*,*) 'base_init.car: cf=',cf,' ics=',ics
     >  ,' fncp=',fncp,' fnbcp=',fnbcp,' Xmax=',Xmax_new,' lx=',lx_new
      write(*,*) 'comment:',comment
      open(9,file=fncp,form='unformatted',status='old',err=111)
      open(7,file=fnbcp,form='unformatted')
      read(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
      write(7)cf,Xmax_new,epsr,lx_new,Jm,lt,nsym
      Im_new=2**lx_new
*
      Npm=1
      call com
      write(*,*)' ***************************************************'
      write(*,*)' *        Number of processors =',Npm,'          *'
      write(*,200) t,dt,Dp,Re,Xmax,epsr,Im,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    t=',1pe10.3,' dt=',e9.2,' Dp=',e9.2,/,
     >'    Re=',e9.2,/,
     >'    Xmax=',e9.2,/,
     >'    epsr=',e9.2,' Im=',i4,' Jm=',i4,' Km=',i4,' Nsym=',i3)
*
      do k=1,Km
        do j=1,Jm
          read(9) (buf(i),i=1,Im)
          ut=buf(ics)
          do i=1,Im_new
            buf(i)=ut
          end do
          write(7) (buf(i),i=1,Im_new)

          read(9) (buf(i),i=1,Im)
          vt=buf(ics)
          do i=1,Im_new
            buf(i)=vt
          end do
          write(7) (buf(i),i=1,Im_new)

          read(9) (buf(i),i=1,Im)
          wt=buf(ics)
          do i=1,Im_new
            buf(i)=wt
          end do
          write(7) (buf(i),i=1,Im_new)
        end do
      end do
      close(9)
      write(7) comment
      close(7)


*
      stop
111   write(*,*)'  cp file ',fncp,' was not found'
      stop
      end
