*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Jmax=65)
      character*12 fnscs
      dimension
     > u(0:Jmax)
     >,v(0:Jmax)
      common
     >/dim/epsr,nsym
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
*
      open(5,file='Pois_init.car')
      read(5,*) Jm,epsr
      read(5,*) lt
      read(5,*) fnscs
      nsym=1
      call com
      write(*,*)' Jm=',Jm,',epsr=',epsr
      write(*,*)' Km=',Km,',nsym doesn`t matter'
*
      do j=1,Jm
        u(j)=1.0-yt(j)**2
        v(j)=0.0
      end do
      open(9,file=fnscs,form='unformatted')
      write(9)epsr,Jm,lt,nsym
      do k=1,Km
        write(9)(u(j),j=1,Jm)
        write(9)(v(j),j=1,Jm)
        write(9)(v(j),j=1,Jm)
      end do
      write(9) 'Pois_init.out'
      close(9)
      write(*,*)' Poisel flow field saved in ',fnscs
*
      stop
      end
