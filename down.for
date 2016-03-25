*
      subroutine down(fncp,u,v,w,t,dt,Dp,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
      open(9,file=fncp,form='unformatted')
      write(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
      do k=1,Km
        do j=1,Jm
          write(9)(u(i,j,k),i=1,Im)
          write(9)(v(i,j,k),i=1,Im)
          write(9)(w(i,j,k),i=1,Im)
        end do
      end do
      close(9)
      return
      end
*
