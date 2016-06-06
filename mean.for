*
      subroutine mean(u,a,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
*
      a = 0.d0
      v = 0.d0
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            v1=hx*yt1(j)*ht*yt(j)
            a=a+u(i,j,k)*v1
            v=v+v1
          end do
        end do
      end do
      a=a/v
      return
      end
