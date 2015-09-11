*
      subroutine divmax(u,v,w,dd,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
*
* Boundary conditions
      do k=1,Km
        do j=1,Jm
          u(0,j,k)=u(Im,j,k)
        end do
      end do
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=0.d0
        end do
      end do
      do j=1,Jm
        do i=1,Im
          w(i,j,0)=0.d0
          w(i,j,Km)=0.d0
        end do
      end do
*
      dd=-1.d0
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      return
      end
