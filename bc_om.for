*
      subroutine bc_om(u,v,w,ox,or,ot,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
*
* Boundary conditions
      do k=1,Km
        do j=1,Jm
          u(Im+1,j,k)=u(1,j,k)
          u(0,j,k)=u(Im,j,k)
          v(Im+1,j,k)=v(1,j,k)
          v(0,j,k)=v(Im,j,k)
          w(Im+1,j,k)=w(1,j,k)
          w(0,j,k)=w(Im,j,k)
        end do
      end do
*
      do i=0,Im
        do j=1,Jm
          u(i,j,0)=u(i,j,1)
          u(i,j,Km+1)=u(i,j,Km)
        end do
        do k=1,Km
          u(i,0,k)=u(i,1,k)
          u(i,Jm+1,k)=-u(i,Jm,k)
        end do
      end do
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=0.d0
          v(i,0,k)=0.d0
        end do
      end do
      do j=1,Jm
        do i=1,Im
          v(i,j,0)=v(i,j,1)
          v(i,j,Km+1)=v(i,j,Km)
        end do
      end do
      do j=1,Jm
        do i=0,Im+1
          w(i,j,0)=0.d0
          w(i,j,Km)=0.d0
        end do
      end do
      do k=0,Km
        do i=1,Im
          w(i,0,k)=w(i,1,k)
          w(i,Jm+1,k)=-w(i,Jm,k)
        end do
      end do
*
* Vorticities
      do i=1,Im
        do j=0,Jm
          do k=0,Km
            w0=w(i,j,k)
            w1=w(i,j+1,k)
            v0=v(i,j,k)
            v1=v(i,j,k+1)
            ox(i,j,k)=(w1-w0)/rt1(j)-(v1-v0)/ht
          end do
        end do
      end do
      do k=0,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j,k+1)
            w0=w(i,j,k)
            w1=w(i+1,j,k)
            or(i,j,k)=(u1-u0)/ht-(w1-w0)/hx
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j+1,k)
            v0=v(i,j,k)
            v1=v(i+1,j,k)
            ot(i,j,k)=(v1-v0)/hx-(u1-u0)/rt1(j)
          end do
        end do
      end do
*
      return
      end
