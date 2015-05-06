* 
      subroutine rp_base(t,u,v,w,ox,or,ot,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,ox(0:Jmax,0:*)
     >,or(0:Jmax,0:*)
     >,ot(0:Jmax,0:*)
      common
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
*
* Boundary conditions
      do j=1,Jm
        u(j,0)=u(j,1)
        u(j,Km+1)=u(j,Km)
        v(j,0)=v(j,1)
        v(j,Km+1)=v(j,Km)
        w(j,0)=0.d0
        w(j,Km)=0.d0
      end do
      do k=0,Km
        u(Jm+1,k)=-u(Jm,k)
        v(Jm,k)=0.d0
        w(Jm+1,k)=-w(Jm,k)*yt(Jm)/yt(Jm+1)
      end do
*
* Vorticities
      do j=1,Jm
        do k=0,Km
          w0=w(j,k)
          w1=w(j+1,k)
          v0=v(j,k)
          v1=v(j,k+1)
          ox(j,k)=((yt(j+1)*w1-yt(j)*w0)/rt1(j)-(v1-v0)/ht)/rt(j)
        end do
      end do
      do k=0,Km
        ox(0,k)=0.d0
      end do
      do k=0,Km
        do j=1,Jm
          u0=u(j,k)
          u1=u(j,k+1)
          or(j,k)=(u1-u0)/(yt(j)*ht)
        end do
      end do
      do k=1,Km
        do j=1,Jm
          u0=u(j,k)
          u1=u(j+1,k)
          ot(j,k)=-(u1-u0)/rt1(j)
        end do
      end do
*
      return
      end
