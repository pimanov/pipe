*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,or,ot,Imax,Jmax)
*      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ut(0:Imax,0:Jmax,0:*)
     >,vt(0:Imax,0:Jmax,0:*)
     >,wt(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
* Boundary conditions
      do k=1,Km
        do j=1,Jm
          u(0,j,k)=u(Im,j,k)
          u(Im+1,j,k)=u(1,j,k)
        end do
      end do
      do i=0,Im
        do j=1,Jm
          u(i,j,0)=u(i,j,Km)
          u(i,j,Km+1)=u(i,j,1)
        end do
        do k=1,Km
          u(i,Jm+1,k)=-u(i,Jm,k)
        end do
      end do
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=0.
        end do
      end do
      do j=1,Jm
        do k=1,Km
          v(0,j,k)=v(Im,j,k)
          v(Im+1,j,k)=v(1,j,k)
        end do
        do i=1,Im
          v(i,j,0)=v(i,j,Km)
          v(i,j,Km+1)=v(i,j,1)
        end do
      end do
      do j=1,Jm
        do i=1,Im
          w(i,j,0)=w(i,j,Km)
        end do
      end do
      do k=0,Km
        do i=1,Im
          w(i,Jm+1,k)=-w(i,Jm,k)*yt(Jm)/yt(Jm+1)
        end do
        do j=1,Jm
          w(0,j,k)=w(Im,j,k)
          w(Im+1,j,k)=w(1,j,k)
        end do
      end do
*
* Vorticities
      do i=1,Im
        do j=1,Jm
          do k=0,Km
            w0=w(i,j,k)
            w1=w(i,j+1,k)
            v0=v(i,j,k)
            v1=v(i,j,k+1)
            ox(i,j,k)=((yt(j+1)*w1-yt(j)*w0)/rt1(j)
     >                -(v1-v0)/ht)/rt(j)
          end do
        end do
        j=0
          sw=0.d0
          do k=1,Km
            sw=sw+w(i,j+1,k)
          end do
          sw=sw/Km
          do k=0,Km
*            ox(i,j,k)=so
            ox(i,j,k)=2.*sw/yt(1)
          end do
      end do
      do k=0,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j,k+1)
            w0=w(i,j,k)
            w1=w(i+1,j,k)
            or(i,j,k)=(u1-u0)/(yt(j)*ht)
     >               -(w1-w0)/hx
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
            ot(i,j,k)=(v1-v0)/hx
     >               -(u1-u0)/rt1(j)
          end do
        end do
      end do
*
* Nonlinear terms
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            v0=0.5d0*(v(i,j-1,k)+v(i+1,j-1,k))
            v1=0.5d0*(v(i,j,k)+v(i+1,j,k))
            ot0=rt(j-1)*rt1(j-1)*ot(i,j-1,k)
            ot1=rt(j)*rt1(j)*ot(i,j,k)
            w0=0.5d0*(w(i,j,k-1)+w(i+1,j,k-1))
            w1=0.5d0*(w(i,j,k)+w(i+1,j,k))
            or0=or(i,j,k-1)
            or1=or(i,j,k)
            ut(i,j,k)=
     >           0.5d0*((v0*ot0+v1*ot1)/(yt(j)*yt1(j))
     >               -(w0*or0+w1*or1))
          end do
        end do
      end do
      do k=1,Km
        do i=1,Im
          do j=1,Jm-1
            w0=0.5d0*(w(i,j,k-1)+w(i,j+1,k-1))
            w1=0.5d0*(w(i,j,k)+w(i,j+1,k))
            ox0=ox(i,j,k-1)
            ox1=ox(i,j,k)
            u0=0.5d0*(u(i-1,j,k)+u(i-1,j+1,k))
            u1=0.5d0*(u(i,j,k)+u(i,j+1,k))
            ot0=ot(i-1,j,k)
            ot1=ot(i,j,k)
            vt(i,j,k)=
     >           0.5d0*((w0*ox0+w1*ox1)
     >                 -(u0*ot0+u1*ot1))
          end do
          vt(i,Jm,k)=0.d0
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u0=0.5d0*(u(i-1,j,k)+u(i-1,j,k+1))
            u1=0.5d0*(u(i,j,k)+u(i,j,k+1))
            or0=or(i-1,j,k)
            or1=or(i,j,k)
            v0=0.5d0*(v(i,j-1,k)+v(i,j-1,k+1))
            v1=0.5d0*(v(i,j,k)+v(i,j,k+1))
            ox0=rt(j-1)*rt1(j-1)*ox(i,j-1,k)
            ox1=rt(j)*rt1(j)*ox(i,j,k)
            wt(i,j,k)=
     >           0.5d0*((u0*or0+u1*or1)
     >                 -(v0*ox0+v1*ox1)/(yt(j)*yt1(j)))
          end do
        end do
      end do
*
* Viscous terms
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            ot0=ot(i,j-1,k)
            ot1=ot(i,j,k)
            or0=or(i,j,k-1)
            or1=or(i,j,k)
            ut(i,j,k)=ut(i,j,k)
     >               -((rt(j)*ot1-rt(j-1)*ot0)/yt1(j)
     >               -(or1-or0)/ht)/yt(j)/Re
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm-1
          do i=1,Im
            ox0=ox(i,j,k-1)
            ox1=ox(i,j,k)
            ot0=ot(i-1,j,k)
            ot1=ot(i,j,k)
            vt(i,j,k)=vt(i,j,k)
     >               -((ox1-ox0)/(rt(j)*ht)
     >                -(ot1-ot0)/hx)/Re
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            ox0=ox(i,j-1,k)
            ox1=ox(i,j,k)
            or0=or(i-1,j,k)
            or1=or(i,j,k)
            wt(i,j,k)=wt(i,j,k)
     >               -((or1-or0)/hx
     >                -(ox1-ox0)/yt1(j))/Re
          end do
        end do
      end do
      return
      end