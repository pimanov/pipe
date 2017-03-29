*
      subroutine add_nl(t,u,v,w,ox,or,ot,ut,vt,wt,Imax,Jmax)
      implicit real*8 (a-h,o-z)
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
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
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
            ut(i,j,k)=ut(i,j,k)+
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
            vt(i,j,k)=vt(i,j,k)+
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
            wt(i,j,k)=wt(i,j,k)+
     >           0.5d0*((u0*or0+u1*or1)
     >                 -(v0*ox0+v1*ox1)/(yt(j)*yt1(j)))
          end do
        end do
      end do
      return
      end
