*
      subroutine visc(ox,or,ot,ut,vt,wt,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
     >,ut(0:Imax,0:Jmax,0:*)
     >,vt(0:Imax,0:Jmax,0:*)
     >,wt(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
* Viscous terms
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            ot0=ot(i,j-1,k)
            ot1=ot(i,j,k)
            or0=or(i,j,k-1)
            or1=or(i,j,k)
            ut(i,j,k)=-((ot1-ot0)/yt1(j)-(or1-or0)/ht)/Re
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
            vt(i,j,k)=-((ox1-ox0)/ht-(ot1-ot0)/hx)/Re
          end do
        end do
      end do
      do k=1,Km-1
        do j=1,Jm
          do i=1,Im
            ox0=ox(i,j-1,k)
            ox1=ox(i,j,k)
            or0=or(i-1,j,k)
            or1=or(i,j,k)
            wt(i,j,k)=-((or1-or0)/hx-(ox1-ox0)/yt1(j))/Re
          end do
        end do
      end do
*
      return
      end
