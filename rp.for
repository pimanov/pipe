* 
      subroutine rp(t,ub,vb,wb,obx,obr,obt,u,v,w,ut,vt,wt,ox,or,ot,Jmax)
      implicit real*8 (a-h,o-z)
      complex*16 u,v,w,ut,vt,wt,ox,or,ot,c0,ci,cu0,cu1,cv0,cv1,cw0,cw1,
     > cv,cw,cox0,cox1,cor0,cor1,cot0,cot1,cor,cot
      dimension
     > ub(0:Jmax,0:*)
     >,vb(0:Jmax,0:*)
     >,wb(0:Jmax,0:*)
     >,obx(0:Jmax,0:*)
     >,obr(0:Jmax,0:*)
     >,obt(0:Jmax,0:*)
     >,u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,ut(0:Jmax,0:*)
     >,vt(0:Jmax,0:*)
     >,wt(0:Jmax,0:*)
     >,ox(0:Jmax,0:*)
     >,or(0:Jmax,0:*)
     >,ot(0:Jmax,0:*)
      common
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/alpha/alpha
*
      c0=(0.d0,0.d0)
      ci=(0.d0,1.d0)
*
* Boundary conditions
      do j=1,Jm
        u(j,0)=u(j,1)
        u(j,Km+1)=u(j,Km)
        v(j,0)=v(j,1)
        v(j,Km+1)=v(j,Km)
        w(j,0)=c0
        w(j,Km)=c0
      end do
      do k=0,Km
        u(Jm+1,k)=-u(Jm,k)
        v(Jm,k)=c0
        w(Jm+1,k)=-w(Jm,k)*yt(Jm)/yt(Jm+1)
      end do
*
* Vorticities
      do j=1,Jm
        do k=0,Km
          cw0=w(j,k)
          cw1=w(j+1,k)
          cv0=v(j,k)
          cv1=v(j,k+1)
          ox(j,k)=((yt(j+1)*cw1-yt(j)*cw0)/rt1(j)-(cv1-cv0)/ht)/rt(j)
        end do
      end do
      do k=0,Km
        ox(0,k)=c0
      end do
      do k=0,Km
        do j=1,Jm
          cu0=u(j,k)
          cu1=u(j,k+1)
          cw=w(j,k)
          or(j,k)=(cu1-cu0)/(yt(j)*ht)-ci*alpha*cw
        end do
      end do
      do k=1,Km
        do j=1,Jm
          cu0=u(j,k)
          cu1=u(j+1,k)
          cv=v(j,k)
          ot(j,k)=ci*alpha*cv-(cu1-cu0)/rt1(j)
        end do
      end do
*
* Nonlinear terms [om x U]
      do j=1,Jm
        do k=1,Km
          v0=vb(j-1,k)
          v1=vb(j,k)
          cot0=rt(j-1)*rt1(j-1)*ot(j-1,k)
          cot1=rt(j)*rt1(j)*ot(j,k)
          w0=wb(j,k-1)
          w1=wb(j,k)
          cor0=or(j,k-1)
          cor1=or(j,k)
          ut(j,k)=
     >           0.5d0*((v0*cot0+v1*cot1)/(yt(j)*yt1(j))
     >               -(w0*cor0+w1*cor1))
        end do
      end do
      do k=1,Km
        do j=1,Jm-1
          w0=0.5d0*(wb(j,k-1)+wb(j+1,k-1))
          w1=0.5d0*(wb(j,k)+wb(j+1,k))
          cox0=ox(j,k-1)
          cox1=ox(j,k)
          u0=0.5d0*(ub(j,k)+ub(j+1,k))
          cot0=ot(j,k)
          vt(j,k)=0.5d0*(w0*cox0+w1*cox1)-u0*cot0
        end do
        vt(Jm,k)=c0
      end do
      do k=1,Km
        do j=1,Jm
          u0=0.5d0*(ub(j,k)+ub(j,k+1))
          cor0=or(j,k)
          v0=0.5d0*(vb(j-1,k)+vb(j-1,k+1))
          v1=0.5d0*(vb(j,k)+vb(j,k+1))
          cox0=rt(j-1)*rt1(j-1)*ox(j-1,k)
          cox1=rt(j)*rt1(j)*ox(j,k)
          wt(j,k)=u0*cor0-0.5d0*(v0*cox0+v1*cox1)/(yt(j)*yt1(j))
        end do
      end do
*
* Nonlinear terms [OM x u]
      do j=1,Jm
        do k=1,Km
          cv0=v(j-1,k)
          cv1=v(j,k)
          ot0=rt(j-1)*rt1(j-1)*obt(j-1,k)
          ot1=rt(j)*rt1(j)*obt(j,k)
          cw0=w(j,k-1)
          cw1=w(j,k)
          or0=obr(j,k-1)
          or1=obr(j,k)
          ut(j,k)=ut(j,k)+
     >         0.5d0*((cv0*ot0+cv1*ot1)/(yt(j)*yt1(j))
     >             -(cw0*or0+cw1*or1))
        end do
      end do
      do k=1,Km
        do j=1,Jm-1
          cw0=0.5d0*(w(j,k-1)+w(j+1,k-1))
          cw1=0.5d0*(w(j,k)+w(j+1,k))
          ox0=obx(j,k-1)
          ox1=obx(j,k)
          cu1=0.5d0*(u(j,k)+u(j+1,k))
          ot1=obt(j,k)
          vt(j,k)=vt(j,k)+
     >         0.5d0*(cw0*ox0+cw1*ox1)-cu1*ot1
        end do
        vt(Jm,k)=0.d0
      end do
      do k=1,Km
        do j=1,Jm
          cu0=0.5d0*(u(j,k)+u(j,k+1))
          or0=obr(j,k)
          cv0=0.5d0*(v(j-1,k)+v(j-1,k+1))
          cv1=0.5d0*(v(j,k)+v(j,k+1))
          ox0=rt(j-1)*rt1(j-1)*obx(j-1,k)
          ox1=rt(j)*rt1(j)*obx(j,k)
          wt(j,k)=wt(j,k)+
     >       cu0*or0-0.5d0*(cv0*ox0+cv1*ox1)/(yt(j)*yt1(j))
        end do
      end do
*
* Viscous terms -rot(om)/Re
      do j=1,Jm
        do k=1,Km
          cot0=ot(j-1,k)
          cot1=ot(j,k)
          cor0=or(j,k-1)
          cor1=or(j,k)
          ut(j,k)=ut(j,k)-
     >             ((rt(j)*cot1-rt(j-1)*cot0)/yt1(j)
     >             -(cor1-cor0)/ht)/yt(j)/Re
        end do
      end do
      do k=1,Km
        do j=1,Jm-1
          cox0=ox(j,k-1)
          cox1=ox(j,k)
          cot=ot(j,k)
          vt(j,k)=vt(j,k) -
     >       ( (cox1-cox0)/(rt(j)*ht) - ci*alpha*cot )/Re
        end do
      end do
      do k=1,Km
        do j=1,Jm
          cox0=ox(j-1,k)
          cox1=ox(j,k)
          cor=or(j,k)
          wt(j,k)=wt(j,k) -
     >        ( ci*alpha*cor - (cox1-cox0)/yt1(j) )/Re
        end do
      end do
      return
      end
