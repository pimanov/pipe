*
      subroutine lin(tau,u,v,w,Jmax)
      implicit real*8 (a-h,o-z)
      complex*8 u,v,w
      dimension
     > u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,ap(256),bp(256),cp(256)
     >,d1(256),e1(256),d2(256),e2(256)
      common
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/pry/apy(128),bpy(128),cpy(128)
     >/prv/apv(128),bpv(128),cpv(128)
     >/prw/apw(128),bpw(128),cpw(128)
     >/Re/Re
      ct=Re/tau
*

*
* At
      do j=1,Jm
        do k=1,Km
          ap(k)=1.d0/ht**2
          cp(k)=ap(k)
          bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
          d1(k)=-ct*real(u(j,k))*yt(j)**2
          d2(k)=-ct*aimag(u(j,k))*yt(j)**2
        end do
        bp(1)=bp(1)+ap(1)
        bp(Km)=bp(Km)+cp(Km)
        call prog3(ap,bp,cp,d1,e1,Km)
        call prog3(ap,bp,cp,d2,e2,Km)
        do k=1,Km
          u(j,k)=cmplx(e1(k),e2(k))
        end do
        do k=1,Km
          ap(k)=1.d0/ht**2
          cp(k)=ap(k)
          bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
          d1(k)=-ct*real(v(j,k))*yt(j)**2
          d2(k)=-ct*aimag(v(j,k))*yt(j)**2
        end do
        bp(1)=bp(1)+ap(1)
        bp(Km)=bp(Km)+cp(Km)
        call prog3(ap,bp,cp,d1,e1,Km)
        call prog3(ap,bp,cp,d2,e2,Km)
        do k=1,Km
          v(j,k)=cmplx(e1(k),e2(k))
        end do
        do k=1,Km-1
          ap(k)=1.d0/ht**2
          cp(k)=ap(k)
          bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
          d1(k)=-ct*real(w(j,k))*yt(j)**2
          d2(k)=-ct*aimag(w(j,k))*yt(j)**2
        end do
        call prog3(ap,bp,cp,d1,e1,Km-1)
        call prog3(ap,bp,cp,d2,e2,Km-1)
        do k=1,Km-1
          w(i,j,k)=cmplx(e1(k),e2(k))
        end do
      end do
*
* Ar
      do k=1,Km
        do j=1,Jm
          bp(j)=bpy(j)-ct
          d1(j)=-ct*real(u(j,k))
          d2(j)=-ct*aimag(u(j,k))
        end do
        bp(Jm)=bp(Jm)-cpy(Jm)
        call prog3(apy,bp,cpy,d1,e1,Jm)
        call prog3(apy,bp,cpy,d2,e2,Jm)
        do j=1,Jm
          u(i,j,k)=cmplx(e1(j),e2(j))
        end do
        do j=1,Jm-1
          bp(j)=bpv(j)-ct
          d1(j)=-ct*rt(j)*real(v(j,k))
          d2(j)=-ct*rt(j)*aimag(v(j,k))
        end do
        call prog3(apv,bp,cpv,d1,e1,Jm-1)
        call prog3(apv,bp,cpv,d2,e2,Jm-1)
        do j=1,Jm-1
          v(i,j,k)=cmplx(e1(j),e2(j))/rt(j)
        end do
        do j=1,Jm
          bp(j)=bpy(j)-ct
          d1(j)=-ct*yt(j)*real(w(j,k))
          d2(j)=-ct*yt(j)*aimag(w(j,k))
        end do
        bp(Jm)=bp(Jm)-cpy(Jm)
        call prog3(apy,bp,cpy,d1,e1,Jm)
        call prog3(apy,bp,cpy,d2,e2,Jm)
        do j=1,Jm
          w(i,j,k)=cmplx(e1(j),e2(j))/yt(j)
        end do
      end do
      return
      end
