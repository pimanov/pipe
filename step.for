*
      subroutine step(t,dt,tol,ub,vb,wb,obx,obr,obt,
     > u,v,w,u1,v1,w1,u2,v2,w2,u3,v3,w3,ox,or,ot,p,q,Jmax)
      implicit real*8 (a-h,o-z)
      complex*8 u,v,w,u1,v1,w1,u2,v2,w2,u3,v3,w3,ox,or,ot,p,q,c0,ci,
     > uu,vv,ww
      dimension
     > ub(0:Jmax,0:*)
     >,vb(0:Jmax,0:*)
     >,wb(0:Jmax,0:*)
     >,ox(0:Jmax,0:*)
     >,or(0:Jmax,0:*)
     >,ot(0:Jmax,0:*)
     >,u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,u1(0:Jmax,0:*)
     >,v1(0:Jmax,0:*)
     >,w1(0:Jmax,0:*)
     >,u2(0:Jmax,0:*)
     >,v2(0:Jmax,0:*)
     >,w2(0:Jmax,0:*)
     >,u3(0:Jmax,0:*)
     >,v3(0:Jmax,0:*)
     >,w3(0:Jmax,0:*)
     >,p(0:Jmax,0:*)
     >,q(0:Jmax,0:*)
      common
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/cf/cf
     >/alpha/alpha
*
      c12=1.d0/2.d0
      c13=1.d0/3.d0
      c23=2.d0/3.d0
      c14=1.d0/4.d0
      c34=3.d0/4.d0
      c32=3.d0/2.d0
      c38=3.d0/8.d0
      c58=5.d0/8.d0
      facmin=0.5d0
      facmax=1.5d0
      ci=(0.d0,1.d0)
      c0=(0.d0,0.d0)
*
1     continue
      tau=c13*dt
      t1=t+c23*dt
      t2=t+dt
      dt23=c23*dt
      dt13=c13*dt
      dt34=c34*dt
* 1-st step
      do k=1,Km
        do j=1,Jm
          ux=ci*alpha*u(j,k)
          vx=ci*alpha*v(j,k)
          wx=ci*alpha*w(j,k)
          u1(j,k)=dt23*(u1(j,k)+cf*ux)
          v1(j,k)=dt23*(v1(j,k)+cf*vx)
          w1(j,k)=dt23*(w1(j,k)+cf*wx)
        end do
      end do
      call lin(tau,u1,v1,w1,Jmax)
      do k=1,Km
        do j=1,Jm
          u1(j,k)=u1(j,k)+u(j,k)
          v1(j,k)=v1(j,k)+v(j,k)
          w1(j,k)=w1(j,k)+w(j,k)
        end do
      end do
* 2-nd step
*      call rp(t1,u1,v1,w1,u2,v2,w2,ox,or,ot,buf,Imax,Jmax)
      call rp(t1,ub,vb,wb,obx,obr,obt,u1,v1,w1,u2,v2,w2,ox,or,ot,Jmax)
      p(0,1)=c0
      call pres(u2,v2,w2,p,Jmax)
      do k=1,Km
        do j=1,Jm
            ux=ci*alpha*u1(j,k)
            vx=ci*alpha*v1(j,k)
            wx=ci*alpha*w1(j,k)
            u2(j,k)=dt13*(u2(j,k)+cf*ux)-(u1(j,k)-u(j,k))
            v2(j,k)=dt13*(v2(j,k)+cf*vx)-(v1(j,k)-v(j,k))
            w2(j,k)=dt13*(w2(j,k)+cf*wx)-(w1(j,k)-w(j,k))
        end do
      end do
      call lin(tau,u2,v2,w2,Jmax)
      do k=1,Km
        do j=1,Jm
            u2(j,k)=u2(j,k)+c32*u1(j,k)-c12*u(j,k)
            v2(j,k)=v2(j,k)+c32*v1(j,k)-c12*v(j,k)
            w2(j,k)=w2(j,k)+c32*w1(j,k)-c12*w(j,k)
        end do
      end do
* 3-rd step
      do k=1,Km
        do j=1,Jm
          u3(j,k)=c34*(u2(j,k)-u(j,k))
          v3(j,k)=c34*(v2(j,k)-v(j,k))
          w3(j,k)=c34*(w2(j,k)-w(j,k))
        end do
      end do
      call lin(tau,u3,v3,w3,Jmax)
      do k=1,Km
        do j=1,Jm
          u3(j,k)=u3(j,k)+c32*u2(j,k)-c34*u1(j,k)+c14*u(j,k)
          v3(j,k)=v3(j,k)+c32*v2(j,k)-c34*v1(j,k)+c14*v(j,k)
          w3(j,k)=w3(j,k)+c32*w2(j,k)-c34*w1(j,k)+c14*w(j,k)
        end do
      end do
* 4-th step
*      call rp(t1,u2,v2,w2,u1,v1,w1,ox,or,ot,buf,Imax,Jmax)
      call rp(t1,ub,vb,wb,obx,obr,obt,u2,v2,w2,u1,v1,w1,ox,or,ot,Jmax)
      call gradp(u1,v1,w1,p,Jmax)
      do k=1,Km
        do j=1,Jm
          ux=ci*alpha*u2(j,k)
          vx=ci*alpha*v2(j,k)
          wx=ci*alpha*w2(j,k)
          u1(j,k)=dt34*(u1(j,k)+cf*ux)-(u3(j,k)-c38*u2(j,k)-c58*u(j,k))
          v1(j,k)=dt34*(v1(j,k)+cf*vx)-(v3(j,k)-c38*v2(j,k)-c58*v(j,k))
          w1(j,k)=dt34*(w1(j,k)+cf*wx)-(w3(j,k)-c38*w2(j,k)-c58*w(j,k))
        end do
      end do
      call lin(tau,u1,v1,w1,Jmax)
      do k=1,Km
        do j=1,Jm
          u1(j,k)=u1(j,k)+c12*u3(j,k)+c34*u2(j,k)-c14*u(j,k)
          v1(j,k)=v1(j,k)+c12*v3(j,k)+c34*v2(j,k)-c14*v(j,k)
          w1(j,k)=w1(j,k)+c12*w3(j,k)+c34*w2(j,k)-c14*w(j,k)
        end do
      end do
      q(0,1)=c0
      call pres(u1,v1,w1,q,Jmax)
* Accuracy estimation
      error=0.d0
      do k=1,Km
        do j=1,Jm
          uu=u1(j,k)-u3(j,k)
          vv=v1(j,k)-v3(j,k)
          ww=w1(j,k)-w3(j,k)
          error=max(error,cabs(uu),cabs(vv),cabs(ww))
        end do
      end do
      fac=(tol/error)**c13
      if(fac.lt.facmin) then
        dt=dt*fac
        if(Np.eq.0)write(*,*)'  STEP:  fac=',fac,'  dt=',dt
*        call rp(t,u,v,w,u1,v1,w1,ox,or,ot,buf,Imax,Jmax)
        call rp(t,ub,vb,wb,obx,obr,obt,u,v,w,u1,v1,w1,ox,or,ot,Jmax)
        p(0,1)=c0
        call pres(u1,v1,w1,p,Jmax)
        goto 1
      end if
      fac=min(fac,facmax)
      do k=1,Km
        do j=1,Jm
          u(j,k)=u1(j,k)
          v(j,k)=v1(j,k)
          w(j,k)=w1(j,k)
        end do
      end do
      t=t+dt
      dt=fac*dt
*      call rp(t,u,v,w,u1,v1,w1,ox,or,ot,buf,Imax,Jmax)
      call rp(t,ub,vb,wb,obx,obr,obt,u,v,w,u1,v1,w1,ox,or,ot,Jmax)
      p(0,1)=c0
      call pres(u1,v1,w1,p,Jmax)
      return
      end
