*
      subroutine step(t,dt,tol,u,v,w,u1,v1,w1,u2,v2,w2,u3,v3,w3
     >     ,ox,or,ot,p,q,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,u1(0:Imax,0:Jmax,0:*)
     >,v1(0:Imax,0:Jmax,0:*)
     >,w1(0:Imax,0:Jmax,0:*)
     >,u2(0:Imax,0:Jmax,0:*)
     >,v2(0:Imax,0:Jmax,0:*)
     >,w2(0:Imax,0:Jmax,0:*)
     >,u3(0:Imax,0:Jmax,0:*)
     >,v3(0:Imax,0:Jmax,0:*)
     >,w3(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,q(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Dp/Dp
     >/cf/cf
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
          do i=1,Im
            ux=(u(i+1,j,k)-u(i-1,j,k))/(2.*hx)
            vx=(v(i+1,j,k)-v(i-1,j,k))/(2.*hx)
            wx=(w(i+1,j,k)-w(i-1,j,k))/(2.*hx)
            u1(i,j,k)=dt23*(u1(i,j,k)+cf*ux)
            v1(i,j,k)=dt23*(v1(i,j,k)+cf*vx)
            w1(i,j,k)=dt23*(w1(i,j,k)+cf*wx)
          end do
        end do
      end do
      call lin(tau,u1,v1,w1,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u1(i,j,k)=u1(i,j,k)+u(i,j,k)
            v1(i,j,k)=v1(i,j,k)+v(i,j,k)
            w1(i,j,k)=w1(i,j,k)+w(i,j,k)
          end do
        end do
      end do
* 2-nd step
      call rp(t1,u1,v1,w1,u2,v2,w2,ox,or,ot,Imax,Jmax)
      p(0,0,1)=0.d0
      call pres(u2,v2,w2,p,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            ux=(u1(i+1,j,k)-u1(i-1,j,k))/(2.*hx)
            vx=(v1(i+1,j,k)-v1(i-1,j,k))/(2.*hx)
            wx=(w1(i+1,j,k)-w1(i-1,j,k))/(2.*hx)
            u2(i,j,k)=dt13*(u2(i,j,k)+cf*ux)
     >               -(u1(i,j,k)-u(i,j,k))
            v2(i,j,k)=dt13*(v2(i,j,k)+cf*vx)
     >               -(v1(i,j,k)-v(i,j,k))
            w2(i,j,k)=dt13*(w2(i,j,k)+cf*wx)
     >               -(w1(i,j,k)-w(i,j,k))
          end do
        end do
      end do
      call lin(tau,u2,v2,w2,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u2(i,j,k)=u2(i,j,k)+c32*u1(i,j,k)-c12*u(i,j,k)
            v2(i,j,k)=v2(i,j,k)+c32*v1(i,j,k)-c12*v(i,j,k)
            w2(i,j,k)=w2(i,j,k)+c32*w1(i,j,k)-c12*w(i,j,k)
          end do
        end do
      end do
* 3-rd step
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u3(i,j,k)=c34*(u2(i,j,k)-u(i,j,k))
            v3(i,j,k)=c34*(v2(i,j,k)-v(i,j,k))
            w3(i,j,k)=c34*(w2(i,j,k)-w(i,j,k))
          end do
        end do
      end do
      call lin(tau,u3,v3,w3,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u3(i,j,k)=u3(i,j,k)+c32*u2(i,j,k)-c34*u1(i,j,k)+c14*u(i,j,k)
            v3(i,j,k)=v3(i,j,k)+c32*v2(i,j,k)-c34*v1(i,j,k)+c14*v(i,j,k)
            w3(i,j,k)=w3(i,j,k)+c32*w2(i,j,k)-c34*w1(i,j,k)+c14*w(i,j,k)
          end do
        end do
      end do
* 4-th step
      call rp(t1,u2,v2,w2,u1,v1,w1,ox,or,ot,Imax,Jmax)
      call gradp(u1,v1,w1,p,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            ux=(u2(i+1,j,k)-u2(i-1,j,k))/(2.*hx)
            vx=(v2(i+1,j,k)-v2(i-1,j,k))/(2.*hx)
            wx=(w2(i+1,j,k)-w2(i-1,j,k))/(2.*hx)
            u1(i,j,k)=dt34*(u1(i,j,k)+cf*ux)
     >          -(u3(i,j,k)-c38*u2(i,j,k)-c58*u(i,j,k))
            v1(i,j,k)=dt34*(v1(i,j,k)+cf*vx)
     >          -(v3(i,j,k)-c38*v2(i,j,k)-c58*v(i,j,k))
            w1(i,j,k)=dt34*(w1(i,j,k)+cf*wx)
     >          -(w3(i,j,k)-c38*w2(i,j,k)-c58*w(i,j,k))
          end do
        end do
      end do
      call lin(tau,u1,v1,w1,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u1(i,j,k)=u1(i,j,k)+c12*u3(i,j,k)+c34*u2(i,j,k)-c14*u(i,j,k)
            v1(i,j,k)=v1(i,j,k)+c12*v3(i,j,k)+c34*v2(i,j,k)-c14*v(i,j,k)
            w1(i,j,k)=w1(i,j,k)+c12*w3(i,j,k)+c34*w2(i,j,k)-c14*w(i,j,k)
          end do
        end do
      end do
      q(0,0,1)=c12
      call pres(u1,v1,w1,q,Imax,Jmax)
* Accuracy estimation
      err=0.
      ss=0.
      eru=0.
      erv=0.
      erw=0.
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            uu=u1(i,j,k)-u3(i,j,k)
            vv=v1(i,j,k)-v3(i,j,k)
            ww=w1(i,j,k)-w3(i,j,k)
              if(abs(uu).gt.eru)then
                eru=abs(uu)
                iu=i
                ju=j
                ku=k
              end if
              if(abs(vv).gt.erv)then
                erv=abs(vv)
                iv=i
                jv=j
                kv=k
              end if
              if(abs(ww).gt.erw)then
                erw=abs(ww)
                iw=i
                jw=j
                kw=k
              end if
*            err=err+(uu**2+vv**2+ww**2)*yt(j)*yt1(j)
*            ss=ss+yt(j)*yt1(j)
            err=max(err,abs(uu),abs(vv),abs(ww))
          end do
        end do
      end do
          write(*,*)' eru=',eru,' iu=',iu,' ju=',ju,' ku=',ku
          write(*,*)' erv=',erv,' iv=',iv,' jv=',jv,' kv=',kv
          write(*,*)' erw=',erw,' iw=',iw,' jw=',jw,' kw=',kw
*      err=sqrt(err/ss)
      fac=(tol/err)**c13
      if(fac.lt.facmin) then
        dt=dt*fac
        write(*,*)'  STEP:  fac=',fac,'  dt=',dt
        call rp(t,u,v,w,u1,v1,w1,ox,or,ot,Imax,Jmax)
        p(0,0,1)=0.d0
        call pres(u1,v1,w1,p,Imax,Jmax)
        goto 1
      end if
      fac=min(fac,facmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u1(i,j,k)
            v(i,j,k)=v1(i,j,k)
            w(i,j,k)=w1(i,j,k)
          end do
        end do
      end do
      t=t+dt
      dt=fac*dt
      call rp(t,u,v,w,u1,v1,w1,ox,or,ot,Imax,Jmax)
      p(0,0,1)=0.d0
      call pres(u1,v1,w1,p,Imax,Jmax)
      Dp=p(0,0,0)
      return
      end
