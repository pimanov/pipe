*
      subroutine prt(t,dt,u,v,w,p,Jmax)
      implicit real*8 (a-h,o-z)
      complex*8 u,v,w,p,d,c0,ci
      dimension
     > u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,p(0:Jmax,0:*)
      common
     >/dim/Xmax,epsr,nsym
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/cf/cf
     >/alpha/alpha
*
      Dp=p(0,0,0)
      c0=(0.d0,0.d0)
      ci=(0.d0,1.d0)

      Ss=0.d0
      ubulk=c0
      dd=0.d0
      do k=1,Km
        do j=1,Jm
          Ss=Ss+yt(j)*ht*yt1(j)
          ubulk=ubulk+u(j,k)*yt(j)*ht*yt1(j)
          call div(j,k,u,v,w,d,Jmax)
          dd=max(dd,cabs(d))
        end do
      end do
      ubulk=ubulk/Ss

      amp1=0.d0
      amp2=0.d0
      do j=1,Jm
        do k=1,Km
          amp1=amp1+yt(j)*ht*yt1(j)*hx*real(u(j,k))**2
     >             +rt(j)*ht*rt1(j)*hx*real(v(j,k))**2
     >             +yt(j)*ht*yt1(j)*hx*real(w(j,k))**2
          amp2=amp2+yt(j)*ht*yt1(j)*hx*aimag(u(j,k))**2
     >             +rt(j)*ht*rt1(j)*hx*aimag(v(j,k))**2
     >             +yt(j)*ht*yt1(j)*hx*aimag(w(j,k))**2
        end do
      end do
      amp=sqrt(amps*2*nsym)

      write(8,120)t,dt,amp1,amp2,real(Dp),aimag(Dp),cf,
     > real(ubulk),aimag(ubulk),dd
      write(*,110)t,dt,amp1,amp2,real(Dp),aimag(Dp),cf,
     > real(ubulk),aimag(ubulk),dd
      return
120   format(15e25.15)
110   format('t=',f10.4,',dt=',f10.4,',Re(amp)=',e12.4,',Im(amp)',e12.4,
     > ',Re(Dp)=',e12.4,',Im(Dp)=',e12.4,',cf=',e12.4,',Re(ub)=',e12.4,
     > ',Im(ub)=',e12.4,',dd=',e12.4)
      end
