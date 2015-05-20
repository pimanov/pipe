*
      subroutine prt(t,dt,u,v,w,p,Jmax)
      implicit real*8 (a-h,o-z)
      complex*16 u,v,w,p,d,c0,ci,ubulk,u33,om,ampl,om1,amp
      dimension
     > u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,p(0:Jmax,0:*)
      common
     >/dim/epsr,nsym
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/cf/cf
     >/alpha/alpha
     >/om/u33,tl,ampl,amp1l
*
      p(1,1)=p(1,1)
*
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
          dd=max(dd,abs(d))
        end do
      end do
      ubulk=ubulk/Ss

      amp=c0
      do j=1,Jm
        do k=1,Km
          amp=amp+yt(j)*ht*yt1(j)*u(j,k)**2
     >           +rt(j)*ht*rt1(j)*v(j,k)**2
     >           +yt(j)*ht*yt1(j)*w(j,k)**2
        end do
      end do
      amp=sqrt(amp*nsym)
 
      om=(log(amp)-log(ampl))/(t-tl)
      ampl=amp

      om1=(log(u(Jm/2,Km/2))-log(u33))/(t-tl)
      tl=t
      u33=u(Jm/2,Km/2)

      write(8,120)t,dt,abs(amp),cf,ubulk,dd,om,om1
      write(*,110)t,dt,abs(amp),cf,ubulk,dd,om,om1
      return
120   format(15e25.15)
110   format('t=',f10.4,',  dt=',f10.4,',  amp=',e12.4,',  cf=',f8.4,
     > ',  ub=(',e12.4,', ',e12.4,')',',  dd=',e12.4,
     > ', om=(',f10.4,', ',f10.4,')',', om1=(',f10.4,', ',f10.4,')')
      end
