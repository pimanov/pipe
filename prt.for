*
      subroutine prt(t,dt,u,v,w,p,Imax,Jmax)
*      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
     >/dim/Xmax,epsr
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/Dp/Dp
     >/cf/cf,icf
*
      ubulk=0.
      Ss=0.
      amp=0.
      dd=0.
      do j=1,Jm
        Ss=Ss+yt(j)*yt1(j)
        u0=0.
        do k=1,Km
          ubulk=ubulk+u(1,j,k)*yt(j)*yt1(j)
          do i=1,Im
            u0=u0+u(i,j,k)
          end do
        end do
        u0=u0/(Im*Km)
        do k=1,Km
          do i=1,Im
            amp=amp
     >       +((u(i,j,k)-u0)**2+w(i,j,k)**2+v(i,j,k)**2)*yt(j)*yt1(j)
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      ubulk=ubulk/(Ss*Km)
      amp=sqrt(amp/(Ss*Im*Km))
      ucl=0.
      do k=1,Km
        do i=1,Im
          ucl=ucl+u(i,1,k)
        end do
      end do
      ucl=ucl/(Im*Km)
*
      write(8,120)t,dt,Dp,amp,ucl,dd,cf
      write(*,110)t,dt,Dp,amp,ucl,dd,cf
120   format(1pe14.6,15e12.4)
110   format('    t=',f10.4,'  dt=',f10.4,'  Dp=',1pe12.4,
     > '  amp=',e12.4,'  Ucl=',e12.4,'  Div=',e12.4,'  cf=',e12.4)
      return
      end
