*
      subroutine prt(t,dt,u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
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
     >/cf/cf
*
      Ss=0.d0
      ubulk=0.d0
      volume=0.d0
      enrg=0.d0
      dd=0.d0
      do k=1,Km
        do j=1,Jm
          Ss=Ss+yt(j)*ht*yt1(j)
          ubulk=ubulk+u(Im,j,k)*yt(j)*ht*yt1(j)
          do i=1,Im
            uP=1.d0-yt(j)**2
            volume=volume+yt(j)*ht*yt1(j)*hx
            enrg=enrg+(u(i,j,k)-uP)**2*yt(j)*ht*yt1(j)*hx
     >                    +v(i,j,k)**2*rt(j)*ht*rt1(j)*hx
     >                    +w(i,j,k)**2*yt(j)*ht*yt1(j)*hx
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      ubulk=ubulk/Ss
      enrg=enrg/2

      amp=0.d0
      do j=1,Jm
        u0=0.d0
        do k=1,Km
          do i=1,Im
            u0=u0+u(i,j,k)
          end do
        end do
        u0=u0/(Im*Km)
        do k=1,Km
          do i=1,Im
            amp=amp+yt(j)*ht*yt1(j)*hx*(u(i,j,k)-u0)**2
     >             +rt(j)*ht*rt1(j)*hx*v(i,j,k)**2
     >             +yt(j)*ht*yt1(j)*hx*v(i,j,k)**2
          end do
        end do
      end do
      amp=sqrt(amp)

      ucl=0.d0
      do k=1,Km
        do i=1,Im
          ucl=ucl+u(i,1,k)
        end do
      end do
      ucl=ucl/(Im*Km)
*
      write(8,120)t,dt,amp,enrg,ucl,Dp,cf,ubulk,dd
      write(*,110)t,dt,amp,enrg,ucl,Dp,cf,ubulk,dd
120   format(15e25.15)
110   format('t=',f10.4,',dt=',f10.4,',amp=',e12.4,',enr=',e12.4,',Ucl='
     > ,e12.4,',Dp=',e12.4,',cf=',e12.4,',ub=',e12.4,',dd=',e12.4)
      return
      end
