*
*     program init
      implicit real*8 (a-h,o-z)
      parameter (Imax=1025, Jmax=65, Kmax=129)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/cf/cf
*
      t=0.0
      open(5,file='init.car')
      read(5,*) Xmax, lx
      read(5,*) epsr, Jm
      read(5,*) nsym, lt
      read(5,*) Re
      read(5,*) dt
      read(5,*) a
      read(5,*) fncp
      close(5)
      call com
      write(*,*)' ***************************************************'
      write(*,200) t,dt,a,Re,Xmax,epsr,Imm,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    t=',1pe10.3,' dt=',e9.2,' amp=',e9.2,/,
     >'    Re=',e9.2,/,
     >'    Xmax=',e9.2,/,
     >'    epsr=',e9.2,' Im=',i4,' Jm=',i4,' Km=',i4,' Nsym=',i3)
      if(Im.ge.Imax) stop 'Im >= Imax'
      if(Im.gt.2048) stop 'Im > 2048 (fft,rlx)'
      if(Jm.ge.Jmax) stop 'Jm >= Jmax'
      if(Jm.gt.128) stop 'Jm > 128 (dimr,..)'
      if(Km.ge.Kmax) stop 'Km >= Kmax'
      if(Km.gt.256) stop 'Km > 256 (rlt)'
*
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=sin(0.1*i+0.01*j+0.3*k+14.0)
            v(i,j,k)=cos(0.02*i+0.1*j+0.2*k+22.0)
            w(i,j,k)=sin(0.05*i+0.05*j+0.1*k+10.0)
          end do
        end do
      end do
      p(0,0,1)=0.d0
      call pres(u,v,w,p,Imax,Jmax)
*
      dd=0.0
      do k=2,Km-1
        do j=2,Jm-1
          do i=2,Im-1
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      write(*,*) 'disturb central div = ',dd

      amp=0.0
      vol=0.0
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            amp=amp+u(i,j,k)**2*yt(j)*yt1(j)*ht*hx 
     >             +v(i,j,k)**2*rt(j)*rt1(j)*ht*hx 
     >             +w(i,j,k)**2*yt(j)*yt1(j)*ht*hx 
            vol=vol+yt(j)*yt1(j)*ht*hx
          end do
        end do
      end do
      amp=sqrt(amp/vol)
      write(*,*) 'disturb amp = ', amp

      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u(i,j,k)/amp*a
            v(i,j,k)=v(i,j,k)/amp*a
            w(i,j,k)=w(i,j,k)/amp*a
          end do
        end do
      end do

      amp=0.0
      vol=0.0
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            amp=amp+u(i,j,k)**2*yt(j)*yt1(j)*ht*hx 
     >             +v(i,j,k)**2*rt(j)*rt1(j)*ht*hx 
     >             +w(i,j,k)**2*yt(j)*yt1(j)*ht*hx 
            vol=vol+yt(j)*yt1(j)*ht*hx
          end do
        end do
      end do
      amp=sqrt(amp/vol)
      write(*,*) 'disturb amp = ', amp

      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u(i,j,k)+1.0-yt(j)**2
          end do
        end do
      end do
      p(0,0,1)=0.5d0
      call pres(u,v,w,p,Imax,Jmax)

      dd=0.0
      do k=2,Km-1
        do j=2,Jm-1
          do i=2,Im-1
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      write(*,*) 'Pois+dist central div = ',dd

      open(9,file=fncp,form='unformatted')
      Dp=p(0,0,0)
      write(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
      do k=1,Km
        do j=1,Jm
          write(9)(u(i,j,k),i=1,Im)
          write(9)(v(i,j,k),i=1,Im)
          write(9)(w(i,j,k),i=1,Im)
        end do
      end do
      close(9)
*
      stop
      end
