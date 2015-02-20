*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Imax=513, Jmax=33, Kmax=65)
      character*12 fncp,fndat
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/pi/pi
*
      open(5,file='init.car')
      read(5,*) Xmax
      read(5,*) lx
      read(5,*) Jm, epsr
      read(5,*) lt
      read(5,*) Re
      read(5,*) amp
      read(5,*) dt
      read(5,*) fncp
*
      call com
      open(9,file=fncp,form='unformatted',status='new',err=1)
      t=0.d0
      do k=1,Km
        tt=(1.d0*k)/Km
        do j=1,Jm
          r=yt(j)
          do i=1,Im
            x=(1.d0*i)/Im
            u(i,j,k)=(1.-r**2)*r*(sin(2.*pi*(x-tt))
     >       +r*sin(2.*pi*(x+2.*tt)))
            v(i,j,k)=0.
            w(i,j,k)=0.
          end do
        end do
      end do
      dd=0.
      do k=1,Km
      do j=1,Jm
      do i=1,Im
      call div(i,j,k,u,v,w,d,Imax,Jmax)
      dd=max(dd,abs(d))
      end do
      end do
      end do
      write(*,*)' Div0=',dd
      p(0,0,1)=0.d0
      call pres(u,v,w,p,Imax,Jmax)
      dd=0.
      a=0.
      ss=0.
      do k=1,Km
      do j=1,Jm
      do i=1,Im
      call div(i,j,k,u,v,w,d,Imax,Jmax)
      dd=max(dd,abs(d))
      a=a+yt(j)*yt1(j)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
      ss=ss+yt(j)*yt1(j)
      end do
      end do
      end do
      a=sqrt(a/ss)
      write(*,*)' Div =',dd,'  a=',a
      Dp=4.d0/Re
      write(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=1.-yt(j)**2+u(i,j,k)*amp/a
            v(i,j,k)=v(i,j,k)*amp/a
            w(i,j,k)=w(i,j,k)*amp/a
          end do
          write(9)(u(i,j,k),i=1,Im)
          write(9)(v(i,j,k),i=1,Im)
          write(9)(w(i,j,k),i=1,Im)
        end do
      end do
      close(9)
      stop
1     write(*,*)'  File ',fncp,' already exists'
      stop
      end
