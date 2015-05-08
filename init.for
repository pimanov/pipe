*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Jmax=129, Kmax=129)
      character*12 fncp
      complex*16 u,v,w,p,c0,ci,Dp
      dimension
     > u(0:Jmax,0:Kmax)
     >,v(0:Jmax,0:Kmax)
     >,w(0:Jmax,0:Kmax)
     >,p(0:Jmax,0:Kmax)
      common
     >/dim/epsr,nsym
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/alpha/alpha
*
      c0=(0.d0,0.d0)
      ci=(0.d0,0.d0)
*
      open(5,file='init.car')
      read(5,*) Jm,epsr
      read(5,*) lt
      read(5,*) nsym
      read(5,*) Re
      read(5,*) amp
      read(5,*) dt
      read(5,*) alpha
      read(5,*) fncp
      close(5)
     
      call com
      if(Jm.gt.Jmax-1.or.Jm.gt.128) then
        write(*,*)'  Jm=',Jm,'  is greater than   Jmax-1=',Jmax-1
        goto 111
      end if
      if(Km.gt.Kmax-1.or.Km.gt.256) then
        write(*,*)'  Km=',Km,'  is greater than   Kmax=',Kmax
        goto 111
      end if


      write(*,*)' ***************************************************'
      write(*,200) dt,Re,epsr,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    dt=',e9.2,'    Re=',e9.2,/,
     >'    epsr=',e9.2,' Jm=',i4,' Km=',i4,' Nsym=',i3)
*

      do j=1,Jm
        do k=1,Km
          u(j,k)=sin(2.0*j/Jm)+ci*cos(3.0*k/Km)
          v(j,k)=sin(2.0*k/Km)+ci*sin(3.0*j/Jm)
          w(j,k)=sin(3.0*j/Jm)*cos(3.0*k/Km)+ci*0.5
        end do
      end do

      p(0,1)=c0
      call pres(u,v,w,p,Imax)

      do j=1,Jm
        do k=1,Km
          a=(abs(u(j,k))+abs(w(j,k)))*ht*yt(j)*yt1(j)
     >                   +abs(v(j,k))*ht*rt(j)*rt1(j)
        end do
      end do
      write(*,*) 'a=',a

      do j=1,Jm
        do k=1,Km
          u(j,k)=u(j,k)*amp/a
          v(j,k)=v(j,k)*amp/a
          w(j,k)=w(j,k)*amp/a
        end do
      end do

      Dp=c0
      t=0.d0
      open(9,file=fncp,form='unformatted')
      write(9)t,dt,Dp,Re,epsr,Jm,lt,nsym,alpha
      do k=1,Km
        write(9)(u(j,k),j=1,Jm)
        write(9)(v(j,k),j=1,Jm)
        write(9)(w(j,k),j=1,Jm)
      end do
      close(9)

111   stop
      end
