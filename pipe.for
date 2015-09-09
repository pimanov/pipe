*
*     program pipeSym
      implicit real*8 (a-h,o-z)
      parameter (Imax=1025, Jmax=65, Kmax=129)
      character*24 fncp,fndat
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,u1(0:Imax,0:Jmax,0:Kmax)
     >,v1(0:Imax,0:Jmax,0:Kmax)
     >,w1(0:Imax,0:Jmax,0:Kmax)
     >,u2(0:Imax,0:Jmax,0:Kmax)
     >,v2(0:Imax,0:Jmax,0:Kmax)
     >,w2(0:Imax,0:Jmax,0:Kmax)
     >,u3(0:Imax,0:Jmax,0:Kmax)
     >,v3(0:Imax,0:Jmax,0:Kmax)
     >,w3(0:Imax,0:Jmax,0:Kmax)
     >,ox(0:Imax,0:Jmax,0:Kmax)
     >,or(0:Imax,0:Jmax,0:Kmax)
     >,ot(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
     >,q(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/cf/cf
*
      open(5,file='pipe.car')
      read(5,*) tol
      read(5,*) nprt
      read(5,*) nwrt
      read(5,*) tmax
      read(5,*) dtmax
      read(5,*) cf
      read(5,*) fncp
      read(5,*) fndat
      close(5)
      open(9,file=fncp,form='unformatted',status='old',err=111)
      read(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
      call com
      write(*,*)' ***************************************************'
      write(*,200) t,dt,Dp,Re,Xmax,epsr,Im,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    t=',1pe10.3,' dt=',e9.2,' Dp=',e9.2,/,
     >'    Re=',e9.2,/,
     >'    Xmax=',e9.2,/,
     >'    epsr=',e9.2,' Im=',i4,' Jm=',i4,' Km=',i4,' Nsym=',i3)
*
      if(Im.ge.Imax) stop 'Im >= Imax'
      if(Im.gt.2048) stop 'Im > 2048 (fft,rlx)'
      if(Jm.ge.Jmax) stop 'Jm >= Jmax'
      if(Jm.gt.128) stop 'Jm > 128 (dimr,..)'
      if(Km.ge.Kmax) stop 'Km >= Kmax'
      if(Km.gt.256) stop 'Km > 256 (rlt)'
*
      do k=1,Km
        do j=1,Jm
          read(9)(u(i,j,k),i=1,Im)
          read(9)(v(i,j,k),i=1,Im)
          read(9)(w(i,j,k),i=1,Im)
        end do
      end do
      close(9)
      dt=min(dt,dtmax)
*
      call rp(t,u,v,w,u1,v1,w1,ox,or,ot,Imax,Jmax)
      p(0,0,1)=0.d0
      call pres(u1,v1,w1,p,Imax,Jmax)
*
      open(8,file=fndat,access='append')
*
      call servis(t,u,v,w,ox,or,ot,p,0,Imax,Jmax)
      call prt(t,dt,u,v,w,p,Imax,Jmax)
      lprt=0
      lwrt=0
10    continue
      call step(t,dt,tol,u,v,w,u1,v1,w1,u2,v2,w2,
     >          u3,v3,w3,ox,or,ot,p,q,Imax,Jmax)
      dt=min(dt,dtmax)
      lprt=lprt+1
      lwrt=lwrt+1
      call servis(t,u,v,w,ox,or,ot,p,0,Imax,Jmax)
      if(lprt.ge.nprt.or.t.ge.tmax) then
        lprt=0
        call servis(t,u,v,w,ox,or,ot,p,1,Imax,Jmax)
        call prt(t,dt,u,v,w,p,Imax,Jmax)
      end if
      if(lwrt.ge.nwrt.or.t+0.01*dt.ge.tmax) then
        lwrt=0
        call servis(t,u,v,w,ox,or,ot,p,2,Imax,Jmax)
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
        close(8)
        write(*,*)'*'
        open(8,file=fndat,access='append')
      end if 
      if(t+0.01*dt.lt.tmax) goto 10
*
      stop
111   stop 'file.scp was not found'
      end
