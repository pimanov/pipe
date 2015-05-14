*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Jmax=1025, Kmax=1025)
      character*24 fncp,fnbcp,fndat
      complex*16 u,v,w,u1,v1,w1,u2,v2,w2,ox,or,ot,p,q,c0,ci
      dimension
     > ub(0:Jmax,0:Kmax)
     >,vb(0:Jmax,0:Kmax)
     >,wb(0:Jmax,0:Kmax)
     >,obx(0:Jmax,0:Kmax)
     >,obr(0:Jmax,0:Kmax)
     >,obt(0:Jmax,0:Kmax)
     >,u(0:Jmax,0:Kmax)
     >,v(0:Jmax,0:Kmax)
     >,w(0:Jmax,0:Kmax)
     >,u1(0:Jmax,0:Kmax)
     >,v1(0:Jmax,0:Kmax)
     >,w1(0:Jmax,0:Kmax)
     >,u2(0:Jmax,0:Kmax)
     >,v2(0:Jmax,0:Kmax)
     >,w2(0:Jmax,0:Kmax)
     >,u3(0:Jmax,0:Kmax)
     >,v3(0:Jmax,0:Kmax)
     >,w3(0:Jmax,0:Kmax)
     >,ox(0:Jmax,0:Kmax)
     >,or(0:Jmax,0:Kmax)
     >,ot(0:Jmax,0:Kmax)
     >,p(0:Jmax,0:Kmax)
     >,q(0:Jmax,0:Kmax)
      common
     >/dim/epsr,nsym
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/cf/cf
     >/alpha/alpha
*
      c0=(0.d0,0.d0)
      ci=(0.d0,1.d0)
*
      open(5,file='pipe.car')
      read(5,*) tol
      read(5,*) nprt
      read(5,*) nwrt
      read(5,*) tmax
      read(5,*) dtmax
      read(5,*) cf
      read(5,*) fncp
      read(5,*) fnbcp
      read(5,*) fndat
      close(5)
     
      open(9,file=fnbcp,form='unformatted',status='old',err=333)
      read(9)epsr,Jm,lt,nsym
      call com
      if(Jm.gt.Jmax-1) then
        write(*,*)'  Jm=',Jm,'  is greater than   Jmax-1=',Jmax-1
        goto 111
      end if
      if(Km.gt.Kmax-1.or.Km.gt.2048) then ! fft <= 2048
        write(*,*)'  Km=',Km,'  is greater than   Kmax=',Kmax
        goto 111
      end if
      do k=1,Km
        read(9)(ub(j,k),j=1,Jm)
        read(9)(vb(j,k),j=1,Jm)
        read(9)(wb(j,k),j=1,Jm)
      end do
      close(9)

      open(9,file=fncp,form='unformatted',status='old',err=444)
      read(9)t,dt,Re,epsr1,Jm1,lt1,nsym1,alpha
      if(epsr1.ne.epsr) goto 555
      if(Jm1.ne.Jm) goto 555
      if(lt1.ne.lt) goto 555
      if(nsym1.ne.nsym) goto 555
      do k=1,Km
        read(9)(u(j,k),j=1,Jm)
        read(9)(v(j,k),j=1,Jm)
        read(9)(w(j,k),j=1,Jm)
      end do
      close(9)
      dt=min(dt,dtmax)

      write(*,*)' ***************************************************'
      write(*,200) t,dt,Re,epsr,Jm,Km,nsym,alpha
      write(*,*)' ***************************************************'
200   format('    t=',e10.3,' dt=',e9.2,' Re=',e9.2,/,
     >'    epsr=',e9.2,' Jm=',i4,' Km=',i4,' Nsym=',i3,' alpha=',e9.3)
*
      call rp_base(t,ub,vb,wb,obx,obr,obt,Jmax)
      call rp(t,ub,vb,wb,obx,obr,obt,u,v,w,u1,v1,w1,ox,or,ot,Jmax)
      call pres(u1,v1,w1,p,Jmax)

      open(8,file=fndat,access='append')

      call servis(t,u,v,w,ox,or,ot,p,0,Jmax)
      call prt(t,dt,u,v,w,p,Jmax)
      lprt=0
      lwrt=0
10    continue
      call step(t,dt,tol,ub,vb,wb,obx,obr,obt,
     > u,v,w,u1,v1,w1,u2,v2,w2,u3,v3,w3,ox,or,ot,p,q,Jmax)
      lprt=lprt+1
      lwrt=lwrt+1
      dt=min(dt,dtmax)
      call servis(t,u,v,w,ox,or,ot,p,0,Jmax)
      if(lprt.ge.nprt.or.t.ge.tmax) then
        lprt=0
        call servis(t,u,v,w,ox,or,ot,p,1,Jmax)
        call prt(t,dt,u,v,w,p,Jmax)
      end if
      if(lwrt.ge.nwrt.or.t+0.01*dt.ge.tmax) then
        lwrt=0
        call servis(t,u,v,w,ox,or,ot,p,2,Jmax)
        open(9,file=fncp,form='unformatted')
        write(9)t,dt,Re,epsr,Jm,lt,nsym,alpha
        do k=1,Km
          write(9)(u(j,k),j=1,Jm)
          write(9)(v(j,k),j=1,Jm)
          write(9)(w(j,k),j=1,Jm)
        end do
        close(9)
        close(8)
        write(*,*)'*************************************************'
        write(*,*)'*               Control Point                   *'
        write(*,*)'*************************************************'
        open(8,file=fndat,access='append')
      end if 
      if(t+0.01*dt.lt.tmax) goto 10
*
      stop
333   write(*,*)'  File ',fnbcp,' was not found'
      stop
444   write(*,*)'  File ',fncp,' was not found'
      stop
555   write(*,*)'  Grids are differ'
111   stop
      end
