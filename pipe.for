*
*     program an
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      parameter (Imax=513, Jmax=129, Kmax=129)
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
     >,buf(2*Imax*Jmax*Kmax)
     >,g1(0:Imax,0:Jmax,0:Kmax)
     >,g2(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
     >/cf/cf
*
      c0=0.d0
      c12=0.5d0
*
      call MPI_INIT(ier)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,Npm,ier)
      call MPI_COMM_RANK(MPI_COMM_WORLD,Np,ier)
*
      if(Np.eq.0) then
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
      end if
      call MPI_BCAST(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(nprt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(nwrt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(tmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(dtmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(cf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
*
      call init_cp(fncp,t,dt,u,v,w,buf,Imax,Jmax,Kmax)
      dt=min(dt,dtmax)
*
      if(Np.eq.0) then
        write(*,*)' ***************************************************'
        write(*,*)' *        Number of processors =',Npm,'          *'
        write(*,200) t,dt,Re,Xmax,epsr,Imm,Jm,Km,dsym
        write(*,*)' ***************************************************'
      endif
200   format('    t=',1pe10.3,' dt=',e9.2,/,
     >'    Re=',e9.2,' Xmax=',e9.2,' epsr=',e9.2,/,
     >'    Im=',i4,' Jm=',i4,' Km=',i4,' nsym=',e9.2)
*
      call rp(t,u,v,w,u1,v1,w1,ox,or,ot,buf,Imax,Jmax)
      call pres(u1,v1,w1,p,c0,buf,Imax,Jmax)
*
      if(Np.eq.0) open(8,file=fndat,access='append')
      call servis(t,g1,g2,u,v,w,ox,or,ot,p,0,buf,Imax,Jmax)
      call prt(t,dt,u,v,w,p,Imax,Jmax)
      lprt=0
      lwrt=0
*
10    continue
        call step(t,dt,tol,u,v,w,u1,v1,w1,u2,v2,w2
     >   ,u3,v3,w3,ox,or,ot,p,q,buf,Imax,Jmax)
        lprt=lprt+1
        lwrt=lwrt+1
        dt=min(dt,dtmax)
        call servis(t,g1,g2,u,v,w,ox,or,ot,p,0,buf,Imax,Jmax)
        if(lprt.ge.nprt.or.t+0.01*dt.ge.tmax) then
          lprt=0
          call servis(t,g1,g2,u,v,w,ox,or,ot,p,1,buf,Imax,Jmax)
          call prt(t,dt,u,v,w,p,Imax,Jmax)
        end if
        if(lwrt.ge.nwrt.or.t+0.01*dt.ge.tmax) then
          lwrt=0
          call servis(t,g1,g2,u,v,w,ox,or,ot,p,2,buf,Imax,Jmax)
          Dp=p(0,0,0)
          call write_cp(t,dt,Dp,u,v,w,fncp,buf,Imax,Jmax)
          if(Np.eq.0)then
            close(8)
            write(*,*)'************************************************'
            write(*,*)'*              Control Point                   *'
            write(*,*)'************************************************'
            open(8,file=fndat,access='append')
          end if
        end if
      if(t+0.01*dt.lt.tmax) goto 10
*
      call MPI_FINALIZE(ier)
      stop
      end
