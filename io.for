*
*
      subroutine init_cp(fncp,t,dt,Dp,u,v,w,buf,Imax,Jmax,Kmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,buf(*)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
*
      if(Np.eq.0)then
        open(9,file=fncp,form='unformatted',status='old',err=111)
        read(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,dsym
      end if
      call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Dp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Re,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Xmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(epsr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(lx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Jm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(lt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(dsym,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
*
      call com
*
      if(Np.eq.0)then
        if(Im.gt.Imax-1.or.Im*Npm.gt.2048) then
          write(*,*)'  Im=',Im,'  is greater than   Imax-1=',Imax-1
          goto 333
        end if
        if(Jm.gt.Jmax-1.or.Jm.gt.128) then
          write(*,*)'  Jm=',Jm,'  is greater than   Jmax-1=',Jmax-1
          goto 333
        end if
        if(Km.gt.Kmax-1.or.Km.gt.256) then
          write(*,*)'  Km=',Km,'  is greater than   Kmax=',Kmax
          goto 333
        end if
      end if
*
      do k=1,Km
        do j=1,Jm
          if(Np.eq.0) then
            read(9)(buf(i),i=1,Imm)
            read(9)(buf(i+Imm),i=1,Imm)
            read(9)(buf(i+2*Imm),i=1,Imm)
            do i=1,Im
              u(i,j,k)=buf(i)
              v(i,j,k)=buf(i+Imm)
              w(i,j,k)=buf(i+2*Imm)
            end do
            do n=1,Npm-1
              call MPI_SEND(buf(n*Im+1),Im,MPI_DOUBLE_PRECISION
     >                     ,n,1,MPI_COMM_WORLD,ier)
              call MPI_SEND(buf(n*Im+1+Imm),Im,MPI_DOUBLE_PRECISION
     >                     ,n,2,MPI_COMM_WORLD,ier)
              call MPI_SEND(buf(n*Im+1+2*Imm),Im,MPI_DOUBLE_PRECISION
     >                     ,n,3,MPI_COMM_WORLD,ier)
            end do
          else
            call MPI_RECV(u(1,j,k),Im,MPI_DOUBLE_PRECISION,0,1
     >                   ,MPI_COMM_WORLD,status,ier)
            call MPI_RECV(v(1,j,k),Im,MPI_DOUBLE_PRECISION,0,2
     >                   ,MPI_COMM_WORLD,status,ier)
            call MPI_RECV(w(1,j,k),Im,MPI_DOUBLE_PRECISION,0,3
     >                   ,MPI_COMM_WORLD,status,ier)
          end if
          call MPI_BARRIER(MPI_COMM_WORLD,ier)
        end do
      end do
      if(Np.eq.0) close(9)
*
      return
111   write(*,*) 'file ', fncp,' was not found'
      call MPI_ABORT(111)
333   call MPI_ABORT(333)
      end subroutine
*
*
      subroutine write_cp(t,dt,Dp,u,v,w,fncp,buf,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,buf(*)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
*
      if(Np.eq.0) then
        open(9,file=fncp,form='unformatted')
        write(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,dsym
      end if
*
      do k=1,Km
        do j=1,Jm
          if(Np.eq.0) then
            do i=1,Im
              buf(i)=u(i,j,k)
              buf(i+Imm)=v(i,j,k)
              buf(i+2*Imm)=w(i,j,k)
            end do
            do n=1,Npm-1
              call MPI_RECV(buf(n*Im+1),Im,MPI_DOUBLE_PRECISION
     >                 ,n,1,MPI_COMM_WORLD,status,ier)
              call MPI_RECV(buf(n*Im+1+Imm),Im,MPI_DOUBLE_PRECISION
     >                 ,n,2,MPI_COMM_WORLD,status,ier)
              call MPI_RECV(buf(n*Im+1+2*Imm),Im,MPI_DOUBLE_PRECISION
     >                 ,n,3,MPI_COMM_WORLD,status,ier)
            end do
            write(9)(buf(i),i=1,Imm)
            write(9)(buf(i+Imm),i=1,Imm)
            write(9)(buf(i+2*Imm),i=1,Imm)
          else
            call MPI_SEND(u(1,j,k),Im,MPI_DOUBLE_PRECISION,0,1
     >                 ,MPI_COMM_WORLD,ier)
            call MPI_SEND(v(1,j,k),Im,MPI_DOUBLE_PRECISION,0,2
     >                 ,MPI_COMM_WORLD,ier)
            call MPI_SEND(w(1,j,k),Im,MPI_DOUBLE_PRECISION,0,3
     >                 ,MPI_COMM_WORLD,ier)
          end if
          call MPI_BARRIER(MPI_COMM_WORLD,ier)
        end do
      end do
*
      if(Np.eq.0) close(9)
      return
      end subroutine
