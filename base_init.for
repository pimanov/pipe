*
*     program an
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      parameter (Imax=1025, Jmax=129, Kmax=129)
      character*12 fncp,fnbcp
      character*256 comment
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,buf(2*Imax*Jmax*Kmax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
     >/cf/cf
*
      integer status(MPI_STATUS_SIZE)
      call MPI_INIT(ier)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,Npm,ier)
      call MPI_COMM_RANK(MPI_COMM_WORLD,Np,ier)
 
      if(Np.eq.0)then
        open(5,file='base_init.car')
        read(5,*) cf
        read(5,*) fncp
        read(5,*) fnbcp
        read(5,*, err=500) comment
500     write(*,*) 'base_init.car: cf=',cf,' fncp=',fncp,' fnbcp=',fnbcp
        write(*,*) 'comment:',comment
        istop=0
        open(9,file=fncp,form='unformatted',status='old',err=111)
      end if
222   call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(istop.ne.0) goto 333
      if(Np.eq.0) read(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
*
      call MPI_BCAST(cf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
*
      call MPI_BCAST(Xmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(epsr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(lx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Jm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(lt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(nsym,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
*
      call com
      if(Np.eq.0) then      
      write(*,*)' ***************************************************'
      write(*,*)' *        Number of processors =',Npm,'          *'
      write(*,200) t,dt,Dp,Re,Xmax,epsr,Imm,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    t=',1pe10.3,' dt=',e9.2,' Dp=',e9.2,/,
     >'    Re=',e9.2,/,
     >'    Xmax=',e9.2,/,
     >'    epsr=',e9.2,' Im=',i4,' Jm=',i4,' Km=',i4,' Nsym=',i3)
*
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
      if(Np.eq.0) then
        close(9)
        open(9,file=fnbcp,form='unformatted')
        write(9)cf,Xmax,epsr,lx,Jm,lt,nsym
      end if
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
     >                   ,n,1,MPI_COMM_WORLD,status,ier)
              call MPI_RECV(buf(n*Im+1+Imm),Im,MPI_DOUBLE_PRECISION
     >                   ,n,2,MPI_COMM_WORLD,status,ier)
              call MPI_RECV(buf(n*Im+1+2*Imm),Im,MPI_DOUBLE_PRECISION
     >                   ,n,3,MPI_COMM_WORLD,status,ier)
            end do
            write(9)(buf(i),i=1,Imm)
            write(9)(buf(i+Imm),i=1,Imm)
            write(9)(buf(i+2*Imm),i=1,Imm)
          else
            call MPI_SEND(u(1,j,k),Im,MPI_DOUBLE_PRECISION,0,1
     >                   ,MPI_COMM_WORLD,ier)              
            call MPI_SEND(v(1,j,k),Im,MPI_DOUBLE_PRECISION,0,2
     >                   ,MPI_COMM_WORLD,ier)              
            call MPI_SEND(w(1,j,k),Im,MPI_DOUBLE_PRECISION,0,3
     >                   ,MPI_COMM_WORLD,ier)   
          end if
          call MPI_BARRIER(MPI_COMM_WORLD,ier)
        end do
      end do
      if(Np.eq.0) then
        write(9) comment
        close(9)
      end if

*
333   call MPI_FINALIZE(ier)
      stop
111   write(*,*)'  cp file ',fncp,' was not found'
      istop=1
      goto 222
      end
