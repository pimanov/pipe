*
*     program an
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      parameter (Imax=513, Jmax=129, Kmax=129)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
     >,buf(2*Imax*Jmax*Kmax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
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
        open(5,file='init.car')
        read(5,*) Xmax, lx
        read(5,*) epsr, Jm
        read(5,*) nsym, lt
        read(5,*) Re
        read(5,*) dt
        read(5,*) a
        read(5,*) fncp
        t=0.0
        close(5)
      end if
*
      call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Re,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Xmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(epsr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(lx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(Jm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(lt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(nsym,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(a,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
*
      call com
      if(Np.eq.0) then      
      write(*,*)' ***************************************************'
      write(*,*)' *        Number of processors =',Npm,'          *'
      write(*,200) t,dt,a,Re,Xmax,epsr,Imm,Jm,Km,nsym
      write(*,*)' ***************************************************'
200   format('    t=',1pe10.3,' dt=',e9.2,' amp=',e9.2,/,
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
          do i=1,Im
            u(i,j,k)=sin(0.1*i+0.01*j+0.3*k+Np*3.0+14.0)
            v(i,j,k)=cos(0.02*i+0.1*j+0.2*k+Np*2.0+22.0)
            w(i,j,k)=sin(0.05*i+0.05*j+0.1*k+Np*0.02+10.0)
          end do
        end do
      end do
      p(0,0,1)=0.d0
      call pres(u,v,w,p,buf,Imax,Jmax)
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
      call MPI_ALLREDUCE(dd,dds,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     >  MPI_COMM_WORLD,ier)
      if(Np.eq.0) write(*,*) 'disturb central div = ',dds

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
      call MPI_ALLREDUCE(amp,amps,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >  MPI_COMM_WORLD,ier)
      call MPI_ALLREDUCE(vol,vols,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >  MPI_COMM_WORLD,ier)
      amp=sqrt(amps/vols)
      if(Np.eq.0) write(*,*) 'disturb amp = ', amp

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
      call MPI_ALLREDUCE(amp,amps,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >  MPI_COMM_WORLD,ier)
      call MPI_ALLREDUCE(vol,vols,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >  MPI_COMM_WORLD,ier)
      amp=sqrt(amps/vols)
      if(Np.eq.0) write(*,*) 'disturb amp = ', amp

      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u(i,j,k)+1.0-yt(j)**2
          end do
        end do
      end do
      p(0,0,1)=0.5d0
      call pres(u,v,w,p,buf,Imax,Jmax)

      dd=0.0
      do k=2,Km-1
        do j=2,Jm-1
          do i=2,Im-1
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      call MPI_ALLREDUCE(dd,dds,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     >  MPI_COMM_WORLD,ier)
      if(Np.eq.0) write(*,*) 'Pois+dist central div = ',dds


      if(Np.eq.0) then
        open(9,file=fncp,form='unformatted')
        Dp=p(0,0,0)
        write(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
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
     >                     ,n,1,MPI_COMM_WORLD,status,ier)
              call MPI_RECV(buf(n*Im+1+Imm),Im,MPI_DOUBLE_PRECISION
     >                     ,n,2,MPI_COMM_WORLD,status,ier)
              call MPI_RECV(buf(n*Im+1+2*Imm),Im,MPI_DOUBLE_PRECISION
     >                     ,n,3,MPI_COMM_WORLD,status,ier)
            end do
            write(9)(buf(i),i=1,Imm)
            write(9)(buf(i+Imm),i=1,Imm)
            write(9)(buf(i+2*Imm),i=1,Imm)
          else
            call MPI_SEND(u(1,j,k),Im,MPI_DOUBLE_PRECISION,0,1
     >                     ,MPI_COMM_WORLD,ier)              
            call MPI_SEND(v(1,j,k),Im,MPI_DOUBLE_PRECISION,0,2
     >                     ,MPI_COMM_WORLD,ier)              
            call MPI_SEND(w(1,j,k),Im,MPI_DOUBLE_PRECISION,0,3
     >                     ,MPI_COMM_WORLD,ier)   
          end if
        end do
      end do
*
333   call MPI_FINALIZE(ier)
      stop
      end
