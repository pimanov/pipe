*
      subroutine base_rp(t,u,v,w,ox,or,ot,buf,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
     >,buf(1:*)
      integer nstat(MPI_STATUS_SIZE,12),nreq(12)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
*
      t=t
* Boundary conditions
      Np0=Np-1
      if(Np.eq.0)Np0=Npm-1
      Np1=Np+1
      if(Np.eq.Npm-1)Np1=0
      lng=Jm*Km
      l=0
      do k=1,Km
        do j=1,Jm
          l=l+1
          buf(l)=u(Im,j,k)
          buf(l+lng)=v(1,j,k)
          buf(l+2*lng)=v(Im,j,k)
          buf(l+3*lng)=w(1,j,k)
          buf(l+4*lng)=w(Im,j,k)
          buf(l+10*lng)=u(1,j,k)
        end do
      end do
      call MPI_ISEND(buf(1),lng,MPI_DOUBLE_PRECISION,Np1,1
     >              ,MPI_COMM_WORLD,nreq(1),ier)
      call MPI_ISEND(buf(1+lng),lng,MPI_DOUBLE_PRECISION,Np0,2
     >              ,MPI_COMM_WORLD,nreq(2),ier)
      call MPI_ISEND(buf(1+2*lng),lng,MPI_DOUBLE_PRECISION,Np1,3
     >              ,MPI_COMM_WORLD,nreq(3),ier)
      call MPI_ISEND(buf(1+3*lng),lng,MPI_DOUBLE_PRECISION,Np0,4
     >              ,MPI_COMM_WORLD,nreq(4),ier)
      call MPI_ISEND(buf(1+4*lng),lng,MPI_DOUBLE_PRECISION,Np1,5
     >              ,MPI_COMM_WORLD,nreq(5),ier)
      call MPI_ISEND(buf(1+10*lng),lng,MPI_DOUBLE_PRECISION,Np0,6
     >              ,MPI_COMM_WORLD,nreq(11),ier)
*
      call MPI_IRECV(buf(1+5*lng),lng,MPI_DOUBLE_PRECISION,Np0,1
     >              ,MPI_COMM_WORLD,nreq(6),ier)
      call MPI_IRECV(buf(1+6*lng),lng,MPI_DOUBLE_PRECISION,Np1,2
     >              ,MPI_COMM_WORLD,nreq(7),ier)
      call MPI_IRECV(buf(1+7*lng),lng,MPI_DOUBLE_PRECISION,Np0,3
     >              ,MPI_COMM_WORLD,nreq(8),ier)
      call MPI_IRECV(buf(1+8*lng),lng,MPI_DOUBLE_PRECISION,Np1,4
     >              ,MPI_COMM_WORLD,nreq(9),ier)
      call MPI_IRECV(buf(1+9*lng),lng,MPI_DOUBLE_PRECISION,Np0,5
     >              ,MPI_COMM_WORLD,nreq(10),ier)
      call MPI_IRECV(buf(1+11*lng),lng,MPI_DOUBLE_PRECISION,Np1,6
     >              ,MPI_COMM_WORLD,nreq(12),ier)
*
      call MPI_WAITALL(12,nreq,nstat,ier)  
*
      l=0
      do k=1,Km
        do j=1,Jm
          l=l+1
          u(0,j,k)=buf(l+5*lng)
          v(Im+1,j,k)=buf(l+6*lng)
          v(0,j,k)=buf(l+7*lng)
          w(Im+1,j,k)=buf(l+8*lng)
          w(0,j,k)=buf(l+9*lng)
          u(Im+1,j,k)=buf(l+11*lng)
        end do
      end do
*
      do i=0,Im
        do j=1,Jm
          u(i,j,0)=u(i,j,1)
          u(i,j,Km+1)=u(i,j,Km)
        end do
        do k=1,Km
          u(i,Jm+1,k)=-u(i,Jm,k)
        end do
      end do
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=0.d0
        end do
      end do
      do j=1,Jm
        do i=1,Im
          v(i,j,0)=v(i,j,1)
          v(i,j,Km+1)=v(i,j,Km)
        end do
      end do
      do j=1,Jm
        do i=0,Im+1
          w(i,j,0)=0.d0
          w(i,j,Km)=0.d0
        end do
      end do
      do k=0,Km
        do i=1,Im
          w(i,Jm+1,k)=-w(i,Jm,k)*yt(Jm)/yt(Jm+1)
        end do
      end do
*
* Vorticities
      do i=1,Im
        do j=1,Jm
          do k=0,Km
            w0=w(i,j,k)
            w1=w(i,j+1,k)
            v0=v(i,j,k)
            v1=v(i,j,k+1)
            ox(i,j,k)=((yt(j+1)*w1-yt(j)*w0)/rt1(j)
     >                -(v1-v0)/ht)/rt(j)
          end do
        end do
        j=0
          do k=0,Km
            ox(i,j,k)=0.d0
          end do
      end do
      do k=0,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j,k+1)
            w0=w(i,j,k)
            w1=w(i+1,j,k)
            or(i,j,k)=(u1-u0)/(yt(j)*ht)
     >               -(w1-w0)/hx
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j+1,k)
            v0=v(i,j,k)
            v1=v(i+1,j,k)
            ot(i,j,k)=(v1-v0)/hx
     >               -(u1-u0)/rt1(j)
          end do
        end do
      end do
*
      return
      end
