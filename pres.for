*
      subroutine pres(u,v,w,p,buf,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,buf(1:*)
     >,a1(2048),a2(2048)
     >,b1(2048),b2(2048)
     >,bp(2048),dp(2048),ep(2048)
      integer nstat(MPI_STATUS_SIZE,2048),nreq(2048)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Kmm,lt
     >/rlx/rlx(2048)
     >/rlt/rlt(256)
     >/pry/apy(128),bpy(128),cpy(128)
     >/proc/Np,Npm
*
      Km=Kmm/Npm
*
* Boundary conditions
      Np0=Np-1
      if(Np.eq.0)Np0=Npm-1
      Np1=Np+1
      if(Np.eq.Npm-1)Np1=0
      l=0
      do k=1,Kmm
        do j=1,Jm
          l=l+1
          buf(l)=u(Im,j,k)
        end do
      end do
      lng=Jm*Kmm
      call MPI_ISEND(buf(1),lng,MPI_DOUBLE_PRECISION,Np1,1
     >              ,MPI_COMM_WORLD,nreq(1),ier)
      call MPI_IRECV(buf(1+lng),lng,MPI_DOUBLE_PRECISION,Np0,1
     >              ,MPI_COMM_WORLD,nreq(2),ier)
      call MPI_WAITALL(2,nreq,nstat,ier)  
      l=lng
      do k=1,Kmm
        do j=1,Jm
          l=l+1
          u(0,j,k)=buf(l)
        end do
      end do
*
      Im2=Imm/2
      cik=4.d0/(Im2*Kmm)
*
      do i=1,Im
        do j=1,Jm
          w(i,j,0)=0.d0
          w(i,j,Kmm)=0.d0
        end do
        do k=1,Kmm
          v(i,Jm,k)=0.d0
          do j=1,Jm
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            p(i,j,k)=d
          end do
        end do
      end do
*
*   FFT in tt-direction
      do j=1,Jm
        do i=1,Im 
          do k=1,Kmm
            a1(k)=p(i,j,k)
          end do
          call ftc05d(a1,b1,lt)
          do k=1,Kmm
            p(i,j,k)=b1(k)
          end do
        end do
      end do
*
*  Data transfer
      l=0
      do k=1,Kmm
        do j=1,Jm
          do i=1,Im
            l=l+1
            buf(l)=p(i,j,k)
          end do
        end do
      end do
      lng=Im*Jm*Kmm
      length=Im*Jm*Km
      n1=0
      do n=0,Npm-1
        if(n.eq.Np) then           
          l=n*length
          do k=1,Km
            do j=1,Jm
              do i=1,Im
                l=l+1
                buf(lng+l)=buf(l)
              end do
            end do
          end do
        else
          n1=n1+1
          call MPI_ISEND(buf(1+n*length),length,MPI_DOUBLE_PRECISION
     >              ,n,1,MPI_COMM_WORLD,nreq(n1),ier)
          n1=n1+1
          call MPI_IRECV(buf(1+n*length+lng),length,MPI_DOUBLE_PRECISION
     >              ,n,1,MPI_COMM_WORLD,nreq(n1),ier)
        end if
      end do
      call MPI_WAITALL(n1,nreq,nstat,ier)  
      l=lng
      do k=1,Kmm
        do j=1,Jm
          do i=1,Im
            l=l+1
            p(i,j,k)=buf(l)
          end do
        end do
      end do
*
*  FFT in x-direction
      do j=1,Jm
        do k=1,Km
          ii=0
          do n=0,Npm-1
            kk=n*Km+k
            do i=1,Im
              ii=ii+1
              b1(ii)=p(i,j,kk)
            end do
          end do 
          do i=1,Im2
            i1=Imm+1-i
            a1(i)=0.5d0*(b1(i)+b1(i1))
            a2(i)=0.5d0*(b1(i)-b1(i1))
          end do
          call ftc05d(a1,b1,lx-1)
          call fts05d(a2,b2,lx-1)
          do i=1,Im2
            i1=Imm+1-i
            a1(i)=b1(i)*cik
            a1(i1)=b2(i)*cik
          end do
          ii=0
          do n=0,Npm-1
            kk=n*Km+k
            do i=1,Im
              ii=ii+1
              p(i,j,kk)=a1(ii)
            end do
          end do           
        end do
      end do
*
*   Solution in wall-normal coordinate
      k0=Np*Km
      do n=0,Npm-1
        do k=1,Km
          k1=n*Km+k
          kl=k0+k
          do i=1,Im
            il=n*Im+i
            do j=1,Jm
              bp(j)=bpy(j)-rlx(il)-rlt(kl)/yt(j)**2
              dp(j)=p(i,j,k1)
            end do
            bp(1)=bp(1)+apy(1)
            bp(Jm)=bp(Jm)+cpy(Jm)
            Jm1=Jm
            ep(Jm)=0.d0
            if(rlx(il).eq.0.d0.and.rlt(kl).eq.0.d0)Jm1=Jm-1
            call prog3(apy,bp,cpy,dp,ep,Jm1)
            do j=1,Jm
              p(i,j,k1)=ep(j)
            end do
          end do 
        end do
      end do
*
*
*  Inverse FFT in x-direction
      do j=1,Jm
        do k=1,Km
          ii=0
          do n=0,Npm-1
            kk=n*Km+k
            do i=1,Im
              ii=ii+1
              b1(ii)=p(i,j,kk)
            end do
          end do 
          do i=1,Im2
            i1=Imm+1-i
            a1(i)=b1(i)
            a2(i)=b1(i1)
          end do
          call ftc05b(a1,b1,lx-1)
          call fts05b(a2,b2,lx-1)
          do i=1,Im2
            i1=Imm+1-i
            a1(i)=b1(i)+b2(i)
            a1(i1)=b1(i)-b2(i)
          end do
          ii=0
          do n=0,Npm-1
            kk=n*Km+k
            do i=1,Im
              ii=ii+1
              p(i,j,kk)=a1(ii)
            end do
          end do           
        end do
      end do
*
*  Inverse Data transfer
      l=0
      do k=1,Kmm
        do j=1,Jm
          do i=1,Im
            l=l+1
            buf(l)=p(i,j,k)
          end do
        end do
      end do
      lng=Im*Jm*Kmm
      length=Im*Jm*Km
      n1=0
      do n=0,Npm-1
        if(n.eq.Np) then           
          l=n*length
          do k=1,Km
            do j=1,Jm
              do i=1,Im
                l=l+1
                buf(l+lng)=buf(l)
              end do
            end do
          end do
        else
          n1=n1+1
          call MPI_ISEND(buf(1+n*length),length,MPI_DOUBLE_PRECISION
     >              ,n,1,MPI_COMM_WORLD,nreq(n1),ier)
          n1=n1+1
          call MPI_IRECV(buf(1+n*length+lng),length,MPI_DOUBLE_PRECISION
     >              ,n,1,MPI_COMM_WORLD,nreq(n1),ier)
        end if
      end do
      call MPI_WAITALL(n1,nreq,nstat,ier)  
      l=lng
      do k=1,Kmm
        do j=1,Jm
          do i=1,Im
            l=l+1
            p(i,j,k)=buf(l)
          end do
        end do
      end do
*
*   Inverse FFT in tt-direction
      do j=1,Jm
        do i=1,Im 
          do k=1,Kmm
            a1(k)=p(i,j,k)
          end do
          call ftc05b(a1,b1,lt)
          do k=1,Kmm
            p(i,j,k)=b1(k)
          end do
        end do
      end do
*
* Boundary conditions
      Np0=Np-1
      if(Np.eq.0)Np0=Npm-1
      Np1=Np+1
      if(Np.eq.Npm-1)Np1=0
      l=0
      do k=1,Kmm
        do j=1,Jm
          l=l+1 
          buf(l)=p(1,j,k)
        end do
      end do
      lng=Jm*Kmm
      call MPI_ISEND(buf(1),lng,MPI_DOUBLE_PRECISION,Np0,1
     >              ,MPI_COMM_WORLD,nreq(1),ier)
      call MPI_IRECV(buf(1+lng),lng,MPI_DOUBLE_PRECISION,Np1,1
     >              ,MPI_COMM_WORLD,nreq(2),ier)
      call MPI_WAITALL(2,nreq,nstat,ier)  
      l=lng
      do k=1,Kmm
        do j=1,Jm
          l=l+1
          p(Im+1,j,k)=buf(l)
        end do
      end do
*
*  Mean pressure gradient
      Ub=p(0,0,1)*(0.5+curv/6.d0) 
      ss=0.d0
      su=0.d0
      do j=1,Jm
        ss=ss+yt(j)*yt1(j)
        ssu=0.d0
        do k=1,Kmm
          do i=1,Im
            ssu=ssu+u(i,j,k)
          end do
        end do
        su=su+ssu*yt(j)*yt1(j)
      end do
      call MPI_ALLREDUCE(su,sus,1,MPI_DOUBLE_PRECISION,MPI_SUM
     >               ,MPI_COMM_WORLD,ier)
      Dpp=Ub-sus/(Imm*Kmm*ss)
      p(0,0,0)=Dpp
*
      call gradp(u,v,w,p,Imax,Jmax)
      return
      end
*
      subroutine gradp(u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
      Dp=p(0,0,0)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u(i,j,k)-(p(i+1,j,k)-p(i,j,k))/hx+Dp
          end do
        end do
      end do
*
      do k=1,Km
        do j=1,Jm-1
          do i=1,Im
            v(i,j,k)=v(i,j,k)-(p(i,j+1,k)-p(i,j,k))/rt1(j)
          end do
        end do
      end do
*
      do j=1,Jm
        do i=1,Im
          do k=1,Km-1
            w(i,j,k)=w(i,j,k)-(p(i,j,k+1)-p(i,j,k))/(yt(j)*ht)
          end do
        end do
      end do
      return
      end
*
      subroutine div(i,j,k,u,v,w,d,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
      d=(u(i,j,k)-u(i-1,j,k))/hx
     > +(rt(j)*v(i,j,k)-rt(j-1)*v(i,j-1,k))/(yt(j)*yt1(j))
     > +(w(i,j,k)-w(i,j,k-1))/(yt(j)*ht)
      return
      end
