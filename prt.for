*
      subroutine prt(t,dt,u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
     >/cf/cf
*
      Dp=p(0,0,0)

      Ss=0.d0
      ubulk=0.d0
      volume=0.d0
      enrg=0.d0
      dd=0.d0
      do k=1,Km
        do j=1,Jm
          Ss=Ss+ht*yt1(j)
          ubulk=ubulk+u(1,j,k)*ht*yt1(j)
          do i=1,Im
            uP=0.d0
            volume=volume+ht*yt1(j)*hx
            enrg=enrg+(u(i,j,k)-uP)**2*ht*yt1(j)*hx
     >                    +v(i,j,k)**2*ht*rt1(j)*hx
     >                    +w(i,j,k)**2*ht*yt1(j)*hx
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      ubulk=ubulk/Ss
      enrg=enrg/2

      call MPI_ALLREDUCE(volume,vols,1,MPI_DOUBLE_PRECISION,MPI_SUM
     >               ,MPI_COMM_WORLD,ier)
      volume=vols*2*dsym
      call MPI_ALLREDUCE(enrg,enrgs,1,MPI_DOUBLE_PRECISION,MPI_SUM
     >               ,MPI_COMM_WORLD,ier)
      enrg=enrgs*2*dsym
      call MPI_ALLREDUCE(dd,dds,1,MPI_DOUBLE_PRECISION,MPI_MAX
     >               ,MPI_COMM_WORLD,ier)
      dd=dds

      amp=0.d0
      do j=1,Jm
        u0=0.d0
        do k=1,Km
          do i=1,Im
            u0=u0+u(i,j,k)
          end do
        end do
        u0=u0/(Im*Km)
        call MPI_ALLREDUCE(u0,u0s,1,MPI_DOUBLE_PRECISION,MPI_SUM
     >               ,MPI_COMM_WORLD,ier)
        u0=u0s/Npm

        do k=1,Km
          do i=1,Im
            amp=amp+ht*yt1(j)*hx*(u(i,j,k)-u0)**2
     >             +ht*rt1(j)*hx*v(i,j,k)**2
     >             +ht*yt1(j)*hx*w(i,j,k)**2
          end do
        end do
      end do
      call MPI_ALLREDUCE(amp,amps,1,MPI_DOUBLE_PRECISION,MPI_SUM
     >               ,MPI_COMM_WORLD,ier)
      amp=sqrt(amps*2*dsym)

      ucl=0.d0
      do k=1,Km
        do i=1,Im
          ucl=ucl+u(i,1,k)
        end do
      end do
      ucl=ucl/(Im*Km)
      call MPI_ALLREDUCE(ucl,ucls,1,MPI_DOUBLE_PRECISION,MPI_SUM
     >               ,MPI_COMM_WORLD,ier)
      ucl=ucls/Npm
*
      u1=u(Im/2,Jm/2,Km/2)
      v1=v(Im/4,Jm/4,Km/4)
      if(Np.eq.0)then
        write(8,120)t,dt,amp,enrg,ucl,Dp,cf,ubulk,dd,u1,v1
        write(*,110)t,dt,amp,enrg,ucl,Dp,cf,ubulk,dd,u1,v1
      end if 
120   format(15e25.15)
110   format('t=',f10.4,',dt=',f10.4,',amp=',e12.4,',enr=',e12.4,',Ucl='
     > ,e12.4,',Dp=',e12.4,',cf=',e12.4,',ub=',e12.4,',dd=',e12.4
     > ,',u*=',e12.4,',v*=',e12.4)
      return
      end
