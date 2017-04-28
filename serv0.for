*
      subroutine servis(t,g1,g2,g3,g4,g5,
     >u,v,w,ox,or,ot,p,ii,buf,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      parameter (Imx=257, Jmx=129, Kmx=129)
      character*24 fngcp1, fngcp2
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,g1(0:Imax,0:Jmax,0:*)
     >,g2(0:Imax,0:Jmax,0:*)
     >,g3(0:Imax,0:Jmax,0:*)
     >,g4(0:Imax,0:Jmax,0:*)
     >,g5(0:Imax,0:Jmax,0:*)
     >,buf(*)
     >,wt(0:Imx,0:Jmx,0:Kmx)
     >,wt1(0:Imx,0:Jmx,0:Kmx)
     >,wt2(0:Imx,0:Jmx,0:Kmx)
     >,vt(0:Imx,0:Jmx,0:Kmx)
     >,vt1(0:Imx,0:Jmx,0:Kmx)
     >,vt2(0:Imx,0:Jmx,0:Kmx)
     >,ox1(0:Imx,0:Jmx,0:Kmx)
     >,ox2(0:Imx,0:Jmx,0:Kmx)
     >,visc(0:Imx,0:Jmx,0:Kmx)
      common
*     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/servst/iserv
     >/servsn/nserv
     >/cf/cf
*
      fngcp1="OXgen1.scp"
      fngcp2="OXgen2.scp"
*
      if (iserv.eq.0) then
        iserv = 1
        nserv = 0
        do k=1,Km
          do j=1,Jm
            do i=1,Im
              g1(i,j,k)=0.d0
              g2(i,j,k)=0.d0
              g3(i,j,k)=0.d0
              g4(i,j,k)=0.d0
              g5(i,j,k)=0.d0
            end do
          end do
        end do
        return
      end if
*
      if(ii.eq.1) then
        call get_nl_part(u,v,w,ox,or,ot,vt1,vt2,wt1,wt2,Imax,Jmax)
        call get_rotx(vt1,wt1,ox1,Imax,Jmax)   
        call get_rotx(vt2,wt2,ox2,Imax,Jmax)   
        call get_visc(ox,or,ot,vt,wt,Imax,Jmax)
        call get_rotx(vt,wt,visc,Imax,Jmax)   
        do k=1,Km
          do j=1,Jm
            do i=1,Im
      u_fnn=(u(i-1,j+1,k+1)+u(i-1,j+1,k)+u(i-1,j,k+1)+u(i-1,j,k)+
     >u(i,j+1,k+1)+u(i,j+1,k)+u(i,j,k+1)+u(i,j,k))/8
      oxx_fnn=(ox(i+1,j,k)-ox(i-1,j,k))/(2*hx)
      ox_fnn=ox(i,j,k)
      ux_fnn=((u(i,j+1,k+1)+u(i,j+1,k)+u(i,j,k+1)+u(i,j,k))- 
     >(u(i-1,j+1,k+1)+u(i-1,j+1,k)+u(i-1,j,k+1)+u(i-1,j,k)))/(4*hx)
              cux=-u_fnn*oxx_fnn
              d1x=ox_fnn*ux_fnn

              dtx=ox1(i,j,k)-cux
              ctx=ox2(i,j,k)-d1x

              g1(i,j,k)=g1(i,j,k)+cux
              g2(i,j,k)=g2(i,j,k)+d1x
              g3(i,j,k)=g3(i,j,k)+ctx
              g4(i,j,k)=g4(i,j,k)+dtx
              g5(i,j,k)=g5(i,j,k)+visc(i,j,k)
            end do
          end do
        end do
        nserv=nserv+1
      end if
*
      if(ii.eq.2) then
        do k=1,Km
          do j=1,Jm
            do i=1,Im
              g1(i,j,k)=g1(i,j,k)/nserv
              g2(i,j,k)=g2(i,j,k)/nserv
              g3(i,j,k)=g3(i,j,k)/nserv
              g4(i,j,k)=g4(i,j,k)/nserv
              g5(i,j,k)=g5(i,j,k)/nserv
            end do
          end do
        end do
        dt=0.25
        Dp=0.0
        call write_cp(t,dt,Dp,g1,g2,g3,fngcp1,buf,Imax,Jmax)
        call write_cp(t,dt,Dp,g4,g5,w,fngcp2,buf,Imax,Jmax)
*
      end if
      return
      end
