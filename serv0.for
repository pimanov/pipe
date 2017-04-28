*
      subroutine servis(t,g1,g2,g3,g4,u,v,w,ox,or,ot,p,ii,buf,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      character*24 fngcp
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
     >,buf(*)
     >,wt(0:257,0:41,0:33),wt1(0:257,0:41,0:33),wt2(0:257,0:41,0:33)
     >,vt(0:257,0:41,0:33),vt1(0:257,0:41,0:33),vt2(0:257,0:41,0:33)
     >,ox1(0:257,0:41,0:33),ox2(0:257,0:41,0:33),visc(0:257,0:41,0:33)
      common
*     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/servst/iserv
     >/servsn/nserv
     >/cf/cf
*
      fngcp="OXgen1.scp"
      fngcp="OXgen2.scp"
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
            end do
          end do
        end do
        return
      end if
*
      if(ii.eq.1) then
        call get_nl_part(u,v,w,ox,or,ot,vt1,vt2,wt1,wt2,Imax,Jmax)
        return
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
              cxx=cux+d1x

              g1(i,j,k)=g1(i,j,k)+cxx
              g2(i,j,k)=g2(i,j,k)+ctx
              g3(i,j,k)=g3(i,j,k)+dtx
              g4(i,j,k)=g4(i,j,k)+visc(i,j,k)
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
            end do
          end do
        end do
        dt=0.25
        Dp=0.0
        call write_cp(t,dt,Dp,g1,g2,g3,fngcp2,buf,Imax,Jmax)
        call write_cp(t,dt,Dp,g4,v,w,fngcp1,buf,Imax,Jmax)
*
      end if
      return
      end
