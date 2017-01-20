*
      subroutine servis(t,OXgen,u,v,w,ox,or,ot,p,ii,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      character*24   fngcp
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,OXgen(0:Imax,0:Jmax,0:*)
      common
*     >/dim/Xmax,epsr
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/servst/iserv
     >/servsn/nserv
*
      fngcp = "OXgen.scp"
*
      if (iserv.eq.0) then
        iserv = 1
        nserv = 0
        do k=1,Km
          do j=1,Jm
            do i=1,Im
              OXgen(i,j,k) = 0.d0
            end do
          end do
        end do
      end if
*
      if(ii.eq.0) then
        do k=1,Km
          do j=1,Jm
            do i=1,Im
      u_fnn=(u(i-1,j+1,k+1)+u(i-1,j+1,k)+u(i-1,j,k+1)+u(i-1,j,k)+
     >u(i,j+1,k+1)+u(i,j+1,k)+u(i,j,k+1)+u(i,j,k))/8
      oxx_fnn=(ox(i+1,j,k)-ox(i-1,j,k))/(2*hx)
      ox_fnn=ox(i,j,k)
      ux_fnn=((u(i,j+1,k+1)+u(i,j+1,k)+u(i,j,k+1)+u(i,j,k))- 
     >(u(i-1,j+1,k+1)+u(i-1,j+1,k)+u(i-1,j,k+1)+u(i-1,j,k)))/(4*hx)
      OXgen(i,j,k)=OXgen(i,j,k)-u_fnn*oxx_fnn+ox_fnn*ux_fnn
            end do
          end do
        end do
        nserv = nserv + 1
      end if
*
      if(ii.eq.2) then
        do k=1,Km
          do j=1,Jm
            do i=1,Im
              OXgen(i,j,k) = OXgen(i,j,k) / nserv
            end do
          end do
        end do
        call write_cp(0.0,0.0,0.0,OXgen,OXgen,OXgen,fngcp,buf,Imax,Jmax)
      end if
      return
      end
