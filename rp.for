*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,or,ot,buf,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ut(0:Imax,0:Jmax,0:*)
     >,vt(0:Imax,0:Jmax,0:*)
     >,wt(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
     >,buf(1:*)
*
      call bc_om(t,u,v,w,ox,or,ot,buf,Imax,Jmax)
      call visc(t,ut,vt,wt,ox,or,ot,Imax,Jmax)
      call add_nl(t,u,v,w,ox,or,ot,ut,vt,wt,Imax,Jmax)
*
      return
      end
