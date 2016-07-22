*
      subroutine rp(t,ub,vb,wb,bx,br,bt,
     > u,v,w,ut,vt,wt,ox,or,ot,buf,Imax,Jmax)
*
      call bc_om(t,u,v,w,ox,or,ot,buf,Imax,Jmax)
      call visc(t,ox,or,ot,ut,vt,wt,Imax,Jmax)
      call add_nl(t,ub,vb,wb,ox,or,ot,ut,vt,wt,Imax,Jmax)
      call add_nl(t,u,v,w,bx,br,bt,ut,vt,wt,Imax,Jmax)
*
      return
      end
