*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,or,ot,buf,Imax,Jmax)
*
      call bc_om(t,u,v,w,ox,or,ot,buf,Imax,Jmax)
      call visc(t,ox,or,ot,ut,vt,wt,Imax,Jmax)
      call add_nl(t,u,v,w,ox,or,ot,ut,vt,wt,Imax,Jmax)
*
      return
      end
