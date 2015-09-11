*
      subroutine load(fncp,u,v,w,t,dt,Dp,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
      open(9,file=fncp,form='unformatted',status='old',err=111)
      read(9)t,dt,Dp,Re,Xmax_,epsr_,lx_,Jm_,lt_,nsym_
*
      if(lt.ne.lt_)     stop '(lt.ne.lt_)    '
      if(Jm.ne.Jm_)     stop '(Jm.ne.Jm_)    '
      if(lx.ge.lx_)     stop '(lx.ge.lx_)    '
      if(nsym.gt.nsym_) stop '(nsym.gt.nsym_)'
      if(Xmax.ge.Xmax_) stop '(Xmax.ge.Xmax_)'
      if(epsr.gt.epsr_) stop '(epsr.gt.epsr_)'
*
      do k=1,Km
        do j=1,Jm
          read(9)(u(i,j,k),i=1,Im)
          read(9)(v(i,j,k),i=1,Im)
          read(9)(w(i,j,k),i=1,Im)
        end do
      end do
      close(9)
      return
111   stop 'file.scp was not found'
      end
*
