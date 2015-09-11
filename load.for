*
      subroutine load(fncp,u,v,w,t,dt,Dp,Imax,Jmax,Kmax)
      implicit real*8 (a-h,o-z)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr,nsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
      if(Im.ge.Imax)    write(*,*) '(Im.ge.Imax)    '
      if(Jm.ge.Jmax)    write(*,*) '(Jm.ge.Jmax)    '
      if(Km.ge.Kmax)    write(*,*) '(Km.ge.Kmax)    '
*
      open(9,file=fncp,form='unformatted',status='old',err=111)
      read(9)t,dt,Dp,Re,Xmax_,epsr_,lx_,Jm_,lt_,nsym_
*
      write(*,*)t,dt,Dp,Re,Xmax_,epsr_,lx_,Jm_,lt_,nsym_
      if(lt.ne.lt_)     write(*,*) '(lt.ne.lt_)    '
      if(Jm.ne.Jm_)     write(*,*) '(Jm.ne.Jm_)    '
      if(lx.ne.lx_)     write(*,*) '(lx.ne.lx_)    '
      if(nsym.ne.nsym_) write(*,*) '(nsym.ne.nsym_)'
      if(Xmax.ne.Xmax_) write(*,*) '(Xmax.ne.Xmax_)'
      if(epsr.ne.epsr_) write(*,*) '(epsr.ne.epsr_)'
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
111   write(*,*) 'file.scp was not found'
      return
      end
*
