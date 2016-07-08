*
      subroutine read_scp(fncp,u,v,w,t,dt,Dp,Imax,Jmax,Kmax)
      implicit real*8 (a-h,o-z)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
      open(9,file=fncp,form='unformatted',status='old',err=111)
      read(9)t,dt,Dp,Re,Xmax1,epsr1,lx1,Jm1,lt1,nsym
*
      write(*,*)t,dt,Dp,Re,Xmax1,epsr1,lx1,Jm1,lt1,nsym
      if(lt.ne.lt1)     write(*,*) '(lt.ne.lt1)    '
      if(Jm.ne.Jm1)     write(*,*) '(Jm.ne.Jm1)    '
      if(lx.ne.lx1)     write(*,*) '(lx.ne.lx1)    '
      if(dsym.ne.nsym) write(*,*) '(dsym.ne.nsym)'
      if(Xmax.ne.Xmax1) write(*,*) '(Xmax.ne.Xmax1)'
      if(epsr.ne.epsr1) write(*,*) '(epsr.ne.epsr1)'
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
*
      subroutine read_dcp(fncp,u,v,w,t,dt,Dp,Imax,Jmax,Kmax)
      implicit real*8 (a-h,o-z)
      character*24 fncp
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,lx
     >/dimr/rt(0:129),rt1(0:129),yt(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
*
      if(Im.ge.Imax)    write(*,*) '(Im.ge.Imax)    '
      if(Jm.ge.Jmax)    write(*,*) '(Jm.ge.Jmax)    '
      if(Km.ge.Kmax)    write(*,*) '(Km.ge.Kmax)    '
*
      open(9,file=fncp,form='unformatted',status='old',err=111)
      read(9)t,dt,Dp,Re,Xmax1,epsr1,lx1,Jm1,lt1,dsym1
*
      write(*,*)t,dt,Dp,Re,Xmax1,epsr1,lx1,Jm1,lt1,dsym1
      if(lt.ne.lt1)     write(*,*) '(lt.ne.lt1)    '
      if(Jm.ne.Jm1)     write(*,*) '(Jm.ne.Jm1)    '
      if(lx.ne.lx1)     write(*,*) '(lx.ne.lx1)    '
      if(dsym.ne.dsym1) write(*,*) '(dsym.ne.dsym1)'
      if(Xmax.ne.Xmax1) write(*,*) '(Xmax.ne.Xmax1)'
      if(epsr.ne.epsr1) write(*,*) '(epsr.ne.epsr1)'
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
