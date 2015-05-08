*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Imax=2049, Jmax=65, Kmax=129)
      character*12 fnscp,fnscs
      dimension
     > u(0:Jmax,0:Kmax)
     >,v(0:Jmax,0:Kmax)
     >,w(0:Jmax,0:Kmax)
     >,buf(Imax)
*
      open(5,file='base_init.car')
      read(5,*) ics
      read(5,*) fnscp
      read(5,*) fnscs
      write(*,*)' ics=  ',isc
      write(*,*)' fnscp=',fnscp
      write(*,*)' fnscs=',fnscs
*
      open(9,file=fnscp,form='unformatted',status='old',err=1)
      read(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym
      Im=2**lx
      Km=2**lt
      do k=1,Km
        do j=1,Jm
          read(9)(buf(i),i=1,Im)
            u(j,k)=buf(ics)
          read(9)(buf(i),i=1,Im)
            v(j,k)=buf(ics)
          read(9)(buf(i),i=1,Im)
            w(j,k)=buf(ics)
        end do
      end do
      close(9)

      if(ics.gt.Im)then
        write(*,*) 'ics>Im, ics=',ics,', Im=',Im
        goto 2
      endif 

      open(9,file=fnscs,form='unformatted')
      write(9)epsr,Jm,lt,nsym
      do k=1,Km
        write(9)(u(j,k),j=1,Jm)
        write(9)(v(j,k),j=1,Jm)
        write(9)(w(j,k),j=1,Jm)
      end do
      write(9) 'base_init.out: Re=',Re,', Xmax=',Xmax,', Im=',Im,
     > 'ics=',ics,', fnscp: ',fnscp
      close(9)
      
      stop
1     write(*,*)'  File ',fncp,' already exists'
2     stop
      end
*
