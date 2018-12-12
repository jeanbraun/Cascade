c initialize_nodal_geometry

      subroutine initialize_nodal_geometry
     & (nnodemax,nnode,x,y,h,fix,delta,surfscale,
     &  run_name,nrun_name)

c subroutine to read in initial nodal information

c INPUT: nnodemax      = maximum number of nodes that can be read
c        run_name      = name of the subdirectory where the run output
c                        will be stored
c        nrun_name     = length of the character string run_name

c OUTPUT: nnode        = number of nodes read
c         x(nnode)     = x-location of nodes (in m)
c         y(nnode)     = y-location of nodes (in m)
c         h(nnode)     = height of nodes (in m)
c         fix(nnode)   = boundary condition (=0 means node height is fixed;
c                        =1 means node height is free)
c         delta        = mean nodal spacing (in m)
c         surfscale    = mean nodal surface (a nodal surface is the surface
c                        of the part of the landscape associated with a node)
c                        (in m**2)

c subroutines called:
c - debug
c - random
c - iread_but_skip_comment
c - read_but_skip_comment

      common /vocal/ ivocal

      real          x(*),y(*),h(*)
      real          fix(*)
      character     run_name*256

      iread=0

      if (iread.eq.0) then

      nx=128
      ny=128
      nnode=nx*ny
        if (nnode.gt.nnodemax) then
        stop 'Too many nodes...'
        endif

      sidex=250.e3
      sidey=250.e3

      delta=sidex/float(nx-1)
      surfscale=sidex*sidey/nnode

c this is a bit of random noise put on the initial grid
c note that the noise is not used for the nodes along the boundary

      if (ivocal.eq.1) call debug ('random$',0)
      call random (x,nnode)
      call random (y,nnode)
      call random (h,nnode)
      if (ivocal.eq.1) call debug ('initialize_nodal_geometry$',1)

      hinit=1. 
        do j=1,ny
        do i=1,nx
        ishake=1
        if (i.eq.1 .or.i.eq.nx .or.j.eq.1 .or.j.eq.ny) ishake=0
        inode=(j-1)*nx+i
        x(inode)=((x(inode)-.5)/float(nx-1)*ishake
     &            +float(i-1)/float(nx-1))*sidex
        y(inode)=((y(inode)-.5)/float(ny-1)*ishake
     &            +float(j-1)/float(ny-1))*sidey
	fix(inode)=1.
        if (j.eq.1.or.j.eq.ny) fix(inode)=0.
        if (i.eq.1.or.i.eq.nx) fix(inode)=0.
!        if (j.eq.ny) fix(inode)=0.
!        dh=y(inode)/sidey*2.*hinit
!        if (y(inode).gt.sidey/2.) dh=2*hinit-dh
!        h(inode)=h(inode)/1.e3+dh
        r=(x(inode)-sidex)**2+(y(inode)-sidey)**2
        h(inode)=hinit*exp(-r/(sidex/2.)**2)*fix(inode)
!        if (abs(x(inode)-500.e3).lt.400.e3 .and. 
!     &      abs(y(inode)-500.e3).lt.400.e3) h(inode)=h(inode)+1000.
        enddo
        enddo
 
      else

      open (7,file=run_name(1:nrun_name)//
     &      '/cascade.node.in',status='old')

      call iread_but_skip_comment (7,1,nnode)
        if (nnode.gt.nnodemax) then
        stop 'Too many nodes...'
        endif

      call read_but_skip_comment (7,nnode,x)
      call read_but_skip_comment (7,nnode,y)
      call read_but_skip_comment (7,nnode,h)
      call read_but_skip_comment (7,nnode,fix)
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
        do i=1,nnode
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
        ymin=min(ymin,y(i))
        ymax=max(ymax,y(i))
        enddo
      sidex=xmax-xmin
      sidey=ymax-ymin
      snode=sqrt(float(nnode))
      delta=(sidex/snode+sidey/snode)/2.
      surfscale=sidex*sidey/nnode

      close (7)

      endif

        if (nnode.lt.3) then
        stop 'nnode too small...'
        endif

        if (delta.le.0.) then
        stop 'delta must be greater than 0...'
        endif

        if (surfscale.le.0.) then
        stop 'surfscale must be greater than 0...'
        endif

        if (sidex.le.0.) then
        stop 'sidex must be greater than 0...'
        endif

        if (sidey.le.0.) then
        stop 'sidey must be greater than 0...'
        endif

      return
      end
