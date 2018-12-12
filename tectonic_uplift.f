c tectonic_uplift

      subroutine tectonic_uplift(x,y,h,h0,hi,nnode,
     &                           fix,dt,time,
     &                           influx,surface)

c this routine defines the tectonic uplift function
c it may vary in space and time

c INPUT: x,y       = x- and y-nodal coordinates
c        h         = present topography
c        h0        = bedrock-alluvials interface
c        hi        = original topograpy
c        nnode     = number of nodes
c        fix       = boundary conditions
c        dt        = time step length
c        time      = current time
c        surface   = surface associated with each node

c OUTPUT:  h        = update topography
c          hi       = updated original topograpy
c          h0       = updated bedrock-alluvials interface
c          influx   = updated influx of material by tectonic uplift

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     x(nnode),y(nnode),h(nnode),h0(nnode),hi(nnode)
      real     fix(nnode)
      real     surface(nnode),influx

c the example is a square topography that stops growing after a while
c note that we have used the array fix to prevent the base level to move
c with tectonic uplift

!      if (time.gt.1.e6) uplift_rate=0.
!      uplift_rate=0.
        do inode=1,nnode
        uplift_rate=0.
        if (abs(x(inode)-125.e3).lt.25.e3 .and. 
     &      y(inode).lt.125.e3) uplift_rate=1.e-2
        dh=uplift_rate*dt*fix(inode)
        h(inode)=h(inode)+dh
c do not touch the following lines
c they update h0, hi and calculate the influx of material into the
c landscape
        h0(inode)=h0(inode)+dh
        hi(inode)=hi(inode)+dh
        influx=influx+dh*surface(inode)
        enddo

      return
      end
