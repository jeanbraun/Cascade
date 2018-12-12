c diffusion_erosion

      subroutine diffusion_erosion (xkdiff,
     &                              x,y,h,dh,nnode,
     &                              nb,nn,nbmax,
     &                              fix,dt,dhmin,dhmax,taumin,
     &                              sea_level,outflux,surface,sides)

c subroutine to calculate diffusion erosion

c INPUT: xkdiff     = nodal diffusivity
c        x, y       = x- and y- nodal coordinates
c        h          = topography
c        dh         = total height removed/added
c        nnode      = number of nodes
c        nb         = number of natural neighbours per node
c        nn         = list of natural neighbours
c        fix        = boundary condition array
c        dt         = time step length (in yrs)
c        dhmin      = minimum amount eroded/deposited
c        dhmax      = maximum amount eroded/deposited
c        sea_level  = sea-level in meters
c        outflux    = flux out of the system through base level
c        surface    = array containing the surface attached to each node
c        sides       = array containing the lengths of the sides of the
c                     Voronoi cells

c OUTPUT: a range of arrays/variables are updated:
c        h, dh, dhmin, dhmax, outflux

c subroutines called:
c - none

      common /vocal/ ivocal

      real    x(nnode),y(nnode),h(nnode),dh(nnode)
      real    xkdiff(nnode)
      real    fix(nnode)
      real    surface(nnode)
      integer nb(nnode),nn(nbmax,nnode)
      real    sides(nbmax,nnode)

c finds height change

      taumin=1.e6
        do i=1,nnode
c          if (fix(i).gt.0.5) then
          dhh=0.
          xk=xkdiff(i)
          surf=surface(i)
          taup=1.e6
          hc=h(i)
            do j=1,nb(i)
            ic=nn(j,i)
              if (ic.ne.i) then
              xl=sqrt((x(i)-x(ic))**2+(y(i)-y(ic))**2)
              side=sides(j,i)
              slope=(hc-h(ic))/xl
              dhh=dhh+slope*side
              taup=min(taup,xl/side)
              endif
            enddo
          dhh=-dhh*xk*dt/surf
          taumin=min(taumin,taup*surf/xk)
            if (fix(i).lt.0.5) then
            outflux=outflux+dhh*surf
            dhh=0.
            endif
            if (h(i)+dhh.lt.sea_level) then
            dhh=sea_level-h(i)
            endif
          dh(i)=dhh
          h(i)=h(i)+dhh
          dhmin=min(dhmin,dhh)
          dhmax=max(dhmax,dhh)
c          endif
        enddo

      return
      end
