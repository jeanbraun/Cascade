c tectonic_movement

      subroutine tectonic_movement (x,y,nnode,
     &                              fix,dt,itime,
     &                              points,vertices,neighbour,
     &                              surface,newsurface,nt,
     &                              nn,nb,nbmax,nswaps,
     &                              mask,mask_e,xy,pp,aa,bb,surfscale,
     &                              sides)

c this routine allows the points of the grid to be moved around
c due to horizontal tectonic movements

c INPUT: x, y      = x- and y-nodal coordinates
c        nnode     = number of nodes
c        fix       = boundary conditions
c        dt        = time step length
c        itime     = time step
c        points    = working array
c        vertices  = connectivity list
c        neighbour = triangle neighbour list
c        surface   = surface attached to each node
c        newsurface= to determine whether nodal surface has to be reestimated
c        nt        = number of triangles
c        nn        = neighbour list
c        nb        = number of neighbour per node
c        nbmax     = maximum number of neighbours
c        mask      = working array
c        mask_e    = working array
c        xy        = working array
c        pp        = working array
c        aa        = working array
c        bb        = working array

c OUTPUT: x,y       = updated nodal coordinates
c         nt        =number of triangles
c        vertices  = connectivity list
c        neighbour = triangle neighbour list
c        surface   = surface attached to each node
c        nn        = neighbour list
c        nb        = number of neighbour per node
c        nswaps    = number of triangle side swaps operated during the
c                    grid reorganisation
c        sides     = length of the sides of the voronoi cells

c subroutines called:
c - debug
c - del_flip
c - find_surface

      common /vocal/ ivocal

      real*4    x(nnode),y(nnode)
      real*8    points(2,*)
      real      fix(nnode)
      integer   vertices(3,*),neighbour(3,*)
      real*4    surface(nnode)
      real      newsurface(nnode)
      integer   nb(nnode),nn(nbmax,nnode)
      real      sides(nbmax,*)
      logical   mask(*),mask_e(3,*)
      real      xy(2),pp(2,nbmax),aa(nbmax,2),bb(nbmax)

c it is in this loop that x and y are updated
c this is the spinning landscape
      omega=2.*3.1415/1.e6
        do i=1,nnode
        x0=x(i)-50.e3
        y0=y(i)-50.e3
          if (x0**2+y0**2.lt.25.e3**2) then
          x(i)=x(i)-omega*dt*y0
          y(i)=y(i)+omega*dt*x0
          endif
        enddo

c from here on you should not change anything...

        do i=1,nnode
        points(1,i)=x(i)
        points(2,i)=y(i)
        enddo

        do i=1,nnode
        newsurface(i)=0.
        enddo

      if (ivocal.eq.1) call debug ('del_flip$',0)
      call del_flip (points,neighbour,vertices,nt,mask,mask_e,nswaps)
        do it=1,nt
        if (mask(it)) then
        newsurface(vertices(1,it))=1.
        newsurface(vertices(2,it))=1.
        newsurface(vertices(3,it))=1.
        endif
        enddo
      if (ivocal.eq.1) call debug ('tectonic_movement$',1)

c if there were any swapping the natural neighbours and surfaces need to be 
c recalculated...

      if (nswaps.ne.0) then

      call find_neighbour_list (nb,nn,nnode,nbmax,x,y,sides,
     &                          nt,vertices,neighbour)

      if (ivocal.eq.1) call debug ('find_surface$',0)
      call find_surface (nn,nb,surface,nbmax,nnode,
     &                   x,y,xy,pp,aa,bb,newsurface,surfscale)
      if (ivocal.eq.1) call debug ('tectonic_movement$',1)

        do i=1,nnode
        nb(i)=nb(i)+1
          if (nb(i).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i),i)=i
        enddo

      endif

      return
      end
