c find_neighbours

      subroutine find_neighbours (x,y,nn,nb,nnode,nbmax,
     &                            points,vertices,neighbour,nodes,
     &                            vis_tlist,vis_elist,add_tlist,nt,
     &                            surface,newsurface,eps,
     &                            xx,pp,aa,bb,surfscale,sides)

c subroutine to find the list of natural neighbours attached to each
c nodes. It computes the nn and nb lists but also the surfaces
c attached to the nodes.

c INPUT: x,y       = x,y nodal coordinates
c        nnode     = number of nodes
c        nbmax     = maximum number of neighbour per node
c        points    = working array
c        nodes     = working array
c        vis_tlist = working array
c        vis_elist = working array
c        add_tlist = working array
c        newsurface= here it is just used as a working array
c        eps       = precision
c        xx        = working array
c        pp        = working array
c        aa        = working array
c        bb        = working array
c        surfscale = average nodal surface

c OUTPUT: nn        = neighbour array
c         nb        = number of neighbour for each node
c         vertices  = triangle list
c         neighbour = neighbour list
c         nt        = number of triangle
c         surface   = nodal surface (voronoi cell surface area)
c         sides     = length of the sides of the voronoi cells

c subroutines called:
c NONE

      common /vocal/ ivocal

      real              x(*),y(*),surface(*)
      integer           nn(nbmax,*)
      integer           nb(*)
      real              sides(nbmax,*)

      real*8		points(2,*),eps
      integer		vertices(3,*)
      integer		neighbour(3,*)
      integer           nodes(*)
      integer           vis_tlist(*)
      integer           vis_elist(*)
      integer           add_tlist(*)
      integer		ccw

      real              xx(2),pp(2,nbmax),aa(nbmax,2),bb(nbmax)
      real              newsurface(*)
 
        do i=1,nnode
        points(1,i)=dble(x(i))
        points(2,i)=dble(y(i))
        nodes(i)=i
        enddo

      np=nnode

      nv_max=nnode

c sorting by x

      if (ivocal.eq.1) call debug ('indexx$',0)
      call indexx(np,points,nodes)
      if (ivocal.eq.1) call debug ('find_neighbours$',1)

c ensure initial triangle is in ccw order
 
      if(ccw(points(1,nodes(1)),
     &       points(1,nodes(2)),
     &       points(1,nodes(3)),k).eq.-1)then
             itemp = nodes(1)
             nodes(1) = nodes(2)
             nodes(2) = itemp
      end if

c						Call the routine that 
c						does the work

      mode=0
      eps=0.d0
      if (ivocal.eq.1) call debug ('delaun$',0)
      call delaun (points,np,neighbour,vertices,nt,2*np,
     &             vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &             mode,inactive,bfirst,itstart,subset)
      if (ivocal.eq.1) call debug ('find_neighbours$',1)

      nt_max=nnode*3
        if(nt.gt.nt_max)then
        write(6,*)' Error number of triangles calculated is larger'
        write(6,*)' than array sizes.'
        write(6,*)' Number of triangles    = ',nt
        write(6,*)' Maximum value (nt_max) = ',nt_max
        stop
        endif

c find neighbour list

      call find_neighbour_list (nb,nn,nnode,nbmax,x,y,sides,
     &                          nt,vertices,neighbour)

c find Voronoi cell surface area

        do i=1,nnode
        newsurface(i)=1.
        enddo

      if (ivocal.eq.1) call debug ('find_surface$',0)
      call find_surface (nn,nb,surface,nbmax,nnode,
     &                   x,y,xx,pp,aa,bb,newsurface,surfscale)
      if (ivocal.eq.1) call debug ('find_neighbours$',1)

        do i=1,nnode
        nb(i)=nb(i)+1
          if (nb(i).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i),i)=i
        enddo

      return
      end
