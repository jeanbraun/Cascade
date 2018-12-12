c check_mesh

      subroutine check_mesh (x,y,h,h0,hi,memory,nmemory,itime,
     &                       param,nparam,
     &                       nn,nb,nnode,
     &                       nbmax,nnodemax,ntmax,
     &                       points,vertices,neighbour,nodes,
     &                       vis_tlist,vis_elist,add_tlist,nt,
     &                       slope,water,
     &                       xx,pp,aa,bb,
     &                       nadd,itadd,jtadd,
     &                       dt,surfmin,
     &                       kcon,nkcon,ndon,
     &                       nnode0,
     &                       nodelist,tlist,c_list,v_local,n_local,
     &                       inactive,eps,surfscale,sides)

c This routine checks the mesh to update it during dynamic remeshing
c The flux (dh*surface) of material removed by river erosion is used
c in this version but any measure could be used to improve the resolution
c locally.
c It is important that, if this routine is used, that is if dynamic 
c remeshing is turned on, all nodal parameters that you have added (such
c as a new nodal property) be passed here for what is called permutation.
c During permutation, nodes are renumbered (see near bottom of the
c subroutine) and nodal properties have to be updated accordingly.
c  In this new version all the properties and parameters that have to
c be permuted are stored in memory and param

c param(*,1)=fluvial erosion constant
c param(*,2)=bedrock erosion length scale
c param(*,3)=diffusion erosion constant

c memory(*,1)=dhcrit
c memory(*,2)=dhfluvial
c memory(*,3)=dhdiff
c memory(*,4)=hiso
c memory(*,5)=fix
c memory(*,6)=newsurface
c memory(*,7)=surface

c INPUT: x,y      = x-, y-coordinates of the nodes
c        h        = topography
c        h0       = bedrock-alluvion interface
c        hi       = initial topography
c        memory   = properties to be memorized (and thus properly
c                   carried through) over a time step
c        nmemory  = number of properties to be memorized
c        itime    = current time step
c        param    = erosional parameters
c        nparam   = number of erosional parameters
c        nn       = neighbour matrix
c        nb       = number of neighbours per node
c        nnode    = number of nodes
c        nbmax    = maximum number of neighbours per node
c        nnodemax = maximum number of nodes
c        ntmax    = maximum number of Delaunay triangles
c        points   = real*8 array used in Delaunay routine
c        vertices = triangulation array used in Delaunay routine
c        neighbour= neighbour array used in Delaunay routine
c        nodes    = working array used in Delaunay routine
c        vis_tlist= working array used in Delaunay routine
c        vis_elist= working array used in Delaunay routine
c        add_tlist= working array used in Delaunay routine
c        nt       = number of Delaunay triangles
c        slope    = "nodal" slopes (not used in this version)
c        water    = nodal discharge (not used in this version)
c        xx       = working array used in Voronoi cell volume calculation
c        pp       = working array used in Voronoi cell volume calculation
c        aa       = working array used in Voronoi cell volume calculation
c        bb       = working array used in Voronoi cell volume calculation
c        nadd     = number of nodes added
c        itadd    = working array used in the adding algorithm
c        jtadd    = working array used in the adding algorithm
c        dt       = time step
c        surfmin  = minimum surface allowed (to prevent run away)
c        kcon     = connectivity array (list of elements attached to
c                   a given node)
c        nkcon    = number of elements attached to a given node
c        ndon     = donor array
c        nnode0   = initial number of nodes (these nodes cannot be
c                   removed)
c        nodelist = working array used in the node remove algorithm
c        tlist    = working array used in the node remove algorithm
c        c_list   = working array used in the node remove algorithm
c        v_local  = working array used in the node remove algorithm
c        n_local  = working array used in the node remove algorithm
c        inactive = working array used in Delaunay routine
c        eps      = working array used in Delaunay routine
c        surfscale= mean surface attached to a node
c        sides     = length of the sides of the voronoi cells

c OUTPUT: updated mesh
c         the following arrays and variables are updated:
c             memory, param, nnode, nt, x, y, h, h0, hi, vertices
c             neighbour, nn, nb, kcon, nkcon, nadd

c subroutines called:
c - debug
c - check_for_removal
c - nn_remove
c - find_surface
c - delaun

      common /vocal/ ivocal

      real              x(*),y(*)
      real              h(*),h0(*),hi(*)
      real              water(*),slope(*)
      integer           nn(nbmax,*)
      integer           nb(*)
      real              sides(nbmax,*)

      real              param(nnodemax,nparam)
      real              memory(nnodemax,nmemory)

      real*8		points(2,*)
      real*8            eps
      integer		vertices(3,*)
      integer		neighbour(3,*)
      integer           nodes(*)
      integer           vis_tlist(*)
      integer           vis_elist(*)
      integer           add_tlist(*)
      integer           nodelist(*)
      integer           tlist(*)
      logical           c_list(*)
      integer           v_local(3,*)
      integer           n_local(3,*)
      logical           inactive(*)

      real              xx(2),pp(2,nbmax),aa(nbmax,2),bb(nbmax)
      integer           nkcon(*),ndon(*)
      integer           kcon(ntmax,*),itadd(*),jtadd(*)

      if (ivocal.eq.1) write (22,*) 'time step : ',itime

c fluxmax_erosion is the maximum flux allowed before a new node must
c be inserted

      fluxmax_erosion=surfscale/dt*1.e-2

      nv_max=nnodemax

        do i=1,nnodemax
        memory(i,6)=0.
        inactive(i)=.false.
        if (memory(i,1).gt..5) memory(i,1)=memory(i,1)+1.
        enddo

c dhcritmin is the number of time steps that an "added"
c node (from the original distribution) is forced to remain active before
c being considered for removal

      dhcritmin=0.

c first remove nodes where they are not needed any longer

      iloc=1
      numtri=nt
      nrem=0

c note that the first nnode0 nodes cannot be removed
c they form a sort of minimum set

12    continue
        if (ivocal.eq.1) call debug ('check_for_removal$',0)
        do i=nnode0+1,nnode
        call check_for_removal (memory(i,7),memory(i,2),
     &                          nb(i),nn(1,i),fluxmax_erosion,
     &                          surfmin,dt,iremove)
          if(iremove.eq.1) then
          nrem=nrem+1
          np=nnode
          if (ivocal.eq.1) write(22,*) -i
            do k=1,nb(i)
            memory(nn(k,i),6)=1.
            enddo
          call nn_remove (i,np,numtri,points,
     &                    vertices,
     &                    neighbour,iloc,nbmax,nv_max,
     &                    vis_tlist,vis_elist,add_tlist,
     &                    v_local,n_local,c_list,nodelist,
     &                    tlist,.false.)
            do j=i+1,nnode
            points(1,j-1)=points(1,j)
            points(2,j-1)=points(2,j)
            x(j-1)=x(j)
            y(j-1)=y(j)
            h(j-1)=h(j)
            h0(j-1)=h0(j)
            hi(j-1)=hi(j)
              do kpar=1,nparam
              param(j-1,kpar)=param(j,kpar)
              enddo
              do kmem=1,nmemory
              memory(j-1,nmemory)=memory(j,nmemory)
              enddo
            enddo
            do it=1,numtri
              do k=1,3
              if (vertices(k,it).gt.i)
     &            vertices(k,it)=vertices(k,it)-1
              enddo
            enddo
          nnode=nnode-1
          call find_neighbour_list (nb,nn,nnode,nbmax,x,y,sides,
     &                              numtri,vertices,neighbour)
          goto 12
          endif
        enddo

      if (ivocal.eq.1) call debug ('check_mesh$',1)

      if (nrem.gt.0) then

c recompute surfaces around removed nodes

      if (ivocal.eq.1) call debug ('find_surface$',0)
      call find_surface (nn,nb,memory(1,7),nbmax,nnode,
     &                   x,y,xx,pp,aa,bb,memory(1,6),surfscale)
      if (ivocal.eq.1) call debug ('check_mesh$',1)

      do i=1,nnode
      if (memory(i,7).eq.0.) then
      print*,'surface nil at ',i,memory(i,6)
      print*,nb(i),(nn(k,i),k=1,nb(i))
        do k=1,nb(i)
        print*,x(nn(k,i)),y(nn(k,i))
        enddo
      endif
      enddo

c recompute connectivities around removed nodes

        do i=1,nnode
        nkcon(i)=0
        enddo
        do it=1,numtri
        if (vertices(1,it).le.nnode
     &.and. vertices(2,it).le.nnode
     &.and. vertices(3,it).le.nnode) then
          do i=1,3
          ic=vertices(i,it)
          nkcon(ic)=nkcon(ic)+1
            if (nkcon(ic).gt.ntmax) then
            print*,ic,nkcon(ic)
            print*,(kcon(ii,ic),ii=1,nkcon(ic)-1)
            stop 'too tight connectivity...'
            endif
          kcon(nkcon(ic),ic)=it
          enddo
        endif
        enddo

      nt=numtri

      endif

c compute in which triangles points are to be added
c and location of new points

      nadd=0
        do it=1,nt
        jtadd(it)=0
        enddo

        do i=1,nnode
        surf=memory(i,7)
        if (surf.gt.surfmin) then
        flux=min(0.,surf*memory(i,2)/dt)
          if(flux.lt.-fluxmax_erosion) then
          do j=1,nkcon(i)
          it=kcon(j,i)
          i1=vertices(1,it)
          i2=vertices(2,it)
          i3=vertices(3,it)
          if (i1.le.nnode.and.i2.le.nnode.and.i3.le.nnode .and.
     &        jtadd(it).eq.0) then
          jtadd(it)=1
          nadd=nadd+1
            if (nnode+nadd.gt.nnodemax) then
            print*,nnode,nadd,nnode+nadd,nnodemax
            print*,'Too many nodes in routine check_mesh/1'
            stop
	    endif
          itadd(nadd)=it
	  newn=nnode+nadd
          if (ivocal.eq.1) write (22,*) newn
          x(newn)=(x(i1)+x(i2)+x(i3))/3.
          y(newn)=(y(i1)+y(i2)+y(i3))/3.
          h(newn)=(h(i1)+h(i2)+h(i3))/3.
          h0(newn)=(h0(i1)+h0(i2)+h0(i3))/3.
          hi(newn)=(hi(i1)+hi(i2)+hi(i3))/3.
          memory(newn,1)=1.
	  memory(newn,2)=(memory(i1,2)+memory(i2,2)+memory(i3,2))/3.
	  memory(newn,3)=(memory(i1,3)+memory(i2,3)+memory(i3,3))/3.
	  memory(newn,4)=(memory(i1,4)+memory(i2,4)+memory(i3,4))/3.
          memory(newn,5)=1.
          if (memory(i1,5).eq.0. .and. memory(i2,5).eq.0.
     &        .and. memory(i3,5).eq.0.) memory(newn,5)=0.
            do kpar=1,nparam
            param(newn,kpar)=
     &      (param(i1,kpar)+param(i2,kpar)+param(i3,kpar))/3.
            enddo
          points(1,newn)=dble(x(newn))
          points(2,newn)=dble(y(newn))
          memory(newn,1)=1.
          endif
          enddo
          endif
          endif
        enddo

        if (nadd.eq.0) then
        return
        endif

c compute new Delaunay triangulation

      mode=3
      np=nnode+nadd
      itstart=1
      numtri=nt
      if (ivocal.eq.1) call debug ('delaun$',0)
      call delaun (points,np,neighbour,vertices,numtri,2*np,
     &             vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &             mode,inactive,nnode+1,itstart,subset)
      if (ivocal.eq.1) call debug ('check_mesh$',1)

      nt=numtri

      nnode=nnode+nadd

c compute neighbour list

      call find_neighbour_list (nb,nn,nnode,nbmax,x,y,sides,
     &                          nt,vertices,neighbour)

        do i=nnode-nadd+1,nnode
        memory(i,6)=1.
          do j=1,nb(i)
          memory(nn(j,i),6)=1.
          enddo
        enddo

c finds new surfaces where needed

      if (ivocal.eq.1) call debug ('find_surfacE$',0)
      call find_surface (nn,nb,memory(1,7),nbmax,nnode,
     &                   x,y,xx,pp,aa,bb,memory(1,6),surfscale)
      if (ivocal.eq.1) call debug ('check_mesh$',1)

        do i=1,nnode
        nb(i)=nb(i)+1
          if (nb(i).gt.nbmax) then
          stop 'nbmax too small...0'
          endif
        nn(nb(i),i)=i
        enddo

      return
      end
