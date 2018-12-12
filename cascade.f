      program cascade

c  -----------------------------------------------------
c  |                                                   |
c  | CCCCC     A     SSSSS CCCCC     A     DDDD  EEEEE |
c  | C        A A    S     C        A A    D   D E     |
c  | C       AAAAA   SSSSS C       AAAAA   D   D EEEE  |
c  | C      A     A      S C      A     A  D   D E     |
c  | CCCCC A       A SSSSS CCCCC A       A DDDD  EEEEE |
c  |                                                   |
c  -----------------------------------------------------

c The program cascade was developed by:

c        Jean Braun
c        Research School of Earth Sciences
c        Australian National University
c        Canberra, ACT, 0200
c        Australia
c        Tel: +61-2-6249-5512
c        Fax: +61-2-6249-5443
c        email: Jean.Braun@anu.edu.au

c (Canberra, June 1st, 1995)
c (Present version September 27, 1999)
c (Present version January 10, 2000)

c Please report any problem and certainly any improvement you bring to
c cascade.

c cascade is a geomorphic program to compute the evolution of landscapes
c by erosion/deposition and tectonic uplift/subsidence. Two types of processes
c are included: 1) short-range (hillslope) processes, modeled by a simple linear
c                  diffusion equation
c               2) long-range (river) processes, modeled by an equation of
c                  reaction between water flowing in a river network and the
c                  substratum

c The main advantage/difference between cascade and other geomorphic models
c is its ability to handle arbitrary, non-rectangular, grids at the corners
c of which the computations are done. This opens opportunities to vary spatial
c discretization in various parts of the model or to adjust spatial discretiza-
c tion to the results of the model (evolving discretization).

c cascade uses the theory of Delaunay triangulation and Voronoi diagrams
c to compute the possible "connections" between neighbouring nodes on the
c arbitrary grid and the surface area of the part of the landscape the height
c of which is slave to any node of the grid.

c The river network is calculated from the "bucket passing" algorithm (BPA)
c that does not require a complete ordering of the nodes according to their
c heights but rather a local-only ordering that depends on the position of the
c nodes with respect to local maxima in heights. In BPA, each node is given a
c bucket full of water and asked to pass it to its lowest neighbour. After that
c operation, all nodes that have not received anything are local maxima and
c are put at the top of a stack. After the next bucket passing step, those that
c have not received anything are put on the stack, and so on until all nodes
c have been put on the stack. The stack contains an ordering of the nodes that
c is appropriate to calculate the effect of river flow on the landscape. 

c note that because this version of the bucket algorithm can handle local
c minima, runs tend to be slower at the beginning when many local minima
c exist. As the minima disappear due to river erosion and deposition, one
c should notice an increase efficiency per time step

c In this version, the diffusion equation is solved by using a 
c finite volume algorithm with explicit time integration.

c In the 2000 version, fluvial erosion scales properly with grid
c size and erosion rate is proportional to the square root of discharge
c area

c cascade has the possibility of including a stratified crust with
c regards to erosional parameters
c Developed by Peter van der Beek, April '96

c IMPORTANT NOTE:
c****************
c Please, note that this software CANNOT be freely distributed. You must
c obtain Jean Braun's permission to use it (or part of it) or to give to other
c potential users. Please, respect this condition of use. I am trying to 
c protect parts of the Delaunay/Voronoi algorithms that we are using in a
c commercial venture with Malcolm Sambridge. This means that some of our
c "clients" had to pay to use this software commercially.

c Good luck.

c--------------------------------------------------------------------------

c subroutines called:
c - debug
c - initialize_general_parameters
c - initialize_nodal_geometry
c - erosional_properties
c - update_time_step
c - find_neighbours
c - write_output
c - change_sea_level
c - erosional_properties
c - tectonic_uplift
c - tectonic_movement
c - orography
c - find_donors
c - find_catchment
c - fluvial_erosion
c - diffusion_erosion
c - flexure
c - show
c - check_mesh
c - update_bedrock
c - write_output
c - update_time_step

      include 'cascade.h'

      common /vocal/ ivocal

c initializes memory

        do k=1,nmemory
          do i=1,nnodemax
          memory(i,k)=0.
          enddo
        enddo

c ivocal=1 means CASCADE will save debugginh information in debug.out
c should be set to ivocal=0 in production runs

      ivocal=0

      if (ivocal.eq.0) then
      open (89,file='debug.out',
     &      status='unknown')
      write (89,'(a)') 'Debugger not activated'
      write (89,'(a)') 'To activate set ivocal to 1 in cascade'
      close (89)
      else
      open (89,file='debug.out',
     &      status='unknown')
      rewind (89)
      endif

c if ivocal=1 opens a log file to record which nodes are being
c added/removed from the mesh

      if (ivocal.eq.1) open (22,file='check_mesh.out',status=
     &                       'unknown')

c initialize general parameters

      if (ivocal.eq.1) call debug ('initialize_general_parameters$',0)
      call initialize_general_parameters
     & (surfscale,
     &  dt0,iadjust,endtime,ishow,
     &  writetime,nshortwrite,iflux,run_name,nrun_name,
     &  iadapt,surfmin,
     &  ihorizontal,
     &  iflexure,hflex,ixflex,iyflex,thickflex,
     &  ym,pratio,rhocflex,rhoaflex,
     &  oro_length,oro_height,oro_scale,wind_direction,
     &  xlf_AL,sea_level,
     &  ideposition,
     &  idiffusion)
      if (ivocal.eq.1) call debug ('cascade$',1)

c initialize the node geometry

      if (ivocal.eq.1) call debug ('initialize_nodal_geometry$',0)
      call initialize_nodal_geometry
     & (nnodemax,nnode,x,y,h,memory(1,5),delta,surfscale,
     &  run_name,nrun_name)
      if (ivocal.eq.1) call debug ('cascade$',1)
      nnode0=nnode
        do i=1,nnode
        h0(i)=h(i)
        hi(i)=h(i)
        enddo

c opens various output files

        if (nrun_name.eq.0) then
        print*,'No run name available '
        stop
        endif

      open (7,file=run_name(1:nrun_name)//'/topography',
     &      status='unknown',err=997)
      goto 998
997   print*,'You must create a folder/directory named '//
     &       run_name(1:nrun_name)
      stop
998   continue

      open (8,file=run_name(1:nrun_name)//'/connectivity',
     &      status='unknown')
      open (9,file=run_name(1:nrun_name)//'/donors',
     &      status='unknown')
      open (10,file=run_name(1:nrun_name)//'/geometry',
     &      status='unknown')
      open (11,file=run_name(1:nrun_name)//'/properties',
     &      status='unknown')
      open (12,file=run_name(1:nrun_name)//'/discharge',
     &      status='unknown')
      open (13,file=run_name(1:nrun_name)//'/erosion_rate',
     &      status='unknown')
      open (14,file=run_name(1:nrun_name)//'/catchments',
     &      status='unknown')
      open (15,file=run_name(1:nrun_name)//'/lakes',
     &      status='unknown')


      if (iflux.eq.1) 
     & open (20,file=run_name(1:nrun_name)//'/outflux',
     &       status='unknown')

c initialize erosional nodal properties

      if (ivocal.eq.1) call debug ('erosional_properties$',0)
      call erosional_properties (x,y,h,h0,hi,nnode,param,
     &                           nparam,nnodemax)
      if (ivocal.eq.1) call debug ('cascade$',1)

c initialize time step

      if (ivocal.eq.1) call debug ('update_time_step$',0)
      call update_time_step (dt,dt0,dtold,param(1,1),nnode,iadjust,
     &                       delta,0,
     &                       dhminfluv,dhmaxfluv,
     &                       dhmindiff,dhmaxdiff,
     &                       time,shorttime,writetime,endtime,
     &                       itime,dtmin,dtmax)
      if (ivocal.eq.1) call debug ('cascade$',1)
 
c finds Delaunay triangulation and voronoi cell surface areas for initial
c set of nodes
        
      if (ivocal.eq.1) call debug ('find_neighbours$',0)
      call find_neighbours (x,y,nn,nb,nnode,nbmax,
     &                      points,vertices,neighbour,nodes,
     &                      vis_tlist,vis_elist,add_tlist,nt,
     &                      memory(1,7),memory(1,6),eps,
     &                      xy,pp,aa,bb,surfscale,sides)
      nbbmax=1
        do i=1,nnode
        nbbmax=max(nbbmax,nb(i))
        enddo
      if (ivocal.eq.1) call debug ('cascade$',1)

c writes the initial conditions

      if (ivocal.eq.1) call debug ('write_output$',0)
      call write_output
     &     (h,x,y,memory(1,5),nnode,nnodemax,iadapt,1,vertices,nt,
     &      ndon,param,nparam,water,memory(1,2),memory(1,3),
     &      ncat,nlake,0.,dt)
      if (ivocal.eq.1) call debug ('cascade$',1)

c start of time stepping

      istep=0

      do while (time.lt.endtime)

      istep=istep+1

c change in sea level

      if (ivocal.eq.1) call debug ('change_sea_level$',0)
      call change_sea_level (time,sea_level)
      if (ivocal.eq.1) call debug ('cascade$',1)

c initializes the flux of material entering and leaving the grid
c for this time step

      influx=0.
      outflux=0.

c distributes erosional properties to rocks according to the amount eroded
c i.e. makes layered crust (Peter, April '96)

      if (ivocal.eq.1) call debug ('erosional_properties$',0)
      call erosional_properties (x,y,h,h0,hi,nnode,param,
     &                           nparam,nnodemax)
      if (ivocal.eq.1) call debug ('cascade$',1)

c add tectonic uplift/subsidence component to landscape height
c (note that heights are in km)

      if (ivocal.eq.1) call debug ('tectonic_uplift$',0)
      call tectonic_uplift (x,y,h,h0,hi,nnode,
     &                      memory(1,5),dt,time,
     &                      influx,memory(1,7))
      if (ivocal.eq.1) call debug ('cascade$',1)

c add tectonic horizontal movement

      if (ihorizontal.eq.1) then
      if (ivocal.eq.1) call debug ('tectonic_movement$',0)
      call tectonic_movement (x,y,nnode,
     &                        memory(1,5),dt,itime,
     &                        points,vertices,neighbour,
     &                        memory(1,7),memory(1,6),nt,
     &                        nn,nb,nbmax,nswaps,
     &                        mask,mask_e,xy,pp,aa,bb,surfscale,
     &                        sides)
      if (ivocal.eq.1) call debug ('cascade$',1)
      endif

c find rain function (water)

      if (ivocal.eq.1) call debug ('orography$',0)
      call orography (x,y,h,water,memory(1,7),length,work,
     &                ndon,nn,nb,nbmax,
     &                nwork,
     &                oro_length,oro_height,oro_scale,
     &                nnode,wind_direction)
      if (ivocal.eq.1) call debug ('cascade$',1)

c find the donor for each node (its lowest neighbour)

      if (ivocal.eq.1) call debug ('find_donors$',0)
      call find_donors (x,y,h,slope,length,nnode,
     &                  ndon,nb,nn,nbmax,delta)
      if (ivocal.eq.1) call debug ('cascade$',1)

c find catchment for each node

      if (ivocal.eq.1) call debug ('find_catchment$',0)
      call find_catchment (ndon,nwork,ncat,nsill,nempty,nlake,
     &                     nn,nb,nbmax,h,nnode,nself,nselfu,memory(1,5),
     &                     x,y,slope,length)
      if (ivocal.eq.1) call debug ('cascade$',1)

c calculate river erosion

      if (ivocal.eq.1) call debug ('fluvial_erosion$',0)
      call fluvial_erosion (param(1,1),param(1,2),xlf_AL,param(1,4),
     &                      x,y,h,h0,hi,ndon,nnode,
     &                      memory(1,7),slope,length,
     &                      water,sediment,
     &                      ibucket,dt,memory(1,5),
     &                      memory(1,2),
     &                      nb,nn,nbmax,
     &                      iorder,itype_node,
     &                      dhminfluv,dhmaxfluv,
     &                      nlake,
     &                      sea_level,outflux,ideposition)
      if (ivocal.eq.1) call debug ('cascade$',1)

c calculate diffusion erosion

      if (idiffusion.eq.1) then
      if (ivocal.eq.1) call debug ('diffusion_erosion$',0)
      call diffusion_erosion (param(1,3),
     &                        x,y,h,memory(1,3),nnode,
     &                        nb,nn,nbmax,
     &                        memory(1,5),dt,dhmindiff,dhmaxdiff,taumin,
     &                        sea_level,outflux,memory(1,7),sides)
      endif
      if (ivocal.eq.1) call debug ('cascade$',1)

c calculate isostatic rebound

      if (iflexure.eq.1) then
      if (ivocal.eq.1) call debug ('flexure$',0)
      call flexure
     &             (x,y,h,h0,hi,memory(1,4),
     &              nnode,flex,work_flex,nflex,memory(1,7),
     &              points,vertices,neighbour,eps,
     &              hflex,rhocflex,thickflex,ym,pratio,rhoaflex,
     &              hisomin,hisomax,ixflex,iyflex,memory(1,5))
      if (ivocal.eq.1) call debug ('cascade$',1)
      endif

c draw a simple plot

      if (ishow.ne.0 .and. (istep/ishow)*ishow.eq.istep) then
      if (ivocal.eq.1) call debug ('show$',0)
      call show (x,y,h,memory(1,7),ndon,nb,nn,nbmax,
     &           ncat,nlake,nnode,
     &           vertices,nt,
     &           -1,0.,1500.,1,1,istep/ishow,istep/ishow)
      if (ivocal.eq.1) call debug ('cascade$',1)
      endif

c check the grid (in case adaptive grid is allowed)

      if (iadapt.eq.1) then
      if (ivocal.eq.1) call debug ('check_mesh$',0)
      call check_mesh (x,y,h,h0,hi,memory,nmemory,itime,
     &                 param,nparam,
     &                 nn,nb,nnode,nbmax,nnodemax,ntmax,
     &                 points,vertices,neighbour,nodes,
     &                 vis_tlist,vis_elist,add_tlist,nt,
     &                 slope,water,
     &                 xy,pp,aa,bb,
     &                 nadded,itadd,jtadd,
     &                 dt,surfmin,
     &                 kcon,nkcon,ndon,
     &                 nnode0,
     &                 nodelist,tlist,c_list,v_local,n_local,
     &                 inactive,eps,surfscale,sides)
      if (ivocal.eq.1) call debug ('cascade$',1)
      endif

c check for bedrock incision

      if (ivocal.eq.1) call debug ('update_bedrock$',0)
      call update_bedrock (h,h0,nnode)
      if (ivocal.eq.1) call debug ('cascade$',1)

c update time

      time=time+dt
      shorttime=shorttime+dt

c writing output

      if (shorttime.eq.writetime) then
      if (ivocal.eq.1) call debug ('write_output$',0)
!      call write_output
!     &     (h,x,y,memory(1,5),nnode,nnodemax,iadapt,0,vertices,nt,
!     &      ndon,param,nparam,water,memory(1,2),memory(1,3),
!     &      ncat,nlake,time,dt)
      if (ivocal.eq.1) call debug ('cascade$',1)
      shorttime=0.
        if (dtold.ne.0.) then
        dt=dtold
        dtold=0.
        endif
      endif

c write a short line to the screen

      if ((itime/nshortwrite*nshortwrite).eq.itime) then
      hmin=h(1)
      hmax=h(1)
        do i=1,nnode
        hmin=min(hmin,h(i))
        hmax=max(hmax,h(i))
        enddo
      hmedian=(hmax-hmin)/2.
      write (6,*) 'TIME ',istep,' (',time,') hrange=[',
     &             hmin,',',hmedian,',',hmax,']'
      endif

c write outflux to file (Peter, April 96)

      if (iflux.eq.1.and.(itime/nshortwrite*nshortwrite).eq.itime) then
      write(20,*) time,outflux/dt,influx/dt
      endif

c update time step

      if (ivocal.eq.1) call debug ('update_time_step$',0)
      call update_time_step (dt,dt0,dtold,param(1,1),nnode,iadjust,
     &                       delta,1,
     &                       dhminfluv,dhmaxfluv,
     &                       dhmindiff,dhmaxdiff,
     &                       time,shorttime,writetime,endtime,
     &                       itime,dtmin,dtmax)
      if (ivocal.eq.1) call debug ('cascade$',1)

c end of time stepping

      enddo

      if (ishow.ne.0 .and. istep.ge.ishow) then
      if (ivocal.eq.1) call debug ('show$',0)
      call show (x,y,h,memory(1,7),ndon,nb,nn,nbmax,
     &           ncat,nlake,nnode,
     &           vertices,nt,
     &           -1,0.,1500.,1,1,0,istep/ishow)
      if (ivocal.eq.1) call debug ('cascade$',1)
      endif

      end
