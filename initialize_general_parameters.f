c initialize_general_parameter

      subroutine initialize_general_parameters
     & (surfscale,
     &  dt,iadjust,endtime,ishow,
     &  writetime,nshortwrite,iflux,run_name,nrun_name,
     &  iadapt,surfmin,
     &  ihorizontal,
     &  iflexure,hflex,ixflex,iyflex,thickflex,
     &  ym,pratio,rhocflex,rhoaflex,
     &  oro_length,oro_height,oro_scale,wind_direction,
     &  xlf_AL,sea_level,
     &  ideposition,idiffusion)

c subroutine to initialize general parameters

c INPUT: surfscale       = scaling surface (in km**2)

c OUTPUT: dt             = (initial) time step length (in yr)
c         iadjust        = allows for dynamic time stepping (=0
c                          means that time step is fixed; =1 means
c                          that time step is adjusted dynamically)
c         endtime        = final time (in yr)
c         ishow          = frequency of graphic displays (in time steps)
c         writetime      = frequency of output (in yr)
c         nshortwrite    = frequency of short screen output (in time steps)
c         iflux          = frequency of flux saves (=0 no flux save;
c                          otherwise frequency in time steps); this output
c                          is of total flux through the fixed nodes
c         run_name       = name of this run; must also be the name
c                          of an existing folder where all output files
c                          will be stored
c         iadapt         = flag to allow dynamic remeshing (=0 no; =1 yes)
c         surfmin        = minimum surface area that can be reached by 
c                          remeshing/refining
c         ihorizontal    = flag to permit horizontal mesh movement (=0
c                          means no horizontal movement; =1 means horizontal
c                          movements permitted
c         iflexure       = flag to permit flexural isostasy (=0 no flexure;
c                          =1 flexure)
c         hflex          = size (in m) of the square mesh on which the
c                          thin elastic plate calculations are done
c         ixflex         = flag to permit flexure in the x-direction
c                          (=0 no elastic strength in x-direction;
c                          =1 means elastic strength in x-direction)
c         iyflex         = flag to permit flexure in the y-direction
c                          (=0 no elastic strength in y-direction;
c                          =1 means elastic strength in y-direction)
c         thickflex      = elastic thickness in m
c         ym             = Young Modulus (in Pa)
c         pratio         = Poissons's ratio
c         rhocflex       = crustal density (in kg/m**3)
c         rhoaflex       =asthenospheric density (in kg/m**3)
c         oro_length     = orographic length scale (in m)
c                          (=0 means no orographic control on precipitation)
c         oro_height     = orographic height scale (in m)
c         oro_scale      = background precipitation (in adequate units)
c         wind_direction = wind direction (0= along x-axis)
c         xlf_AL         = fluvial erosion length scale for alluvials
c                          (in m)
c         sea_level      = where sea-level is (in m above the 0 datum)
c                          below sea-level no sediment transport is
c                          allowed
c         ideposition    = flag to allow for fluvial deposition (=0 erosion
c                          only; =1 sedimentation allowed)
c         idiffusion     = flag to allow diffusion processes (=0 no
c                          diffusion; =1 diffusion allowed)

c subroutines called:
c - read_but_skip_comment
c -iread_but_skip_comment

      common /vocal/ ivocal

      character*256    run_name

      iread=0

      if (iread.eq.0) then

      dt=1.e2
      iadjust=0
      endtime=1.e7
!      endtime=4.e5
      ishow=1000
      writetime=1.e6
      nshortwrite=100
      iflux=1
        do k=1,256
        run_name(k:k)=' '
        enddo
      run_name(1:4)="RUN1"
        do k=256,1,-1
        if (run_name(k:k).eq.' ') nrun_name=k-1
        enddo
        if (nrun_name.eq.0) then
        print*,'No run name available '
        stop
        endif
      iadapt=0
      surfmin=surfscale/4.
      ihorizontal=0
      iflexure=0
! time change to 3e5
! and bc to one side
      hflex=1.e7
      ixflex=1
      iyflex=1
      thickflex=25.e3
      ym=1.e11
      pratio=0.25
      rhocflex=2750.
      rhoaflex=3300.
      oro_length=0.
      oro_height=1.
      oro_scale=5000000.
      wind_direction=90.
      xlf_AL=10.e3
      sea_level=0.
      ideposition=1
      idiffusion=1

      else

      print*,'Enter run name >'
      read (*,'(a)') run_name

        do k=256,1,-1
        if (run_name(k:k).eq.' ') nrun_name=k-1
        enddo
        if (nrun_name.eq.0) then
        print*,'No run name available '
        stop
        endif

      open (7,file=run_name(1:nrun_name)//
     &       '/cascade.parameters.in',status='old',err=997)
      goto 998
997   print*,'Cant open input file'
      print*,'You must create a folder/directory named: '//
     &       run_name(1:nrun_name)
      print*,'where the input file cascade.parameters.in resides'
      stop
998   continue

      call read_but_skip_comment (7,1,dt)
      call iread_but_skip_comment (7,1,iadjust)
      call read_but_skip_comment (7,1,endtime)
      call iread_but_skip_comment (7,1,ishow)
      call read_but_skip_comment (7,1,writetime)
      call iread_but_skip_comment (7,1,nshortwrite)
      call iread_but_skip_comment (7,1,iflux)
      call iread_but_skip_comment (7,1,iadapt)
      call read_but_skip_comment (7,1,surfmin)
      call iread_but_skip_comment (7,1,ihorizontal)
      call iread_but_skip_comment (7,1,iflexure)
      call read_but_skip_comment (7,1,hflex)
      call iread_but_skip_comment (7,1,ixflex)
      call iread_but_skip_comment (7,1,iyflex)
      call read_but_skip_comment (7,1,thickflex)
      call read_but_skip_comment (7,1,ym)
      call read_but_skip_comment (7,1,pratio)
      call read_but_skip_comment (7,1,rhocflex)
      call read_but_skip_comment (7,1,rhoaflex)
      call read_but_skip_comment (7,1,oro_length)
      call read_but_skip_comment (7,1,oro_height)
      call read_but_skip_comment (7,1,oro_scale)
      call read_but_skip_comment (7,1,wind_direction)
      call read_but_skip_comment (7,1,xlf_AL)
      call read_but_skip_comment (7,1,sea_level)
      call iread_but_skip_comment (7,1,ideposition)
      call iread_but_skip_comment (7,1,idiffusion)

      close (7)

      endif

c if iecho = 1 the input parameters are echoed on the screen

      iecho=1

      if (iecho.eq.1) then
      print*,'dt= ',dt
      print*,'iadjust= ',iadjust
      print*,'endtime= ',endtime
      print*,'ishow= ',ishow
      print*,'writetime= ',writetime
      print*,'nshortwrite= ',nshortwrite
      print*,'iflux= ',iflux
      print*,'run_name(1:4)= ',run_name(1:4)
      print*,'iadapt= ',iadapt
      print*,'surfmin= ',surfmin
      print*,'ihorizontal= ',ihorizontal
      print*,'iflexure= ',iflexure
      print*,'hflex= ',hflex
      print*,'ixflex= ',ixflex
      print*,'iyflex= ',iyflex
      print*,'thickflex= ',thickflex
      print*,'ym= ',ym
      print*,'pratio= ',pratio
      print*,'rhocflex= ',rhocflex
      print*,'rhoaflex= ',rhoaflex
      print*,'oro_length= ',oro_length
      print*,'oro_height= ',oro_height
      print*,'oro_scale= ',oro_scale
      print*,'wind_direction= ',wind_direction
      print*,'xlf_AL= ',xlf_AL
      print*,'sea_level= ',sea_level
      print*,'ideposition= ',ideposition
      print*,'idiffusion= ',idiffusion
      endif

      return
      end
