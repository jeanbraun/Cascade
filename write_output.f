c write_output

      subroutine write_output
     &     (h,x,y,fix,nnode,nnodemax,iadapt,ifirst,vertices,nt,
     &      ndon,param,nparam,water,dhfluvial,dhdiffusion,
     &      ncat,nlake,time,dt)

c subroutine to write output to a series of files in ASCII

c INPUT: h         = current topography
c        x,y       = x- and y-nodal coordinates
c        fix       = boundary conditions
c        nnode     = number of nodes
c        nnodemax  = maximum number of nodes
c        iadapt    = flag to allow dynamic remeshing (=0 no; =1 yes)
c        ifirst    = =1 means first zeroth time step (used to store
c                    initial geometry and parameters)
c        vertices  = triangles connectivity
c        nt        = number of triangles
c        ndon      = donor array
c        param     = erosion parameters
c        nparam    = number of erosion parameters
c        water     = discharge
c        dhfluvial = fluvial erosion increment over the time step
c        dhdiffusion= diffusion erosion increment over the time step
c        ncat      = catchment names
c        nlake     = lake flag
c        time      = current time

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     h(nnode),x(nnode),y(nnode)
      real     dhfluvial(nnode),dhdiffusion(nnode)
      real     param(nnodemax,nparam),water(nnode)
      integer  vertices(3,nt),ndon(nnode),ncat(nnode),nlake(nnode)
      real     fix(nnode)

c topography
      write (7,*) 'TIME ',time,nnode
        do i=1,nnode
        write (7,*) i,h(i)
        enddo

c connections
        if (iadapt.eq.1 .or. ifirst.eq.1) then
        write (8,*) 'TIME ',time,nt
          do k=1,nt
          write (8,*) k,(vertices(i,k),i=1,3)
          enddo
        endif

c donors
        if (ifirst.eq.0) then
        write (9,*) 'TIME ',time,nnode
          do i=1,nnode
          write (9,*) i,ndon(i)
          enddo
        endif

c geometry
        if (iadapt.eq.1 .or. ifirst.eq.1) then
        write (10,*) 'TIME ',time,nnode
          do i=1,nnode
          write (10,*) i,x(i),y(i),fix(i)
          enddo
        endif

c parameters
      write (11,*) 'TIME ',time,nparam,nnode
        do i=1,nnode
        write (11,*) (param(i,k),k=1,nparam)
        enddo

c discharge
        if (ifirst.eq.0) then
        write (12,*) 'TIME ',time,nnode
          do i=1,nnode
          write (12,*) i,water(i)
          enddo
        endif

c erosion rate
        if (ifirst.eq.0) then
        write (13,*) 'TIME ',time,nnode
          do i=1,nnode
          write (13,*) i,dhfluvial(i)/dt,dhdiffusion(i)/dt,
     &                  (dhfluvial(i)+dhdiffusion(i))/dt
          enddo
        endif

c catchements
        if (ifirst.eq.0) then
        mcat=0
        write (14,*) 'TIME ',time,nnode
1111    lcat=0
          do i=1,nnode
            if (ncat(i).gt.0) then
            lcat=1
            ncat0=ncat(i)
            mcat=mcat+1
            write (14,*) 'Cat ',mcat
              do j=i,nnode
                if (ncat(j).eq.ncat0) then
                write (14,*) j
                ncat(j)=-ncat(j)
                endif
              enddo
            endif
          if (lcat.eq.1) goto 1111
          enddo
          do i=1,nnode
          ncat(i)=-ncat(i)
          enddo
        endif

c lakes
        if (ifirst.eq.0) then
        write (15,*) 'TIME ',time,nnode
          do i=1,nnode
          if (nlake(i).eq.1) write (15,*) i
          enddo
        endif

      return
      end
