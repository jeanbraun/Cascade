c erosional_properties

      subroutine erosional_properties(x,y,h,h0,hi,nnode,
     &                                param,nparam,nnodemax)

c subroutine to determine erosional properties (or lithology)
c written by Peter
c (comments by Jean)

c the main ingredients are the present topography (h)
c the present position of the bedrock-alluvion interface (h0)
c and the initial topography (hi)
c for every node

c note that erosional properties are nodal values
c and include param(*,1)=xkf (fluvial constant times mean precipitation rate)
c             param(*,2)=xlf_BR (fluvial length scale)
c             param(*,3)=xkdiff (diffusional constant)
c             param(*,4)=width_c (channel width)

c INPUT: x,y      = x-,y- coordinates of nodes
c        h        = present topography
c        h0       = bedrock-alluvion interface
c        hi       = initial topography (at time 0)
c        nnode    = number of nodes
c        param    = erosional parameters
c        nparam   = numper ob erosional parameters
c        nnodemax = maximum number of nodes

c OUTPUT: param is updated

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     x(nnode),y(nnode),h(nnode)
      real     h0(nnode),hi(nnode)
      real     param(nnodemax,nparam)
      character choice*3

      choice='ref'

      select case (choice)

      case ('ref')

        do inode=1,nnode
        param(inode,1)=2.e-1
        param(inode,2)=1.e5
        param(inode,3)=100.
        param(inode,4)=.1
!          if (h(inode).lt.hi(inode)*.8) then
!          param(inode,1)=2.e-2
!          param(inode,3)=.8
!          elseif (h(inode).lt.hi(inode)*.6) then
!          param(inode,1)=4.e-1
!          param(inode,3)=16.
!          endif
        enddo

      case ('can')

        do inode=1,nnode
        param(inode,1)=2.e0
        param(inode,2)=1.e6
        param(inode,3)=10.
        param(inode,4)=.1
        enddo

      case ('bad')

        do inode=1,nnode
        param(inode,1)=2.e-1
        param(inode,2)=3.e4
        param(inode,3)=10.
        param(inode,4)=.1
        enddo

      case ('hil')

        do inode=1,nnode
        param(inode,1)=3.e-1
        param(inode,2)=1.e5
        param(inode,3)=100.
        param(inode,4)=.1
        enddo

      case ('flx')

        do inode=1,nnode
        param(inode,1)=3.e-1
        param(inode,2)=1.e5
        param(inode,3)=30.
        param(inode,4)=.1
        enddo

      case ('hir')

        do inode=1,nnode
        param(inode,1)=3.e-1
        param(inode,2)=1.e5
        param(inode,3)=10.
        param(inode,4)=.1
        enddo

      case default

        do inode=1,nnode
        param(inode,1)=3.e-1
        param(inode,2)=1.e5
        param(inode,3)=3.
        param(inode,4)=.1
        enddo

      end select

      return
      end
