c orography

      subroutine orography (x,y,h,water,surface,length,tot,
     &                      ndon,nn,nb,nbmax,
     &                      nwork,
     &                      oro_length,oro_height,oro_scale,
     &                      nnode,alpha)

c subroutine orography to calculate the amount of rain that falls on each cell
c taking into account a dominant wind direction (alpha), an orographic
c model (oro) and a finite amount of water to be precipitated (precipitation)

C***  This appears to be a simplified version of the orographic model...  It
C  only allows one kind of rain distribution...  It is probably sufficient for
C  now (though I need to figure out what orog_lenght, oro_height and oro_scale
C  stand for), but in the long run I will need to either 1) get the "original"
C  (full-version) code from Jean, 2) piece in the orographic model from 
C  Drainal, or 3) write my own code.  The easiest is probably to contact Jean.

c INPUT: x, y         = x- and y-nodal coordinates
c        h            = present topography
c        surface      = surface associated with each node
c        length       = length associated with each nodal link
c        tot          = working array
c        ndon         = donor array
c        nn           = natural neighbour list
c        nb           = number of natural neighbours for each node
c        nbmax        = maximum number of natural neighbours
c        nwork        = working array
c        oro_length   = orographic length scale
c        oro_height   = orographic height scale
c        oro_scale    = orgraphic precipitation scale
c        nnode        = total number of nodes
c        alpha        = direction of dominant winds

c OUPUT: water        = precipitation at each node

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     x(nnode),y(nnode),h(nnode),water(nnode)
      real     surface(nnode),length(nnode),tot(nnode)
      integer  ndon(nnode),nb(nnode),nn(nbmax,nnode)
      integer  nwork(nnode)

        if (oro_length.le.0.) then
          do i=1,nnode
          water(i)=surface(i)
          enddo
        return
        endif

      cosa=cos(alpha*3.141592654/180.)
      sina=sin(alpha*3.141592654/180.)

c calculate, for each node, a "donor" node, that is the one of the neighbours
c that is closest to be downwind. Note that a donor node cannot be upwind.

        do i=1,nnode
        x0=x(i)
        y0=y(i)
c        crossmax=sqrt(2.)/2.
        crossmax=-2.
        ndon(i)=i
          do j=1,nb(i)
          k=nn(j,i)
          if (k.ne.i) then
          cross=(x0-x(k))*cosa+(y0-y(k))*sina
            if (cross.gt.crossmax) then
            crossmax=cross
            ndon(i)=k
            endif
          endif
          enddo
        if (ndon(i).gt.nnode) ndon(i)=i
        enddo

c finds those nodes that do not receive anything

        do i=1,nnode
        nwork(i)=0
        tot(i)=-1.
        water(i)=0.
        length(i)=sqrt((x(i)-x(ndon(i)))**2+(y(i)-y(ndon(i)))**2)
        enddo

        do i=1,nnode
        nwork(ndon(i))=1
        enddo

      print*,minval(length),maxval(length)

        do i=1,nnode
          if (nwork(i).eq.0) then
          k=i
          kp=k
          np=0
1         np=np+1
          nwork(np)=k
          water(k)=surface(k)/length(kp)*h(k)/oro_height/oro_length
          kp=k
          k=ndon(kp)
          if (k.ne.kp .and. tot(k).le.0.) goto 1
          if (tot(k).le.0.) tot(k)=oro_scale
          kp=k
            do j=np,1,-1
            k=nwork(j)
            water(k)=tot(kp)*water(k)
            tot(k)=amax1(0.,tot(kp)-water(k))
            enddo
          endif
        enddo

      print*,minval(water),maxval(water)

      return
      end
