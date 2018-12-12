      subroutine find_neighbour_list (nb,nn,nnode,nbmax,x,y,sides,
     &                                nt,vertices,neighbour)

c subroutine to calculate the new neighbour list
c as well as the length of the sides of the Voronoi cells

c INPUT: x,y       = x,y nodal coordinates
c        nnode     = number of nodes
c        nbmax     = maximum number of neighbour per node
c        vertices  = triangle list
c        neighbour = neighbour list
c        nt        = number of triangle

c OUTPUT: nn        = neighbour array
c         nb        = number of neighbour for each node
c         sides     = length of the sides of the voronoi cells

c subroutines called:
c find_centre

      integer     nb(nnode),nn(nbmax,nnode)
      real        sides(nbmax,nnode)
      integer     vertices(3,nt),neighbour(3,nt)
      real        x(nnode),y(nnode)

        do i=1,nnode
        nb(i)=0
        enddo

        do it=1,nt
        i1=vertices(1,it)
        i2=vertices(2,it)
        i3=vertices(3,it)
        call find_centre (xa,ya,x(i1),y(i1),x(i2),y(i2),x(i3),y(i3))
        nb(i1)=nb(i1)+1
          if (nb(i1).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i1),i1)=i2
        itp=neighbour(3,it)
          if (itp.eq.0) then
          nb(i2)=nb(i2)+1
            if (nb(i2).gt.nbmax) then
            stop 'nbmax too small...'
            endif
          nn(nb(i2),i2)=i1
          rad=(xa-x(i1))**2+(ya-y(i1))**2
          xln=((x(i1)-x(i2))**2+(y(i1)-y(i2))**2)/4.
          side=2.*sqrt(rad-xln)
          sides(nb(i2),i2)=side
          sides(nb(i1),i1)=side
          else
          j1=vertices(1,itp)
          j2=vertices(2,itp)
          j3=vertices(3,itp)
          call find_centre (xb,yb,x(j1),y(j1),x(j2),y(j2),x(j3),y(j3))
          sides(nb(i1),i1)=sqrt((xa-xb)**2+(ya-yb)**2)
          endif
        nb(i2)=nb(i2)+1
          if (nb(i2).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i2),i2)=i3
        itp=neighbour(1,it)
          if (itp.eq.0) then
          nb(i3)=nb(i3)+1
            if (nb(i3).gt.nbmax) then
            stop 'nbmax too small...'
            endif
          nn(nb(i3),i3)=i2
          rad=(xa-x(i2))**2+(ya-y(i2))**2
          xln=((x(i2)-x(i3))**2+(y(i2)-y(i3))**2)/4.
          side=2.*sqrt(rad-xln)
          sides(nb(i3),i3)=side
          sides(nb(i2),i2)=side
          else
          j1=vertices(1,itp)
          j2=vertices(2,itp)
          j3=vertices(3,itp)
          call find_centre (xb,yb,x(j1),y(j1),x(j2),y(j2),x(j3),y(j3))
          sides(nb(i2),i2)=sqrt((xa-xb)**2+(ya-yb)**2)
          endif
        nb(i3)=nb(i3)+1
          if (nb(i3).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i3),i3)=i1
        itp=neighbour(2,it)
          if (itp.eq.0) then
          nb(i1)=nb(i1)+1
            if (nb(i1).gt.nbmax) then
            stop 'nbmax too small...'
            endif
          nn(nb(i1),i1)=i3
          rad=(xa-x(i3))**2+(ya-y(i3))**2
          xln=((x(i3)-x(i1))**2+(y(i3)-y(i1))**2)/4.
          side=2.*sqrt(rad-xln)
          sides(nb(i1),i1)=side
          sides(nb(i3),i3)=side
          else
          j1=vertices(1,itp)
          j2=vertices(2,itp)
          j3=vertices(3,itp)
          call find_centre (xb,yb,x(j1),y(j1),x(j2),y(j2),x(j3),y(j3))
          sides(nb(i3),i3)=sqrt((xa-xb)**2+(ya-yb)**2)
          endif
        enddo

      return
      end
