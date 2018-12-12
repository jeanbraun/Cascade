c---
      subroutine find_centre (xc,yc,x1,y1,x2,y2,x3,y3)

c subroutine to find the centre of the circumcircle defined
c about the three vertices of a triangle

c INPUT: x1,x2,x3 = x-coordinates of the thrre vertices
c        y1,y2,y3 = y-coordinates of the thrre vertices

c OUTPUT: xc,yc x- and y- coordinates of the centre

c subroutines called:
c none

      dx2m1 = x2-x1
      dx2p1 = x2+x1
      dy2m1 = y2-y1
      dy2p1 = y2+y1
      dx3m1 = x3-x1
      dx3p1 = x3+x1
      dy3m1 = y3-y1
      dy3p1 = y3+y1
      denom = dx2m1*dy3m1-dx3m1*dy2m1
      xc = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &     -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &     (denom)

      yc = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0
     &     -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/
     &     (denom)

      return
      end

