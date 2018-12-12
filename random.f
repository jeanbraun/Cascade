c random

      subroutine random (x,n)

c subroutine used to generate random numbers (derived from
c numerical recipes; see book for further information)

c INPUT: n    = number of random numbers to be generated

c OUTPUT: x   = array of random numbers of length n

c subroutines called:
c NONE

      common /vocal/ ivocal

      real*4    x(n)

      idum=-1
      xx=ran0(idum1)

      idum=1
          do 100 i=1,n
          x(i)=ran0(idum)
100       continue

      return
      end

c---
      FUNCTION RAN(ISEED)
      PARAMETER(IA=7141,IC=54773,IM=259200)
      ISEED=MOD(ISEED*IA+IC,IM)
      RAN=FLOAT(ISEED)/FLOAT(IM)
      RETURN
      END

c---
      FUNCTION RAN0(IDUM)
      DIMENSION V(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        ISEED=ABS(IDUM)
        IDUM=1
        DO 11 J=1,97
          DUM=RAN(ISEED)
11      CONTINUE
        DO 12 J=1,97
          V(J)=RAN(ISEED)
12      CONTINUE
        Y=RAN(ISEED)
      ENDIF
      J=1+INT(97.*Y)
      IF(J.GT.97.OR.J.LT.1)PAUSE
      Y=V(J)
      RAN0=Y
      V(J)=RAN(ISEED)
      RETURN
      END
