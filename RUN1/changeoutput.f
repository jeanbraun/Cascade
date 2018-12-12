      program ch

      open (7,file='outflux',status='old')
      open (8,file='out.txt',status='unknown')
      do i=1,1000
      read (7,*) x,y,z
      write (8,*) x,y
      enddo
      close (7)
      close (8)

      end
