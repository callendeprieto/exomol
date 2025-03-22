      program reverse
c
      character*80 a
      dimension a(70000000)
c
      i=0
   10 continue
      i=i+1
      read(5,501,end=20) a(i)
      go to 10
   20 continue
  501 format(a80)
      n=i-1
c
      do i=n,1,-1
         write(6,501) a(i)
      end do
c
      end
