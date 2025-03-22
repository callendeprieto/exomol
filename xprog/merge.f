      program mergelist
c
      parameter(nl=50000000)
      dimension alam(nl),anum(nl),gf(nl),excl(nl),gr(nl),gs(nl),gw(nl)     
c
      i=1
   10 continue
      read(5,*,end=20) alam(i),anum(i),gf(i),excl(i),gr(i),gs(i),gw(i)
      i=i+1
      go to 10
   20 continue
      n=i-1
      write(*,*) 'number of lines in the first list ',n
c
      j=0
      k=1
      inewend=0
      ioldend=0
   30 continue
      read(10,*,end=50) al,anu,a,ex,r,s,w
      write(*,*)'-->',al,anu,a,ex,r,s,w
      if(ioldend.eq.1) then
         write(11,611) al,anu,a,ex,r,s,w
         j=j+1
         go to 30
      end if
      if(al.le.alam(k)) then
         write(11,611) al,anu,a,ex,r,s,w
         j=j+1
       else
   40    continue
         if(alam(k).lt.al) then 
             write(6,621) k,n,alam(k),al,j,ioldend,inewend
  621        format(2i10,2f15.4,3i10)
            if(k.gt.n) go to 70
            write(11,611) alam(k),anum(k),gf(k),excl(k),gr(k),gs(k),
     *      gw(k)
            k=k+1
            go to 40
          else
            write(11,611) al,anu,a,ex,r,s,w
            j=j+1
         end if
      end if
      go to 30
   50 continue
      inewend=1
   60 continue
      if(k.gt.n) go to 70
         write(11,611) alam(k),anum(k),gf(k),excl(k),gr(k),gs(k),gw(k)
         k=k+1
         go to 60
   70 continue
      ioldend=1
      if(inewend.eq.0) then
         write(11,611) al,anu,a,ex,r,s,w
         j=j+1
         go to 30
      end if
      write(*,*) 'number of lines in the second list',j
  611 format(f10.4,f10.2,f8.3,f12.3,1p3e10.2)
c
      end
 
