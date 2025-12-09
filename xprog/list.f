      program list
c     ============
c
c     program to extract line data for a particular isotopologue
c     of a molecule from the EXOMOL data and output them in the
c     Synspec format
c
c     Standard Input: 
C      type   - a Kurucz-type label of a molecule, e.g 822.00 for TiO
c      molname - a character variable with the chemical label, e.g. 'TiO'
c      iat1,iso1 - atomic number and isotope of the first atom.
c                  e.g  22, 46 for 46Ti
c      iat2, iso2 - the same for the second atom, e.g.  8  16 for 16O
c      core - character variable with the core name of the EXOMOL files,
c             e.g. '46Ti-16O__Toto'
c      wbeg - lowest wavelength [nm] for a line to be included in the list
c      wend - the same for the longest wavelength
c      gflim - log(gf) limit- lines with log(gf) below gflim are excluded
c     
c     Input files:
c      - the program accumes the the three EXOMOL files are present 
c        in the current directory, e.g.
c        46Ti-16O__Toto.states   - energy states
c        46Ti-16O__Toto.trans    - transitions
c        46Ti-16O__Toto.pf       - partition functions
c      - file "isotops" - solar isotoic fractions of the individual atoms
c      - file with the name [molname].pf, e.g. TiO.pf, located in the
c        directory ../../EXOMOL_pfs (may be changed to suit the user's needs)  
c
c     Output:
c      - standard ouput - a short log file
c      - file fort.10 - the resulting line list in the Synspec format,
c        but in reverse wavelength order  - from the long to short.
c        An accompanied short program "reverse.f" reverses it to the
c        ordering needed by Synspec - short to long wavelengths
c     
c        Notice that wavelengths above 2000A are transformed from the
c        vacuum to air wavelengths.
c        an accompanied short program "reverse.f" reverses it to the
c        ordering needed by Synspec - short to long wavelengths
c
c
      implicit real*8(a-h,o-z)
      parameter (mt=10000,mst=1000000,mato=92,miso=290)
      dimension t(mt),partf(mt)
      dimension fra(mato,miso)

      character*50 core,molname,sfile,tfile,pfile,pafile
      dimension e(mst),g(mst)
c
      con=1.e-14/8./3.1415926/0.02654
      gr=1.6e7
      gs=4.1e-17
      gw=1.e-7
c
      do ia=1,mato
         do is=1,miso
            fra(ia,is)=0.
         end do
      end do
c
      read(5,*) type,molname
      read(5,*) iat1,iso1,iat2,iso2
      read(5,*) core
c     read(5,*) wbeg,wend,gflim,fiso,relu
      read(5,*) wbeg,wend,gflim
      write(*,*) 'species:       ',molname
      write(*,*) 'EXOMOL files  :',core
      write(6,602) iat1,iso1,iat2,iso2
  602 format(' IAT1,ISO1,IAT2,ISO2',2i4,2x,2i4)
      write(6,603) wbeg,wend,gflim
  603 format(' Wbeg,Wend,gflim    ',3f11.2)
c
c     determine the approprite isotopic fraction of the given isotopologue
c
      open(unit=7,file='isotops',status='old')
   10 continue
      read(7,*,end=20) iatm,isot,fr
      fra(iatm,isot)=fr*0.01
      go to 10
   20 continue
      close(7)
c
c     normalize the isotopic fraction
c
      sfr=0.
      do is1=1,miso
         do is2=1,miso
            sfr=sfr+fra(iat1,is1)*fra(iat2,is2)
         end do
      end do
      write(*,*) 'frtot',sfr
c
      fiso=fra(iat1,iso1)*fra(iat2,iso2)
c
      wtbeg=1.d7/max(wbeg,100.)
      wtend=1.d7/wend
c
      do i=len(core),1,-1
         if(core(i:i).ne.' ') go to 30
      end do
   30 continue
      sfile= core(1:i) // '.states'
      tfile= core(1:i) // '.trans'
      pfile= core(1:i) // '.pf'
      write(*,*) 'sfile  ',sfile
      write(*,*) 'tfile  ',tfile
      write(*,*) 'pfile  ',pfile
c
c     read file with partition functions for the given isotopologue
c
      open(unit=7,file=pfile,status='old')
cc    read(7,*) t1,pf1
c     
   35 continue
      read(7,*) t2,pf2
      if(t2.gt.999.9.and.t2.lt.1000.1) go to 37
      go to 35
   37 continue
c
c     read avereaged partition function
c
      do i=len(molname),1,-1
         if(molname(i:i).ne.' ') go to 40
      end do
   40 continue
      im=i
      pafile= molname(1:im) // '.pf'
      open(unit=7,file=pafile,status='old')
cc    read(7,*) ta1,pfa1
   45 continue
      read(7,*) ta2,pfa2
      if(ta2.gt.999.9.and.ta2.lt.1000.1) go to 47
      go to 45
   47 continue
      close(7)
c
c     determine ratio of the partition functions at 1000 K,
c     which is constant with T
c
      relu=1.
      relu1=1.
      relu2=1.
c     relu1=pf1/pfa1
      relu2=pf2/pfa2
      relu=relu2
      fis=fiso
      fiso=log10(fiso/relu)
c
      write(6,601) t2,pf2,ta2,pfa2,relu2
  601 format(' T = ',0pf10.1,1pe11.3,0pf10.1,1pe11.3,0pf10.3/)
      write(6,604) fis,relu,fiso
  604 format(' isot.fract. ',f21.14/
     *       ' U(iso)/U(av)',f21.14/
     *       ' log(cor.fac)',f21.14)   
c           
      open(unit=11,file=sfile,status='old')
      open(unit=12,file=tfile,status='old')
c
c     read data for energy states
c
      i=0
   50 continue
      i=i+1
      read(11,*,err=60,end=60) i,e(i),g(i)
      go to 50
   60 continue
      i=i-1
      write(*,*) 'total number of states',i
      close(11)
c
c     read data for transitions and write line data in the Synspec format
c
      n=0
      ninc=0
   70 continue
      n=n+1
      read(12,*,end=80) nu,nl,a,wtr
c     if(wtr.lt.wtend) go to 70
c     if(wtr.gt.wtbeg) go to 80
      if(e(nu)-e(nl).le.0.) go to 70
      if(e(nu).gt.0..and.e(nu).lt.2.e5.and.
     *   e(nl).gt.0..and.e(nl).lt.2.e5) then
CC    wlam=1.d7/(e(nu)-e(nl))
      wlam=1.d7/wtr
      if(wlam.ge.wend*0.999999) go to 70
      if(wlam.lt.wbeg) go to 70
c
c     transform to air wavelengths above 2000 A
c
c     vac=1.e20
      vac=200.
      wlam0=wlam
      if(wlam.gt.vac) then
         wl0=wlam*10.
         ALM=1.E8/(WL0*WL0)
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
c        if(n.le.2) write(6,641) i,wl0,alm,xn1,wlam
         WL0=WL0/(XN1*1.E-6+1.)
         wlam=wl0*0.1
      END IF
c
      gf=log10(a*wlam0*wlam0*g(nu)*con)+fiso
      if(gf.lt.gflim) go to 70
      ninc=ninc+1
c     write(10,610) wlam,type,gf,e(nl),gr,gs,gw
      write(10,610) wlam,type,gf,e(nl)
c 610 format(f10.4,f10.2,f8.3,f12.3,1p3e10.2)
  610 format(f10.4,f10.2,f8.3,f12.3)
      end if
      go to 70
   80 continue
      close(12)
      close(10)
c
      write(*,*) 'read',n,' transitions; included',ninc
      close(6)
      end
