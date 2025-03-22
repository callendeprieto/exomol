      program isosum
c
      parameter (mt=10000,mato=26,miso=58,mcomp=6)
      dimension tem(mt,mcomp),pf(mt,mcomp),partf(mt),ntem(mcomp)
      dimension fra(mato,miso),t(mt),fract(mcomp),tin(mt),pin(mt)
      dimension dtem(mcomp),pfint(mt)
      character*30 filepf
c
c     read (solar) isotopic fractions
c
      open(unit=7,file='isotops',status='old')
   10 continue
      read(7,*,end=20) iatm,isot,fr
      fra(iatm,isot)=fr*0.01
      go to 10
   20 continue
      close(7)
c
      do i=1,mtemp
         partf(i)=0.
      end do
      fratot=0.
      read(5,*) iat1,iat2
      iso=0
      isomax=1
      framax=0.
   25 continue
c
c        read partition funcion foor all isotopologues
c
         read(5,*,end=50) iso1,iso2,filepf
         frac=fra(iat1,iso1)*fra(iat2,iso2)
         fratot=fratot+frac
         iso=iso+1
         fract(iso)=frac
         if(frac.ge.framax) then
            framax=frac
            isomax=iso
         end if
         open(unit=7,file=filepf,status='old')
         it=0
   30    continue
         it=it+1
         read(7,*,end=40) tem(it,iso),pf(it,iso)
         go to 30
   40    continue
         ntem(iso)=it-1
         dtem(iso)=tem(2,iso)-tem(1,iso)
         close(7)
         write(6,601) iat1,iso1,iat2,iso2,ntem(iso),dtem(iso),
     *                fra(iat1,iso1),fra(iat2,iso2),frac,fratot
  601    format('atom1',2i5,'    atom2',i7,i5,i6,f6.1,4f8.4)
c        write(*,*) 'isomax,framax,fratot',isomax,framax,fratot
      go to 25
   50 continue
      niso=iso
      do isof=1,niso
         fract(isof)=fract(isof)/fratot
      end do
c
c     the case of just one isotopologue
c
      if(niso.eq.1) then
         do it=1,ntem(iso)
            t(it)=tem(it,iso)
            partf(it)=pf(it,iso)
         end do
       else
c
c        more isotopologues
c
         ntmin=ntem(1)
         ntmax=ntem(1)
         tmin=tem(1,1)
         tmax=tem(ntem(1),1)
         do iso=2,niso
            ntmin=min(ntmin,ntem(iso))
            ntmax=max(ntmax,ntem(iso))
            tmin=min(tmin,tem(1,iso))
            tmax=max(tmax,tem(ntem(iso),iso))
         end do
c
c        the case of all the isotop. tables have the same T_min and T_max
c        and also the same numer of temperatures - no interpolation
c
         if((tmax-tmin)/tmin.lt.0.0001.and.ntmax.eq.ntmin) then
            do it=1,ntmax
               t(it)=tem(it,1)
               do iso=1,niso
                  partf(it)=partf(it)+pf(it,iso)*fract(iso)
               end do
            end do
          else
c
c           other cses - for simplicity, one interpolates even if not nercessary
c
            do it=1,ntem(isomax)
               t(it)=tem(it,isomax)
               partf(it)=pf(it,isomax)
c              pf0(it)=partf(it)
            end do
            do iso=1,niso
               if(iso.ne.isomax) then
                  do it=1,ntem(iso)
                     tin(it)=tem(it,iso)
                     pin(it)=pf(it,iso)
                  end do
                  call interp(tin,pin,t,pfint,ntem(iso),ntem(isomax),
     *                        2,0,0)
                  do it=1,ntem(isomax)
                     if(t(i).lt.tin(1)) pfint(i)=0.
                     if(t(i).gt.tin(ntem(iso))) pfint(i)=0.
                     partf(it)=partf(it)+pfint(it)*fract(iso)
                  end do
               end if
            end do   
         end if
      end if
c
      do it=1,ntem(isomax)
         write(10,610) t(it),partf(it)
      end do
  610 format(f8.1,f16.4)
      close(10)
c
      end
      
c----------------


      SUBROUTINE INTERP(X,Y,XX,YY,NX,NXX,NPOL,ILOGX,ILOGY)
C     ====================================================
C
C     General interpolation procedure of the (NPOL-1)-th order
C
C     for  ILOGX = 1  logarithmic interpolation in X
C     for  ILOGY = 1  logarithmic interpolation in Y
C
C     Input:
C      X    - array of original x-coordinates
C      Y    - array of corresponding functional values Y=y(X)
C      NX   - number of elements in arrays X or Y
C      XX   - array of new x-coordinates (to which is to be
C             interpolated
C      NXX  - number of elements in array XX
C     Output:
C      YY   - interpolated functional values YY=y(XX)
C
      parameter(mx=10000)
      DIMENSION X(MX),Y(MX),XX(MX),YY(MX)
C
C     no interpolation for NPOL.LE.0 or NX.le.0
C
      IF(NPOL.LE.0.OR.NX.LE.0) THEN
         N=NX
         IF(NXX.GE.NX) N=NXX
         DO I=1,N
            XX(I)=X(I)
            YY(I)=Y(I)
         END DO
         RETURN
      END IF
C
C     interpolation
C
C     if required, compute logarithms to be interpolated
C
      IF(ILOGX.GT.0) THEN
      DO I=1,NX
         X(I)=LOG10(X(I))
      END DO
      DO I=1,NXX
         XX(I)=LOG10(XX(I))
      END DO
      END IF
      IF(ILOGY.GT.0) THEN
      DO I=1,NX
         Y(I)=LOG10(Y(I))
      END DO
      END IF
C
      NM=(NPOL+1)/2
      NM1=NM+1
      NUP=NX+NM1-NPOL
      DO ID=1,NXX
         XXX=XX(ID)
         DO I=NM1,NUP
            IF(XXX.LE.X(I)) GO TO 10
         END DO
         I=NUP
   10    J=I-NM
         JJ=J+NPOL-1
         YYY=0.
         DO K=J,JJ
            T=1.
            DO 20 M=J,JJ
               IF(K.EQ.M) GO TO 20
               T=T*(XXX-X(M))/(X(K)-X(M))
   20       CONTINUE
            YYY=Y(K)*T+YYY
         END DO
         YY(ID)=YYY
      END DO
C
      IF(ILOGX.GT.0) THEN
      DO I=1,NX
         X(I)=EXP(X(I)*2.30258509299405)
      END DO
      DO I=1,NXX
         XX(I)=EXP(XX(I)*2.30258509299405)
      END DO
      END IF
      IF(ILOGY.GT.0) THEN
      DO I=1,NX
         Y(I)=EXP(Y(I)*2.30258509299405)
      END DO
      DO I=1,NXX
         YY(I)=EXP(YY(I)*2.30258509299405)
      END DO
      END IF
C
      RETURN
      END

