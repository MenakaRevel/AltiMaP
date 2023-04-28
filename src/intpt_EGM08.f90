!==============================================================================
! Original script was subjected to their own license
! This script is modified by Menaka Revel@IIS and is subjected to the same license
! Credit goes to the original authors
!==============================================================================
program int_main
    !      
    !-----------------------------------------------------------------------
    !
    !     INTERP IS A F77 PROGRAM FOR INTERPOLATING FROM A GLOBAL GRID OF
    !     REGULARLY SPACED POINT VALUES TO SELECTED SCATTERED COORDINATES.
    !     INTERP REQUIRES AN INPUT FILE CONTAINING THE GLOBALLY GRIDDED
    !     VALUES, AND AN INPUT FILE CONTAINING THE GEOGRAPHICAL COORDINATES 
    !     OF THE SCATTERED POINTS FOR WHICH INTERPLATED VALUES ARE REQUIRED.
    !
    !     THE INPUT FILE OF GLOBALLY GRIDDED POINT VALUES WILL STORE THESE  
    !     VALUES IN REAL*4 SEQUENTIAL BINARY FORMAT, WITH ONE RECORD FOR 
    !     EACH PARALLEL BAND. THESE INPUT VALUES ARE EQUALLY SPACED IN 
    !     LATITUDE AND LONGITUDE, AND ARE SITUATED AT THE CORNERS OF THEIR 
    !     RESPECTIVE CELLS, SUCH THAT THE TOP-LEFT POINT IN THE GRID HAS A 
    !     LONGITUDE OF ZERO DEGREES AND A LATITUDE OF NINETY DEGREES. THE 
    !     FIRST RECORD CONTAINS THE NORTHERN-MOST PARALLEL, AND THE FIRST 
    !     VALUE IN EACH RECORD IS THE WESTERNMOST VALUE FOR THAT PARALLEL. 
    !     NOTE THAT GRID VALUES SITUATED ON THE ZERO MERIDIAN APPEAR ONLY
    !     ONCE, AS THE FIRST VALUE IN THEIR RESPECTIVE RECORD AT ZERO 
    !     LONGITUDE. THESE VALUES ARE NOT REPEATED AT THE END OF THEIR 
    !     RESPECTIVE RECORDS AT LONGITUDE = 360 DEGREEES. 
    !
    !     THE ASCII INPUT FILE CONTAINING THE GEOGRAPHICAL COORDINATES OF
    !     THE SCATTERED POINTS CONTAINS ONE RECORD FOR EACH POINT. EACH 
    !     RECORD CONTAINS THE GEODETIC LATITUDE AND THEN THE LONGITUDE 
    !     IN DECIMAL DEGREES FOR ITS RESPECTIVE POINT. 
    !
    !     THE ASCII OUTPUT FILE CONTAINING THE INTERPOLATED VALUES FOR THE
    !     SCATTERED POINTS ALSO CONTAINS ONE RECORD FOR EACH POINT. EACH 
    !     RECORD CONTAINS THE GEODETIC LATITUDE AND LONGITUDE IN DECIMAL 
    !     DEGREES, FOLLOWED BY THE INTERPOLATED VALUE FOR THAT POINT.
    !
    !     THE INPUT PARAMETERS FOR INTERP ARE:
    !
    !     path_gd  :  PATH TO THE INPUT FILE OF GLOBALLY GRIDDED VALUES, 
    !     name_gd  :  NAME OF THE INPUT FILE OF GLOBALLY GRIDDED VALUES,
    !     path_pt  :  PATH TO THE INPUT FILE OF INTERPOLATION LATS AND LONS,  
    !     name_pt  :  NAME OF THE INPUT FILE OF INTERPOLATION LATS AND LONS,
    !     path_out :  PATH TO THE OUTPUT FILE OF INTERPOLATED VALUES,  
    !     name_out :  NAME OF THE OUTPUT FILE OF INTERPOLATED VALUES,
    !
    !     dlat     :  LAT SPACING (DEGREES) OF INPUT GRID,
    !     dlon     :  LON SPACING (DEGREES) OF INPUT GRID.
    !      
    !-----------------------------------------------------------------------
    !     ORIGINAL PROGRAM:                           NIKOS PAVLIS
    !     MODIFIED FOR USE BY NGA WITH 30 SEC GRIDS   SIMON HOLMES, AUG 2007
    !     CORNER-CELL REGISTRATION                    SIMON HOLMES, MAY 2008 
    !-----------------------------------------------------------------------
          implicit real*8(a-h,o-z)
          character*120 path_gd,p_gd,name_gd,n_gd,fnu1
          character*120 path_pt,p_pt,name_pt,n_pt,fnu2
          character*120 path_out,pout,name_out,nout,fnu10
    !-----------------------------------------------------------------------
    !
    !     INPUT PARAMETERS
    !
    !-----------------------------------------------------------------------
          parameter( &
         &   path_gd  = './src',                                           &
         &   name_gd  = '/Und_min1x1_egm2008_isw=82_WGS84_TideFree_SE',  &
         &   path_pt  = './',                                           &
         &   name_pt  = 'INPUT.DAT',                                    &
         &   path_out = './',                                           &
         &   name_out = 'OUTPUT_EGM08.DAT')
    !-----------------------------------------------------------------------
          parameter(dlat     = 1.d0/60.d0, &
         &          dlon     = 1.d0/60.d0)
    !-----------------------------------------------------------------------
    !
    !     NON-INPUT PARAMETERS
    !
    !-----------------------------------------------------------------------
          parameter(iwind    = 6,           &
         &          iw       = iwind+1)     
          parameter(nrows    = 10801,       &
         &          ncols    = 21600,       &
         &          nriw2    = nrows+2*iw,  &
         &          nciw2    = ncols+2*iw)
    !      
          real*8 scrd(ncols,2),statd(22),grid(nriw2,nciw2)
          real*4 temp(nciw2) 
    !-----------------------------------------------------------------------
          write(6,5000)
          write(6,400)
      400 format(30x,'Execution')
          write(6,5000)
    !-----------------------------------------------------------------------
    !
    !     OPEN INPUT AND OUTPUT FILES
    !
    !-----------------------------------------------------------------------
          call extract_name_120(path_gd,p_gd,nchr_pin)
          call extract_name_120(name_gd,n_gd,nchr_nin)
          fnu1=p_gd(1:nchr_pin)//n_gd(1:nchr_nin)
    !
          open(1,file=trim(fnu1),form='unformatted', &
          &     status='old',iostat=ios)
    !
          call extract_name_120(path_pt,p_pt,nchr_ppt)
          call extract_name_120(name_pt,n_pt,nchr_npt)
          fnu2=p_pt(1:nchr_ppt)//n_pt(1:nchr_npt)
    !
          open(2,file=trim(fnu2),form='formatted',status='old',iostat=ios)
    !
          call extract_name_120(path_out,pout,nchr_pot)
          call extract_name_120(name_out,nout,nchr_not)
          fnu10=pout(1:nchr_pot)//nout(1:nchr_not)
    !
          open(10,file=trim(fnu10),form='formatted',iostat=ios)
    !
          write(6,410) n_gd
      410 format(' Input file containing globally gridded values : ', &
         &                                                          a120,//)
          write(6,415) n_pt
      415 format(' Input file file containing lats and lons of ', &
         &                               'interpolation points : ' ,a120,//)
          write(6,420) nout
      420 format(' Output file containing lats, lons and interpolated ', &
         &                                                 'values : ',a120)
    !-----------------------------------------------------------------------
    !
    !     READ INPUT GRID AND COMPUTE GLOBAL STATS
    !
    !-----------------------------------------------------------------------
          write(6,5000)
          write(6,700) name_gd
      700 format(10x,'Statistics of Input Grid: ',a120,//)
    !-----------------------------------------------------------------------
          do i = 1, nrows
            ii = nrows+1-i
            read(1) (temp(j+iw), j=1,ncols)
            do j = 1, ncols
              grid(ii+iw,j+iw) = temp(j+iw)
              scrd(j,1) = grid(ii+iw,j+iw)
              scrd(j,2) = 1.d0
            enddo  !  j
            call stats(90.d0,0.d0,i,dlat,dlon,nrows,ncols,scrd,9999.d0,0, &
         &                                                      statd,icell)    
          enddo  !  i  
          close(1)
    !-----------------------------------------------------------------------
    !
    !     PAD EDGES OF INPUT GRID
    !
    !-----------------------------------------------------------------------
          do i = iw+1, nrows+iw
            do j = 1, iw
              grid(i,j) = grid(i,j+ncols)
            enddo  !  j
            do j = ncols+iw+1, ncols+2*iw
              grid(i,j) = grid(i,j-ncols)
            enddo  !  j
          enddo  !  i
    !
          do i = 1, iw
            ii = 2*iw+1-i
            do j = 1, ncols+2*iw
              jj = j + ncols/2
              if (jj.gt.ncols+2*iw) jj = jj-ncols
              grid(i,j) = grid(ii,jj)
            enddo  !  j
          enddo  !  i
    !
          do i = nrows+iw+1, nrows+2*iw
            ii = 2*(nrows+iw) + 1 - i
           do j = 1, ncols+2*iw
              jj = j + ncols/2
              if (jj.gt.ncols+2*iw) jj = jj-ncols
              grid(i,j) = grid(ii,jj)
            enddo  !  j
          enddo  !  i
    !-----------------------------------------------------------------------
    !
    !     INTERPOLATE
    !
    !-----------------------------------------------------------------------
          dmin =  0.d0
          slat = -90.d0 - dlat*iw  !  (lat < -90) OK
          wlon =        - dlon*iw  !  (lon <   0) OK
    !
          npt = 0
    !
     1000 read(2,*,end=2000) flat,flon
            sav_lon = flon
    !
            npt = npt + 1
            if (flon.gt.360.d0) flon = flon-360.d0
            if (flon.lt.  0.d0) flon = flon+360.d0 
    !-----------------------------------------------------------------------
    !
    !       COORDS OK?
    !
    !-----------------------------------------------------------------------
            if (flat.gt.90.d0.or.flat.lt.-90.d0.or. &
         &      flon.gt.360.d0.or.flon.lt.0.d0) then
              write(10,1050) flat,flon
              goto 1000
            endif  !  flat;flon
     1050   format(2(1x,f11.6),2x,'INVALID')    
    !-----------------------------------------------------------------------
            call interp(iwind,dmin,grid,slat,wlon,dlat,dlon,nriw2,nciw2, &
         &              nriw2,nciw2,flat,flon,val)
    !
            write(10,1100) flat,sav_lon,val
     1100   format(2(1x,f11.6),2x,f9.3)         
            goto 1000
    !
     2000 continue
          close(10)
    !-----------------------------------------------------------------------
    !
    !     DONE
    !
    !-----------------------------------------------------------------------
          write(6,5000)
          write(6,2100) npt
          write(6,2300)
          write(6,5000)
    !
     2100 format(10x,'Number of Points Interpolated     = ',i12,/)
     2300  format(27x,'Normal Termination')
    !           
     5000 format(//,'!',71('-'),//)
    !
          stop
          end
    !
    !-----------------------------------------------------------------------
    !
          SUBROUTINE STATS(TOPLAT,WSTLON,I,GRDN,GRDE,is,NCOLS,DATA, &
          &                 EXCLUD,ISIG,STAT,icent)
    !-----------------------------------------------------------------------
          IMPLICIT REAL*8(A-H,O-Z)
          CHARACTER*20 DLABEL(14)
          CHARACTER*20 SLABEL( 8)
          DIMENSION DATA(NCOLS,2)
          DOUBLE PRECISION DEXCLUD,DG,SD,STAT(22)
          SAVE
          data ix/0/
          DATA PI/3.14159265358979323846D+00/
          DATA DLABEL/'    Number of Values','  Percentage of Area', &
         &            '       Minimum Value',' Latitude of Minimum', &
         &            'Longitude of Minimum','       Maximum Value', &
         &            ' Latitude of Maximum','Longitude of Maximum', &
         &            '     Arithmetic Mean','  Area-Weighted Mean', &
         &            '      Arithmetic RMS','   Area-Weighted RMS', &
         &            '   Arithmetic S.Dev.','Area-Weighted S.Dev.'/
          DATA SLABEL/'       Minimum Sigma',' Latitude of Minimum', &
         &            'Longitude of Minimum','       Maximum Sigma', &
         &            ' Latitude of Maximum','Longitude of Maximum', &
         &            'Arithmetic RMS Sigma','Area-wghtd RMS Sigma'/
    !-----------------------------------------------------------------------
          IF(I.EQ.1) THEN
          DEXCLUD=EXCLUD
          DTR=PI/180.D0
          FOURPI=4.D0*PI
          DPR=GRDN*DTR
          DLR=GRDE*DTR
          CAREA=2.D0*DLR*SIN(DPR/2.D0)
          DO 10 K=1,22
          STAT(K)=0.D0
       10 CONTINUE
          STAT( 3)= DEXCLUD
          STAT( 6)=-DEXCLUD
          STAT(15)= DEXCLUD
          STAT(18)= 0.D0
          ENDIF
    !-----------------------------------------------------------------------
    !
          dlat = toplat -(i-1.d0)*grdn - (grdn/2.d0)*icent
          COLATC=(90.D0-DLAT)*DTR
          AREA=CAREA*SIN(COLATC)
    !
    !   The following statement avoids over-estimating the total area
    !   covered by the grid (of POINT values) when the pole(s) are part
    !   of the grid.
    !
          if(abs(dlat).eq.90.d0) area=2.d0*dlr*sin(dpr/4.d0)**2
    !
    !-----------------------------------------------------------------------
          DO 110 J=1,NCOLS
          DLON=WSTLON+(J-1.D0)*GRDE +(GRDE/2.D0)*icent
          DG=DATA(J,1)
          SD=DATA(J,2)
          IF(DG.LT.DEXCLUD) THEN
    !-----------------------------------------------------------------------
          STAT( 1)=STAT( 1)+1.D0
          STAT( 2)=STAT( 2)+AREA
          IF(DG.LE.STAT( 3)) THEN
          STAT( 3)=DG
          STAT( 4)=DLAT
          STAT( 5)=DLON
          ENDIF
          IF(DG.GE.STAT( 6)) THEN
          STAT( 6)=DG
          STAT( 7)=DLAT
          STAT( 8)=DLON
          ENDIF
          STAT( 9)=STAT( 9)+DG
          STAT(10)=STAT(10)+DG*AREA
          STAT(11)=STAT(11)+DG**2
          STAT(12)=STAT(12)+DG**2*AREA
          IF(SD.LE.STAT(15)) THEN
          STAT(15)=SD
          STAT(16)=DLAT
          STAT(17)=DLON
          ENDIF
          IF(SD.GE.STAT(18)) THEN
          STAT(18)=SD
          STAT(19)=DLAT
          STAT(20)=DLON
          ENDIF
          STAT(21)=STAT(21)+SD**2
          STAT(22)=STAT(22)+SD**2*AREA
    !-----------------------------------------------------------------------
          ENDIF
      110 CONTINUE
    !-----------------------------------------------------------------------
          IF(I.NE.is) RETURN
          IF(STAT(1).GT.0.D0) THEN
          STAT( 9)=STAT( 9)/STAT( 1)
          STAT(10)=STAT(10)/STAT( 2)
          STAT(11)=SQRT(STAT(11)/STAT( 1))
          STAT(12)=SQRT(STAT(12)/STAT( 2))
          STAT(13)=SQRT(STAT(11)**2-STAT( 9)**2)
          STAT(14)=SQRT(STAT(12)**2-STAT(10)**2)
          STAT(21)=SQRT(STAT(21)/STAT( 1))
          STAT(22)=SQRT(STAT(22)/STAT( 2))
          STAT( 2)=STAT( 2)/FOURPI*100.D0
          ELSE
          DO 120 J=3,22
          STAT(J)=DEXCLUD
      120 CONTINUE
          ENDIF
    !=======================================================================
          NUM=INT(STAT(1))
          WRITE(6,6001) DLABEL(1),NUM
     6001 FORMAT(5X,A20,3X,I11)
          DO 210 K=2,14
          WRITE(6,6002) DLABEL(K),STAT(K)
     6002 FORMAT(5X,A20,3X,F15.3)
      210 CONTINUE
          WRITE(6,6003)
     6003 FORMAT(' ')
          IF(ISIG.EQ.1) THEN
          DO 220 K=1,8
          WRITE(6,6002) SLABEL(K),STAT(K+14)
      220 CONTINUE
          ENDIF
    !=======================================================================
          RETURN
          END
    !
    !-----------------------------------------------------------------------
    !
          subroutine extract_name_120(old_name,name,n)
          implicit none
          integer*4 i,n
          character old_name*120,name*120,ch*1
    !
          n = 0
          do i = 1, 120
            ch = old_name(i:i)
            if (ch .ne. ' ') then
              n = n + 1
              name(n:n) = ch
            endif  !  ch
          enddo  !  i
    !
          return
          end
    !
    !-----------------------------------------------------------------------
    !
          SUBROUTINE INTERP(iwO,DMIN,H,PHIS,DLAW,DDFI,DDLA,NPHI,NDLA, &
          &                  IPDIM,ILDIM,PHI,DLA,VALINT)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                      !
    !     SUBROUTINE FOR INTERPOLATION OF VALUES FROM A STANDARD DTM-GRID  !
    !     TO INDIVIDUAL STATION LOCATIONS.                                 !
    !                                                                      !
    !                                                                      !
    !     INPUT PARAMETERS...                                              !
    !     ===================                                              !
    !     iwO...    A SPLINE WINDOW OF SIZE 'iwO' X 'iwO' WILL BE !
    !                  USED AROUND EACH STATION. IF 'iwO' IS 0 OR 1,    !
    !                  BILINEAR INTERPOLATION WILL BE USED.                !
    !     DMIN...      MINIMUM ACCEPTABLE DISTANCE FROM THE GRID EDGE IN   !
    !                  KM (USEFUL FOR FFT GRIDS).                          !
    !     H...         2D DATA ARRAY (ELEMENT (1,1) IN SW CORNER).         !
    !     PHIS,DLAW... LATITUDE AND LONGITUDE OF SW GRID POINT.            !
    !     DDFI,DDLA... GRID SPACING IN LATITUDE AND LONGITUDE DIRECTION.   !
    !     NPHI,NDLA... NUMBER OF GRID POINTS IN LATITUDE AND LONGITUDE     !
    !                  DIRECTION.                                          !
    !     IPDIM,ILDIM..DIMENSIONS OF 2D DATA ARRAY 'H' AS DECLARED IN THE  !
    !                  CALLING PROGRAM.                                    !
    !     PHI,DLA...   LATITUDE AND LONGITUDE OF INTERPOLATION POINT.      !
    !                                                                      !
    !                                                                      !
    !     OUTPUT PARAMETERS...                                             !
    !     ====================                                             !
    !     VALINT...    INTERPOLATED VALUE.                                 !
    !                                                                      !
    !                                                                      !
    !     EXECUTION TIME ON CDC 990 IS...                                  !
    !     ===============================                                  !
    !     +------------------+-------------------+-------------------+     !
    !     I  INTERPOLATION   I  OPT=LOW          I  OPT=HIGH         I     !
    !     I------------------I-------------------I-------------------I     !
    !     I  BILINEAR        I  1.44 MSEC/STAT.  I  1.44 MSEC/STAT.  I     !
    !     I  3 X 3 SPLINE    I  1.53 MSEC/STAT.  I  1.51 MSEC/STAT.  I     !
    !     I  5 X 5 SPLINE    I  1.70 MSEC/STAT.  I  1.67 MSEC/STAT.  I     !
    !     I  7 X 7 SPLINE    I  2.02 MSEC/STAT.  I  1.74 MSEC/STAT.  I     !
    !     I  9 X 9 SPLINE    I  2.31 MSEC/STAT.  I  2.00 MSEC/STAT.  I     !
    !     +------------------+-------------------+-------------------+     !
    !                                                                      !
    !                                                                      !
    !     PROGRAM CREATION BY...   H. DENKER          MAY 30, 1987         !
    !                              H. DENKER          MARCH 13, 1989       !
    !                                                                      !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          PARAMETER(IPA1=20)
          IMPLICIT REAL*8(A-H,O-Z)
          LOGICAL LODD
          REAL*8 H(IPDIM,ILDIM)
          DIMENSION A(IPA1),R(IPA1),Q(IPA1),HC(IPA1)
          IDIM1=IPA1
           TWOPI=6.28318530717959D0
             RHO=360.D0/TWOPI
          REARTH=6371000.D0
          IF(iwO.LT.2) iwO=2
          IF(iwO.GT.IDIM1) iwO=IDIM1
          ILIM=DMIN*1000.*RHO/(REARTH*DDFI)
          JLIM=DMIN*1000.*RHO/(REARTH*DDLA*COS((PHIS+DDFI*NPHI/2.)/RHO))
          LODD=(iwO/2)*2.NE.iwO
          RI=(PHI-PHIS)/DDFI
          RJ=(DLA-DLAW)/DDLA
          IF(LODD) THEN
            I0=RI-0.5
            J0=RJ-0.5
          ELSE
            I0=RI
            J0=RJ
          ENDIF
          I0=I0-iwO/2+1
          J0=J0-iwO/2+1
          II=I0+iwO-1
          JJ=J0+iwO-1
          IF(I0.LT.0 .OR. II.GE.NPHI .OR. J0.LT.0 .OR. JJ.GE.NDLA) THEN
            WRITE(6,7008) PHI,DLA
            VALINT=999999.
            RETURN
          ELSEIF(I0.LT.ILIM .OR. II.GT.NPHI-ILIM .OR. J0.LT.JLIM .OR. &
          &  JJ.GT.NDLA-JLIM) THEN
            IF(NPOINT.LE.ILIST) WRITE(6,7009) PHI,DLA
            VALINT=999999.
            RETURN
          ENDIF
    7008  FORMAT(' ',2F12.6,' STATION TOO NEAR GRID BOUNDARY  - NO INT.' &
         &,' POSSIBLE|')
    7009  FORMAT(' ',2F12.6,' STATION OUTSIDE ACCEPTABLE AREA - NO INT.' &
         &,' PERFORMED|')
          IF(iwO.GT.2) THEN
            DO 110 I=1,iwO
              DO 111 J=1,iwO
                A(J)=H(I0+I,J0+J)
    111       CONTINUE
              CALL INITSP(A,iwO,R,Q)
              HC(I)=SPLINE(RJ-J0+1.,A,iwO,R)
    110     CONTINUE
            CALL INITSP(HC,iwO,R,Q)
            VALINT=SPLINE(RI-I0+1.,HC,iwO,R)
          ELSE
            VALINT=BILIN(RI+1.,RJ+1.,H,NPHI,NDLA,IPDIM,ILDIM)
          ENDIF
          RETURN
          END
    !
    !-----------------------------------------------------------------------
    !
          FUNCTION BILIN(RI,RJ,A,IMAX,JMAX,IADIM,JADIM)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                      !
    !                           B I L I N                                  !
    !                                                                      !
    !  INTERPOLATES VALUES IN AN ARRAY A USING BILINEAR                    !
    !  (PARABOLIC HYPERBOLOID) INTERPOLATION.                              !
    !                                                                      !
    !----------------------------------------------------------------------!
    !                                                                      !
    !  PARAMETERS:                                                         !
    !                                                                      !
    !  BILIN...       INTERPOLATED VALUE                                   !
    !                                                                      !
    !  RI, RJ...      INTERPOLATION ARGUMENT, (1,1) IN LOWER LEFT CORNER,  !
    !                 (IMAX, JMAX) IN UPPER RIGHT.                         !
    !                                                                      !
    !  A...           INTEGER*2 ARRAY WITH ARGUMENTS                       !
    !                                                                      !
    !  IMAX, JMAX...  NUMBER OF POINTS IN GRID                             !
    !                                                                      !
    !  IADIM, JADIM...DECLARED DIMENSIONS OF 'A'                           !
    !                                                                      !
    !  OUTSIDE AREA COVERED BY 'A' THE FUNCTION RETURNS THE VALUE OF       !
    !  THE NEAREST BOUNDARY POINT.                                         !
    !                                                                      !
    !----------------------------------------------------------------------!
    !                                                                      !
    !  PROGRAMMER:                                                         !
    !  RENE FORSBERG, JULY 1983                                            !
    !                                                                      !
    !  MODIFICATIONS BY:                                                   !
    !  HEINER DENKER, 07/01/1987                                           !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IMPLICIT REAL*8(A-H,O-Z)
          DIMENSION A(IADIM, JADIM)
          IN = IFRAC(RI)
          IE = IFRAC(RJ)
          RN = RI - IN
          RE = RJ - IE
          IF (IN.LT.1) THEN
            IN = 1
            RN = 0.0
          ELSEIF (IN.GE.IMAX) THEN
            IN = IMAX-1
            RN = 1.0
          ENDIF
          IF (IE.LT.1) THEN
            IE = 1
            RE = 0.0
          ELSEIF (IE.GE.JMAX) THEN
            IE = JMAX-1
            RE = 1.0
          ENDIF
          RNM1=1.-RN
          REM1=1.-RE
          BILIN = RNM1*REM1*A(IN,IE) + &
          & RN*REM1*A(IN+1,IE) + RNM1*RE*A(IN,IE+1) + &
          & RN*RE*A(IN+1,IE+1)
          RETURN
          END
    !
    !-----------------------------------------------------------------------
    !
          SUBROUTINE INITSP(Y, N, R, Q)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                      !
    !                      I N I T S P                                     !
    !                                                                      !
    !  INITIALIZATION PROCEDURE FOR FAST 1-DIMENSIONAL EQUIDISTANT         !
    !  SPLINE INTERPOLATION, WITH FREE BOUNDARY END CONDITIONS             !
    !  REFERENCE: JOSEF STOER: EINFUHRUNG IN DIE NUMERISCHE MATHEMATIK     !
    !  I, SPRINGER 1972, PAGE 82 AND 86.                                   !
    !                                                                      !
    !----------------------------------------------------------------------!
    !                                                                      !
    !  PARAMETERS (REAL):                                                  !
    !                                                                      !
    !  Y...   GIVEN VALUES, Y(1), ..., Y(N)                                !
    !                                                                      !
    !  R...   SPLINE MOMENTS (1 ... N), TO BE USED BY FUNCTION 'SPLINE'    !
    !                                                                      !
    !  Q...   WORK-ARRAY, DECLARED AT LEAST 1:N                            !
    !                                                                      !
    !----------------------------------------------------------------------!
    !                                                                      !
    !  RENE FORSBERG, JULY 1983                                            !
    !                                                                      !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IMPLICIT REAL*8(A-H,O-Z)
          DIMENSION Y(N), R(N), Q(N)
          Q(1) = 0.0
          R(1) = 0.0
          DO 11 K = 2, N-1
            P = Q(K-1)/2+2
            Q(K) = -0.5/P
            R(K) = (3*(Y(K+1)-2*Y(K)+Y(K-1)) - R(K-1)/2)/P
       11 CONTINUE
          R(N) = 0.0
          DO 12 K = N-1, 2, -1
            R(K) = Q(K)*R(K+1)+R(K)
       12 CONTINUE
          RETURN
          END
    !
    !-----------------------------------------------------------------------
    !
          FUNCTION SPLINE(X, Y, N, R)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                      !
    !                          S P L I N E                                 !
    !                                                                      !
    !  FAST ONE-DIMENSIONAL EQUIDISTANT SPLINE INTERPOLATION FUNCTION.     !
    !  REFERENCE: JOSEF STOER: EINFUHRUNG IN DIE NUMERISCHE MATHEMATIK     !
    !  I, SPRINGER 1972, PAGE 81.                                          !
    !                                                                      !
    !----------------------------------------------------------------------!
    !                                                                      !
    !  PARAMETERS:                                                         !
    !                                                                      !
    !  X...  INTERPOLATION ARGUMENT (REAL), X = 1 FIRST DATA-POINT,        !
    !        X = N LAST DATA-POINT. OUTSIDE THE RANGE LINEAR EXTRA-        !
    !        POLATION IS USED.                                             !
    !                                                                      !
    !  Y...  REAL*8 ARRAY, 1 .. N : DATA VALUES                            !
    !                                                                      !
    !  R...  DO: SPLINE MOMENTS CALCULATED BY SUBROUTINE 'INITSP'          !
    !                                                                      !
    !----------------------------------------------------------------------!
    !                                                                      !
    !  PROGRAMMER:                                                         !
    !  RENE FORSBERG, JUNE 1983                                            !
    !                                                                      !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IMPLICIT REAL*8(A-H, P-Z)
          DIMENSION Y(N), R(N)
          IF (X.LT.1) THEN
            SPLINE = Y(1) + (X-1)*(Y(2)-Y(1)-R(2)/6)
          ELSEIF (X.GT.N) THEN
            SPLINE = Y(N) + (X-N)*(Y(N)-Y(N-1)+R(N-1)/6)
          ELSE
            J = IFRAC(X)
            XX = X - J
            SPLINE = Y(J) + &
         &           XX * ((Y(J+1)-Y(J)-R(J)/3-R(J+1)/6) + &
         &           XX * (R(J)/2 + &
         &           XX * (R(J+1)-R(J))/6))
          ENDIF
          RETURN
          END
    !
    !-----------------------------------------------------------------------
    !
          FUNCTION IFRAC(R)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                      !
    !                    F U N ! T I O N   I F R A !                       !
    !                  ===============================                     !
    !                                                                      !
    !  SUBROUTINE GIVING TRUE INTEGER PART OF A REAL E.G.                  !
    !                                                                      !
    !    FOR   1. = R < 2.   IS    IFRAC = 1                               !
    !    FOR   0. = R < 1.   IS    IFRAC = 0                               !
    !    FOR  -1. = R < 0.   IS    IFRAC =-1                               !
    !    FOR  -2. = R <-1.   IS    IFRAC =-2                               !
    !                                                                      !
    !  RF, JUNE 1983                                                       !
    !  HD, JANUARY 1987                                                    !
    !                                                                      !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IMPLICIT REAL*8(A-H,O-Z)
          IFRAC=R
          IF (R.GE.0) RETURN
          IF (R.EQ.IFRAC) RETURN
          IFRAC = IFRAC - 1
          RETURN
          END
    !
    !-----------------------------------------------------------------------
    !