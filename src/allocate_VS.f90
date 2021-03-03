program SET_MAP
    !==========================================
    ! convert HydroWeb, HydroSat, IICESat VS to CaMa-Flood grids
    ! Menaka@IIS
    ! 2021.01.22
    !==========================================
    implicit none
    ! index CaMa
    integer                       ::  iXX, iYY !, jXX, jYY, kXX, kYY, uXX, uYY
    integer                       ::  nXX, nYY, nFL                   !! x-y matrix GLOBAL
    real                          ::  gsize, csize                    !! grid size [degree]
    real                          ::  west, north, east, south
    character*128                 ::  buf
    character*128                 ::  camadir, map, tag               !! CaMa dir , map[e.g.:glb_15min], tag[e.g. 1min, 15sec, 3sec]
    ! integer                       ::  ngrid
    integer                       ::  ios
    integer                       ::  cnum, num, mwin
    character*128,allocatable     ::  datanames(:)
    character*128                 ::  dataname

    ! index higer resoultion [1min, 15sec, 3sec]
    integer                       ::  ix, iy, jx, jy, kx, ky, dx, dy
    integer                       ::  nx, ny
    integer                       ::  hres
    real                          ::  west1, south1, north1, east1
    character*7                   ::  cname
    
    ! coarse resolution 
    ! integer                       ::  dXX, dYY
    ! integer                       ::  ix0, iy0 !! allocated station xy
    
    ! input
    real,allocatable              ::  uparea(:,:),elevtn(:,:),nxtdst(:,:) !! drainage area (GRID base),elevation
    integer,allocatable           ::  basin(:,:), biftag(:,:)             !! next grid X
    integer,allocatable           ::  nextXX(:,:), nextYY(:,:)
    ! integer,allocatable           ::  upstXX(:,:,:), upstYY(:,:,:)
    ! real,allocatable              ::  glon(:), glat(:)
    
    ! higher resoultion data
    integer*2,allocatable         ::  dwx1m(:,:), dwy1m(:,:)
    integer*2,allocatable         ::  catmXX(:,:), catmYY(:,:)
    integer*1,allocatable         ::  catmZZ(:,:), visual(:,:)
    real,allocatable              ::  upa1m(:,:), ele1m(:,:)
    real,allocatable              ::  riv1m(:,:), flddif(:,:), hand(:,:)
    ! real,allocatable              ::  lon1m(:), lat1m(:)
    
    ! calculation
    integer                       ::  nn
    real                          ::  lag, lag_now!, upa
    ! real                          ::  err1, slope, threshold
    
    ! Station list
    ! character*128         ::  id
    !integer              ::  id
    character*128                 ::  id, station,river,bsn,country
    character*128                 ::  sat,sday,eday,stime,etime,status
    real                          ::  lat0, lon0, ele0 !, lat, lon,area
    real                          ::  egm08,egm96
    integer                       ::  flag
    
    ! file
    character*128                 ::  hiresmap,regmap, finp
    character*128                 ::  outdir
    character*128                 ::  rfile1, rlist
    character*128                 ::  wfile1
    ! ===============================================
    call getarg(1,buf)
    read(buf,*) west1
    call getarg(2,buf)
    read(buf,*) south1
    call getarg(3,buf)
     read(buf,"(A)") dataname
    call getarg(4,buf)
     read(buf,"(A)") camadir
    call getarg(5,buf)
     read(buf,"(A)") map
    call getarg(6,buf)
     read(buf,"(A)") tag
    call getarg(7,buf)
     read(buf,"(A)") outdir
    
    !print*, west1, south1, trim(dataname)
    !==
    finp=trim(camadir)//"/map/"//trim(map)//"/params.txt"
    open(11,file=finp,form='formatted')
    read(11,*) nXX
    read(11,*) nYY
    read(11,*) nFL
    read(11,*) gsize
    read(11,*) west
    read(11,*) east
    read(11,*) south
    read(11,*) north
    close(11)
    !=======
    ! allocate(datanames(1))
    ! datanames=[character(len=128) :: "GRRATS" ] !"HydroWeb","CGLS","HydroSat","ICESat"]
    !=======
    ! tag
    hres=60
    cnum=240
    csize=1./dble(cnum)
    mwin=30
    if (trim(tag)=="1min") then
        hres=60
        cnum=60
        csize=1./dble(cnum)
        mwin=1
    elseif( trim(tag)=="15sec" )then
        hres=4*60
        cnum=240
        csize=1./dble(cnum)
        mwin=30
    elseif( trim(tag)=="3sec" )then
        hres=20*60
        cnum=1200
        csize=1./dble(cnum)
        mwin=10
    elseif( trim(tag)=="1sec" )then
        hres=60*60
        cnum=3600
        csize=1./dble(cnum)
        mwin=1
    end if
    !
    ! ==========
    regmap=trim(camadir)//"/map/"//trim(map)
    
    ! print *, regmap
    allocate(uparea(nXX,nYY),basin(nXX,nYY),elevtn(nXX,nYY),nxtdst(nXX,nYY))
    allocate(nextXX(nXX,nYY),nextYY(nXX,nYY),biftag(nXX,nYY))

    !================
    ! high resolution
    nx=int( mwin*cnum )
    ny=int( mwin*cnum )
    allocate(upa1m(nx,ny),catmXX(nx,ny),catmYY(nx,ny),catmZZ(nx,ny),dwx1m(nx,ny),dwy1m(nx,ny))
    allocate(flddif(nx,ny),hand(nx,ny),ele1m(nx,ny),riv1m(nx,ny),visual(nx,ny))
    
    rfile1=trim(regmap)//'/uparea.bin'
    !print *, "uparea", rfile1
    open(11, file=rfile1, form='unformatted', access='direct' , action='READ', recl=4*nXX*nYY,status='old',iostat=ios)
    read(11,rec=1) uparea
    close(11)
    
    !print *, "basin"
    rfile1=trim(regmap)//'/basin.bin'
    open(11, file=rfile1, form='unformatted', access='direct' , action='READ', recl=4*nXX*nYY,status='old',iostat=ios)
    read(11,rec=1) basin
    close(11)
    
    !print *, "nextxy"
    rfile1=trim(regmap)//'/nextxy.bin'
    open(11, file=rfile1, form='unformatted', access='direct' , action='READ', recl=4*nXX*nYY,status='old',iostat=ios)
    read(11,rec=1) nextXX
    read(11,rec=2) nextYY
    close(11)
    
    rfile1=trim(regmap)//'/elevtn.bin'
    open(11, file=rfile1, form='unformatted', access='direct' , action='READ', recl=4*nXX*nYY,status='old',iostat=ios)
    read(11,rec=1) elevtn
    close(11)
    
    rfile1=trim(regmap)//'/nxtdst.bin'
    open(11, file=rfile1, form='unformatted', access='direct' , action='READ', recl=4*nXX*nYY,status='old',iostat=ios)
    read(11,rec=1) nxtdst
    close(11)
    
    rfile1=trim(regmap)//'/biftag.bin'
    open(11, file=rfile1, form='unformatted', access='direct' , action='READ', recl=4*nXX*nYY,status='old',iostat=ios)
    read(11,rec=1) biftag
    close(11)

    !====================
    ! ! write to file
    ! ! need to change station character length a40 -> a50 [CGLS] 2021.1.27 {not yet changed}
    ! wfile1=trim(outdir)//"/altimetry_GRRATS_"//trim(map)//".txt"
    ! open(27, file=wfile1, form='formatted')
    ! write(27,'(6x,a2,38x,a7,7x,a8,9x,a3,7x,a3,7x,a2,6x,a2,5x,a8,2a10,8x,a9)') &
    !     &"ID","station","dataname","lon",&
    !     &"lat","ix","iy","ele_diff",&
    !     &"EGM08","EGM96","satellite"
    !=======
    ! do num=1,1
    !     dataname=datanames(num)
    north1=south1+10.0
    east1=west1+10.0
    !====================================
    if (trim(tag) == "1min") then
        cname=trim(tag)
    else 
        call set_name(west1,south1,cname)
    end if
    !===========
    
    !print *, trim(station),lon0,lat0,egm08,egm96

    hiresmap=trim(camadir)//"/map/"//trim(map)//"/"//trim(tag)//"/"

    rfile1=trim(hiresmap)//trim(cname)//'.catmxy.bin'
    !print *, rfile1
    open(21,file=rfile1,form='unformatted',access='direct',recl=2*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    !print *, "=====",rfile1
    read(21,rec=1) catmXX
    read(21,rec=2) catmYY
    close(21)
    !nloc=nloc+1
    else
        print*, "no file, ", rfile1, ".  Please check the soft link"
        print*, "ln -sf /org_dir/ ", trim(tag)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.catmzz.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) catmZZ
    close(21)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.flddif.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) flddif
    close(21)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.hand.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) hand
    close(21)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.elevtn.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) ele1m
    close(21)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.uparea.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) upa1m
    close(21)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.rivwth.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) riv1m 
    close(21)
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.visual.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) visual
    close(21)
    endif
    ! ===============================================
    ! read data 
    ! ===============================================
    ! print *, 'read '//trim(dataname)//' data'
    rlist='./inp/'//trim(dataname)//'Station_list.txt'
    !print *, rlist
    open(11, file=rlist, form='formatted')
    read(11,*)
    !----
1000 continue
    read(11,*,end=1090) id, station, river, bsn, country, lon0, lat0, ele0, egm08, egm96, sat, stime, etime, status
    ! print*, trim(id), trim(station), trim(river)
    !== regional maps
    if (lon0 < west1 .or. lon0 > east1 .or. lat0 < south1 .or. lat0 > north1) then
        goto 1000
    end if
    !print*, trim(id), " ", trim(station), " ", trim(river), lon0, lat0
    !-----------
    ix=int( (lon0 - west1 )*dble(hres) )+1
    iy=int( (north1 - lat0)*dble(hres) )+1
    !-----------
    ! flag identity
    ! 1 = location was directly found
    ! 2 = found the nearest permenat water
    ! 3 = correction for ocean grids
    ! 4 = bifurication location
    flag=-9
    ! print*, trim(station), "ix:", ix, "iy:",iy
    ! get the nearest west and south 
    ! call westsouth(lon0,lat0,mwin,west1,south1)
    !================
    ! iXX=catmXX(ix,iy)
    ! iYY=catmYY(ix,iy)
    ! print*, iXX, iYY
    !================
    !! permanat water
    kx=ix
    ky=iy
    if( riv1m(ix,iy)/=-9999 .and. riv1m(ix,iy)/=0 )then
        kx=ix
        ky=iy
        flag=1
    else
        nn=5
        lag=1.0e20
        do dy=-nn,nn
            do dx=-nn,nn
            jx=ix+dx
            jy=iy+dy
            !! for regional maps
            if ((east-west)==360.0) then
                if( jx<=0 ) jx=jx+nx
                if( jx>nx ) jx=jx-nx
            else
                if( jx<=0 ) jx=1
                if( jx>nx ) jx=nx
            end if
            !-----
            if ((catmXX(jx,jy)<=0) .or. (catmYY(jx,jy)<=0)) cycle
            lag_now=sqrt((real(dx)**2)+(real(dy)**2))
            if (lag_now<lag) then
                if( riv1m(jx,iy)/=-9999 .and. riv1m(jx,iy)/=0 )then
                ! iXX=catmXX(jx,jy)
                ! iyy=catmYY(jx,jy)
                kx=jx
                ky=jy
                lag=lag_now
                end if
            end if
            end do
        end do
        flag=2
    end if
    !===========
    iXX=catmXX(kx,ky)
    iyy=catmYY(kx,ky)
    !! correction for ocean grids
    if (iXX<=0 .or. iYY<=0) then
        nn=5
        lag=1.0e20
        do dy=-nn,nn
            do dx=-nn,nn
                jx=kx+dx
                jy=ky+dy
                !! for regional maps
                if ((east-west)==360.0) then
                if( jx<=0 ) jx=kx+nx
                if( jx>nx ) jx=kx-nx
                else
                if( jx<=0 ) jx=1
                if( jx>nx ) jx=nx
                end if
                !-----
                if ((catmXX(jx,jy)<=0) .or. (catmYY(jx,jy)<=0)) cycle
                lag_now=sqrt((real(dx)**2)+(real(dy)**2))
                if (lag_now<lag) then
                kx=jx
                ky=jy
                lag=lag_now
                end if
            end do
        end do
        flag=3
    end if
    !===========
    iXX=catmXX(kx,ky)
    iyy=catmYY(kx,ky)
    !============
    ! find maximum uparea perpendicular to river
    ! in case of braided river
    ! considering the bifurcation tag
    if (biftag(iXX,iYY) == 1) then
        ix=iXX
        iy=iYY
        call loc_pepnd(ix,iy,nXX,nYY,nextXX,nextYY,uparea,iXX,iYY)
        flag=4
    end if
    ! print*, flag
    if (iXX > 0 .or. iYY > 0) then
        print '(a30,2x,a60,2x,a10,2x,2f10.2,2x,2i8.0,2x,3f10.2,2x,a15,2x,i4.0)', trim(adjustl(id)),&
        &trim(station), trim(dataname), lon0, lat0, iXX, iYY, elevtn(iXX,iYY)-ele1m(kx,ky),&
        &egm08, egm96, trim(sat), flag
    ! else
    !     print*, "no data"
    end if
        ! write(27,'(a14,2x,a40,2x,a10,2x,2f10.2,2x,2i8.0,2x,3f10.2,2x,a15)') trim(adjustl(id)),& 
        ! &trim(station),trim(dataname), lon0, lat0,iXX, iYY,elevtn(iXX,iYY)-ele1m(kx,ky),&
        ! &egm08, egm96, trim(sat)

    goto 1000
1090 continue
    close(11)
    !---
    deallocate(uparea,basin,elevtn,nxtdst,nextXX,nextYY)
    deallocate(upa1m,catmXX,catmYY,catmZZ,dwx1m,dwy1m,flddif,hand,ele1m,riv1m,visual)
    !=================
    end program SET_MAP
    !***************************************************
    function roundx(ix, nx)
    implicit none
    !-- for input -----------
    integer                     ix, nx
    !-- for output ----------
    integer                     roundx
    !------------------------
    if (ix .ge. 1) then
        roundx = ix - int((ix -1)/nx)*nx
    else
        roundx = nx - abs(mod(ix,nx))
    end if 
    return
    end function roundx
    !*****************************************************************
    subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
    implicit none
    !- for input -----------------
    integer                   ix, iy, nx, ny
    !- for output ----------------
    integer                   iix, iiy,roundx
    !-----------------------------
    if (iy .lt. 1) then
        iiy = 2 - iy
        iix = ix + int(nx/2.0)
        iix = roundx(iix, nx)
    else if (iy .gt. ny) then
        iiy = 2*ny -iy
        iix = ix + int(nx/2.0)
        iix = roundx(iix, nx)
    else
        iiy = iy
        iix = roundx(ix, nx)
    end if
    return
    end subroutine ixy2iixy  
    !*****************************************************************
    subroutine upstream(i,j,nx,ny,nextX,nextY,uparea,x,y)
    implicit none 
    ! find the upstream pixel with closest upstream area to the i,j
    !--
    integer                             :: i,j,nx,ny,x,y
    integer,dimension(nx,ny)            :: nextX,nextY !,rivseq
    real,dimension(nx,ny)               :: uparea !nextdst, 
    !--
    real                                :: dA ! area differnce
    integer                             :: ix,iy,tx,ty, d !,iix,iiy,ud
    ! real                                :: length,rl
    !--
    x=-9999
    y=-9999
    d=10 ! look at 10*10 box
    dA=1.0e20 ! area differnce
    !--
    !write(*,*)i,j
    do tx=i-d,i+d
        do ty=j-d,j+d
        !write(*,*)tx,ty
        call ixy2iixy(tx,ty, nx, ny, ix, iy)
        !write(*,*)nextX(ix,iy),nextY(ix,iy),ix,iy,uparea(ix,iy),rivseq(ix,iy)
        if (nextX(ix,iy) == i .and. nextY(ix,iy) == j) then
            !write(*,*)ix,iy
            if (uparea(ix,iy) - uparea(i,j) < dA) then
                dA=uparea(ix,iy)-uparea(i,j)
                x=ix
                y=iy
                !write(*,*)x,y
            end if
        end if
        end do
    end do
    return
    end subroutine upstream
    !*****************************************************************
    subroutine set_name(lon,lat,cname)
    ! ===============================================
    implicit none
    !
    real            ::  lon, lat

    character*1     ::  ew, sn
    character*2     ::  clat
    character*3     ::  clon
    character*7     ::  cname
    ! ===============================================
    if( lon<0 )then
        ew='w'
        write(clon,'(i3.3)') int(-lon)
    else
        ew='e'
        write(clon,'(i3.3)')  int(lon)
    endif

    if( lat<0 )then
        sn='s'
        write(clat,'(i2.2)') int(-lat)
    else
        sn='n'
        write(clat,'(i2.2)')  int(lat)
    endif

    cname=sn//clat//ew//clon

    end subroutine set_name
    !*****************************************************************
    subroutine westsouth(lon,lat,mwin,west,south)
    implicit none
    !
    real           :: lon, lat
    integer        :: mwin
    real           :: west, south
    !
    real           :: lon1, lat1
    ! get the nearest west and south 
    if (lon>0.0) then
    !print*, lon
    west=real(int(lon - mod(lon,real(mwin))))
    else
    !print*, lon
    lon1=abs(lon)
    west=-(real(int(lon1 - mod(lon1,real(mwin)))) + 10.0)
    end if
    !----------------
    if (lat>0.0) then
    south=real(int(lat - mod(lat,real(mwin))))
    else
    lat1=abs(lat)
    south=-(real(int(lat1 - mod(lat1,real(mwin)))) + 10.0)
    end if
    return
    end subroutine westsouth
    !*****************************************************************
    subroutine loc_pepnd(ix,iy,nx,ny,nextX,nextY,uparea,oxx,oyy)
    ! river location perpendicular to the flowing direction
    implicit none
    integer                      :: ix, iy, nx, ny
    integer,dimension(nx,ny)     :: nextX, nextY
    real,dimension(nx,ny)        :: uparea
    integer                      :: oxx, oyy
    integer,dimension(2)         :: xlist, ylist
    integer                      :: k

    integer                      :: iix, iiy, dx, dy, jx, jy, i
    real                         :: tval 
    integer                      :: dval, D8 ! d8 numbering
    real                         :: upa, upn
    !==============================================
    jx=nextX(ix,iy)
    jy=nextY(ix,iy)
    !--------------
    dx=jx-ix 
    dy=jy-iy
    dval=D8(dx,dy)
    k=0
    if (dval==1 .or. dval==5) then
        k=2
        !-------------------------
        iix=ix
        iiy=iy-1
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(1)=iix 
        ylist(1)=iiy
        !-------------------------
        iix=ix
        iiy=iy+1
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(2)=iix 
        ylist(2)=iiy
    elseif (dval==3 .or. dval==7) then
        k=2
        !-------------------------
        iix=ix-1
        iiy=iy
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(1)=iix 
        ylist(1)=iiy
        !-------------------------
        iix=ix+1
        iiy=iy
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(2)=iix 
        ylist(2)=iiy
    elseif (dval==4 .or. dval==8) then
        k=2
        !-------------------------
        iix=ix+1
        iiy=iy-1
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(1)=iix 
        ylist(1)=iiy
        !-------------------------
        iix=ix-1
        iiy=iy+1
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(2)=iix 
        ylist(2)=iiy
    elseif (dval==2 .or. dval==6) then
        k=2
        !-------------------------
        iix=ix-1
        iiy=iy-1
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(1)=iix 
        ylist(1)=iiy
        !-------------------------
        iix=ix+1
        iiy=iy+1
        call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        xlist(2)=iix 
        ylist(2)=iiy
    end if
    !-----------------------------
    upa=uparea(ix,iy)
    upn=uparea(ix,iy)
    oxx=ix
    oyy=iy
    do i=1,k 
        iix=xlist(i)
        iiy=ylist(i)
        upn=uparea(iix,iiy)
        if (upa < upn) then
            oxx=iix 
            oyy=iiy 
            upa=upn 
        end if
    end do 
    return
    end subroutine loc_pepnd
    !*****************************************************************
    function D8(dx,dy)
    implicit none
    integer                         :: dx, dy
    integer                         :: D8, dval
    real                            :: tval
    !-------------------------------
    ! Angle (degree) |  tan value |
    !-------------------------------
    !      22.5      |  0.4142135 |
    !      67.5      |  2.4142135 |
    !-------------------------------
    ! D 8 graph
    !-----------|
    ! 4 | 3 | 2 |
    !----------|
    ! 5 | 0 | 1 |
    !-----------|
    ! 6 | 7 | 8 |
    !-----------|
    if (dx == 0)  then
        if (dy > 0) dval=3
        if (dy < 0) dval=7
    elseif (dy == 0)  then
        if (dx > 0) dval=1
        if (dx < 0) dval=5
    elseif (dx > 0 .and. dy > 0) then
        tval=abs(real(dy)/real(dx))
        if (tval > 0 .and. tval <= 0.4142135) dval=1 
        if (tval > 0.4142135 .and. tval <= 2.4142135) dval=2
        if (tval < 2.4142135) dval=3
    elseif (dx > 0 .and. dy < 0) then
        tval=abs(real(dy)/real(dx))
        if (tval > 0 .and. tval <= 0.4142135) dval=1 
        if (tval > 0.4142135 .and. tval <= 2.4142135) dval=8
        if (tval < 2.4142135) dval=7
    elseif (dx < 0 .and. dy < 0) then
        tval=abs(real(dy)/real(dx))
        if (tval > 0 .and. tval <= 0.4142135) dval=5 
        if (tval > 0.4142135 .and. tval <= 2.4142135) dval=6
        if (tval < 2.4142135) dval=7
    elseif (dx < 0 .and. dy < 0) then
        tval=abs(real(dy)/real(dx))
        if (tval > 0 .and. tval <= 0.4142135) dval=5 
        if (tval > 0.4142135 .and. tval <= 2.4142135) dval=4
        if (tval < 2.4142135) dval=3
    end if
    D8=dval
    return
    end function D8
    !*****************************************************************