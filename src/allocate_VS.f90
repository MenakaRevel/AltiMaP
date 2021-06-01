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
    integer                       ::  ibx, iby ! bifurication locations
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
    integer*1,allocatable         ::  catmZZ(:,:), visual(:,:), flwdir(:,:)
    real,allocatable              ::  upa1m(:,:), ele1m(:,:)
    real,allocatable              ::  riv1m(:,:), flddif(:,:), hand(:,:)
    ! real,allocatable              ::  lon1m(:), lat1m(:)
    
    ! calculation
    integer                       ::  nn
    real                          ::  lag, lag_now!, upa
    real                          ::  lat1, lon1, lat2, lon2
    real                          ::  down_dist, flow_dist, hubeny_real
    ! real                          ::  err1, slope, threshold
    
    ! Station list
    ! character*128         ::  id
    !integer              ::  id
    character*128                 ::  id, station,river,bsn,country
    character*128                 ::  sat,sday,eday,stime,etime,status
    real                          ::  lat0, lon0, ele0 !, lat, lon,area
    real                          ::  egm08,egm96, diffdist
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
    allocate(flddif(nx,ny),hand(nx,ny),ele1m(nx,ny),riv1m(nx,ny),visual(nx,ny),flwdir(nx,ny))
    
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

    rfile1=trim(hiresmap)//trim(cname)//'.flwdir.bin'
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) flwdir
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
    !-------------------------------------
    ix=int( (lon0 - west1 )*dble(hres) )+1
    iy=int( (north1 - lat0)*dble(hres) )+1
    !-----------
    ! flag identity
    ! 1 = location was directly found
    ! 2 = location was on the unit-catchment outlet
    ! 3 = found the nearest permeant water
    ! 4 = correction for ocean grids
    ! 5 = bifurcation location
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
    lat1=lat0
    lon1=lon0
    ! print*, "==========================================================="
    ! print*, "Initial allocation: ",kx, ky, visual(kx,ky)
    ! if( riv1m(ix,iy)/=-9999 .and. riv1m(ix,iy)/=0 )then
    if (visual(kx,ky) == 10) then  !! river channel
        ! print*, "flag: 1  ","river channel"
        kx=ix
        ky=iy
        flag=1
    else if (visual(kx,ky) == 20) then  !! unit-catchment mouth
        ! print*, "flag: 2  ","unit-catchment mouth"
        kx=ix
        ky=iy
        flag=2
    else if (visual(kx,ky) == 2 .or. visual(kx,ky) == 3 .or. visual(kx,ky) == 5 .or. visual(kx,ky) == 7) then   !! correction for land/grid box, boundry to channel
        ! print*, "flag: 3  ","correction for land to channel"
        nn=10
        lag=1.0e20
        do dy=-nn,nn
            do dx=-nn,nn
                jx=ix+dx
                jy=iy+dy
                if ( jx<=0 ) cycle !jx=1
                if ( jx>nx ) cycle !jx=nx
                if ( jy<=0 ) cycle !jy=1
                if ( jy>nx ) cycle !jy=ny
                ! ! !! for regional maps
                ! ! if ( (east-west)==360.0 ) then
                ! !     if ( jx<=0 ) jx=jx+nx
                ! !     if ( jx>nx ) jx=jx-nx
                ! ! else
                ! !     if ( jx<=0 ) jx=1
                ! !     if ( jx>nx ) jx=nx
                ! ! end if
                !-----
                if ( kx == jx .and. ky == jy) cycle
                if ( (catmXX(jx,jy)<=0) .or. (catmYY(jx,jy)<=0) ) cycle
                if ( visual(jx,jy) /= 10 ) cycle
                ! lag_now=sqrt((real(dx)**2)+(real(dy)**2))
                ! print*, "===",kx,ky,jx,jy,"==="
                ! lag_now=flow_dist(kx,ky,jx,jy,west1,south1,csize,flwdir,visual,nx,ny,hiresmap)
                lat2=north1 - (jy-1)*(1/dble(hres))
                lon2=west1 + (jx-1)*(1/dble(hres))
                lag_now=hubeny_real(lat1, lon1, lat2, lon2)
                if ( lag_now == -9999.0 ) cycle
                ! print*, lag, lag_now
                if ( lag_now < lag ) then
                    ! if( riv1m(jx,jy)/=-9999 .and. riv1m(jx,jy)/=0 )then
                    if ( visual(jx,jy) == 10 ) then
                    ! iXX=catmXX(jx,jy)
                    ! iyy=catmYY(jx,jy)
                        kx=jx
                        ky=jy
                        lag=lag_now
                        ! print*, "Found new location: ",flag, kx,ky,lag
                    end if
                end if
            end do
        end do
        if ( kx /= ix .or. ky /= iy ) then
            flag=3
        else
            flag=1
        end if
    else if ( visual(kx,ky)==0 .or. visual(kx,ky)==1 .or. visual(kx,ky)==25 ) then !! correction for ocean grids
        ! print*, "flag: 4  ","correction for ocean grids"
        nn=10
        lag=1.0e20
        do dy=-nn,nn
            do dx=-nn,nn
                jx=ix+dx
                jy=iy+dy
                if ( jx<=0 ) cycle !jx=1
                if ( jx>nx ) cycle !jx=nx
                if ( jy<=0 ) cycle !jy=1
                if ( jy>nx ) cycle !jy=ny
                ! ! !! for regional maps
                ! ! if ( (east-west)==360.0 ) then
                ! !     if( jx<=0 ) jx=kx+nx
                ! !     if( jx>nx ) jx=kx-nx
                ! ! else
                ! !     if( jx<=0 ) jx=1
                ! !     if( jx>nx ) jx=nx
                ! ! end if
                !-----
                if ( kx == jx .and. ky == jy) cycle
                if ((catmXX(jx,jy) <= 0) .or. (catmYY(jx,jy) <= 0)) cycle
                if ( visual(jx,jy) /= 10) cycle
                ! lag_now=sqrt((real(dx)**2)+(real(dy)**2))
                ! print*, "===",kx,ky,jx,jy,"==="
                ! lag_now=flow_dist(kx,ky,jx,jy,west1,south1,csize,flwdir,visual,nx,ny,hiresmap)
                lat2=north1 - (jy-1)*(1/dble(hres))
                lon2=west1 + (jx-1)*(1/dble(hres))
                lag_now=hubeny_real(lat1, lon1, lat2, lon2)
                if ( lag_now == -9999.0 ) cycle
                ! print*, lag, lag_now
                if ( lag_now < lag ) then
                    if ( visual(jx,jy) == 10 ) then
                        kx=jx
                        ky=jy
                        lag=lag_now
                        ! print*, "Found new location: ",flag, kx,ky,lag
                    end if
                end if
            end do
        end do
        if ( kx /= ix .or. ky /= iy ) then
            flag=4
        else
            flag=1
        end if
    else
        flag=9
        ! print*, "flag: 9 ->", kx, ky, visual(kx,ky)
    end if
    !
    !---
    if ( kx < 1 .or. ky < 1 .or. kx > nx .or. ky > ny ) then 
        goto 1000
    end if
    ! print*, kx, ky
    !===========
    iXX=catmXX(kx,ky)
    iyy=catmYY(kx,ky)
    !!============
    if ( iXX < 1 .or. iYY < 1 .or. iXX > nXX .or. iYY > nYY ) then 
        goto 1000
    end if
    !============
    ! find maximum uparea perpendicular to river
    ! in case of braided river
    ! considering the bifurcation tag
    if ( biftag(iXX,iYY) == 1 ) then
        ! call loc_pepnd(ix,iy,nXX,nYY,nextXX,nextYY,uparea,iXX,iYY)
        call loc_pepndD8(kx,ky,nx,ny,flwdir,visual,uparea,west1,south1,hiresmap,ibx,iby)
        if ( ibx/=-9 .and. iby/=-9 ) then
            ! print*, "burification location", ibx, iby
            flag=5
            kx=ibx
            ky=iby
        end if
    end if
    ! print*, flag
    ! print*, "After allocation:   ",kx, ky, visual(kx,ky), flag
    !---
    ! print*, kx, ky
    if ( kx < 1 .or. ky < 1 .or. kx > nx .or. ky > ny ) then 
        goto 1000
    end if
    !--
    if (visual(kx,ky)==10) then
        diffdist=down_dist(kx,ky,west1,south1,csize,flwdir,visual,nx,ny,hiresmap)
    else
        diffdist=0.0
    end if
    !===========
    iXX=catmXX(kx,ky)
    iyy=catmYY(kx,ky)
    !============
    ! print*, trim(station), flag, diffdist*1e-3
    if (iXX > 0 .or. iYY > 0) then
        print '(a30,2x,a65,2x,a10,2x,2f10.2,2x,2i8.0,2x,3f10.2,2x,a15,2x,f13.2,2x,i4.0,2x,2i8.0)', trim(adjustl(id)),&
        &trim(station), trim(dataname), lon0, lat0, iXX, iYY, elevtn(iXX,iYY)-ele1m(kx,ky),&
        &egm08, egm96, trim(sat), diffdist*1e-3, flag, kx, ky
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
    deallocate(upa1m,catmXX,catmYY,catmZZ,dwx1m,dwy1m,flddif,hand,ele1m,riv1m,visual,flwdir)
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

    integer                      :: iix, iiy, dx, dy, jx, jy, i, j
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
    k=5
    if (dval==1 .or. dval==5) then
        j=1
        do i=-k,k
            iix=ix 
            iiy=iy+i
            call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
        ! k=2
        ! !-------------------------
        ! iix=ix
        ! iiy=iy-1
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(1)=iix 
        ! ylist(1)=iiy
        ! !-------------------------
        ! iix=ix
        ! iiy=iy+1
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(2)=iix 
        ! ylist(2)=iiy
    elseif (dval==3 .or. dval==7) then
        j=1
        do i=-k,k
            iix=ix+i 
            iiy=iy
            call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
        ! k=2
        ! !-------------------------
        ! iix=ix-1
        ! iiy=iy
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(1)=iix 
        ! ylist(1)=iiy
        ! !-------------------------
        ! iix=ix+1
        ! iiy=iy
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(2)=iix 
        ! ylist(2)=iiy
    elseif (dval==4 .or. dval==8) then
        j=1
        do i=-k,k
            iix=ix+i 
            iiy=iy-i
            call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
        ! k=2
        ! !-------------------------
        ! iix=ix+1
        ! iiy=iy-1
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(1)=iix 
        ! ylist(1)=iiy
        ! !-------------------------
        ! iix=ix-1
        ! iiy=iy+1
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(2)=iix 
        ! ylist(2)=iiy
    elseif (dval==2 .or. dval==6) then
        j=1
        do i=-k,k
            iix=ix+i 
            iiy=iy+i
            call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
        ! k=2
        ! !-------------------------
        ! iix=ix-1
        ! iiy=iy-1
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(1)=iix 
        ! ylist(1)=iiy
        ! !-------------------------
        ! iix=ix+1
        ! iiy=iy+1
        ! call ixy2iixy(iix,iiy,nx,ny,iix,iiy)
        ! xlist(2)=iix 
        ! ylist(2)=iiy
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
    function down_dist(ix,iy,west,south,csize,flwdir,visual,nx,ny,hiresmap)
    implicit none
    ! calculate the distance from the exact VS location to unit catchment mouth
    real                            :: down_dist
    real                            :: west, south, csize
    integer                         :: nx, ny
    integer*1,dimension(nx,ny)      :: flwdir, visual
    integer                         :: ix, iy, dval, dx,dy
    character*128                   :: hiresmap

    integer                         :: iix, iiy, iix0, iiy0
    real                            :: west0, south0, north0, tval
    real                            :: lon1, lat1, lon2, lat2
    integer*1,dimension(nx,ny)      :: flwdir0, visual0
    real                            :: hubeny_real
    !--------------------
    ! visual
    ! 0  - sea
    ! 1  - land(undefied)
    ! 2  - land(defined in CaMa)
    ! 3  - grid box
    ! 5  - catchment boundry
    ! 10 - channel
    ! 20 - outlet pixel
    ! 25 - river mouth
    !--------------------
    flwdir0=flwdir
    visual0=visual
    west0=west+csize/2.0
    south0=south+csize/2.0
    north0=south+10.0+csize/2.0
    down_dist = 0.0
    iix = ix 
    iiy = iy
    lon1=west0+real(ix)*csize
    lat1=north0-real(iy)*csize
    do while (visual0(iix,iiy) /= 20)
        if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
            ! call got_to_next_tile(iix,iiy,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
            ! west0=west+csize/2.0
            ! south0=south+csize/2.0
            ! north0=south+10.0+csize/2.0
            ! print*, "go to next tile", west0,south0
            ! flag=-9
            exit
        end if
        ! river mouth
        if (flwdir0(iix,iiy) == -9 ) then
            ! print*, "River mouth", visual(iix,iiy)
            exit
        end if
        ! land
        if (visual0(iix,iiy) == 2 ) then
            ! print*, "Not in the river channel", visual(iix,iiy)
            exit
        end if
        ! outlet pixel
        if (visual0(iix,iiy) == 20 ) then
            ! print*, "Outlet pixel", visual(iix,iiy)
            exit
        end if
        dval=flwdir0(iix,iiy)
        call next_D8(dval,dx,dy)
        iix = iix + dx 
        iiy = iiy + dy
        if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
            iix0 = iix
            iiy0 = iiy 
            call got_to_next_tile(iix0,iiy0,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
            !replace west0, south0, flwdir0, visual0
            west0=west0+csize/2.0
            ! south0=south0+csize/2.0
            north0=south0+10.0+csize/2.0
            ! print*, west, south, "to ",west0, south0
            ! print*, iix0, iiy0, "to ", iix, iiy
            ! exit
        end if
        lon2=lon1+real(dx)*csize 
        lat2=lat1+real(dy)*csize
        down_dist=down_dist+hubeny_real(lat1, lon1, lat2, lon2)
        ! print*, iix, iiy, lon1, lat1, down_dist ,visual(iix,iiy), dval
        lon1=lon2
        lat1=lat2
    end do 
    return
    end function down_dist
    !*****************************************************************
    subroutine next_D8(dval,dx,dy)
    implicit none
    integer                         :: dval
    integer                         :: dx, dy
    real                            :: tval
    ! -----------|
    !  D 8 graph
    !|-----------|
    !| 8 | 1 | 2 |
    !|-----------|
    !| 7 | 0 | 3 |
    !|-----------|
    !| 6 | 5 | 4 |
    !|-----------|
    if (dval == 1) then
        dx = 0
        dy = -1
    end if
    if (dval == 2) then
        dx = 1
        dy = -1
    end if
    if (dval == 3) then
        dx = 1
        dy = 0
    end if
    if (dval == 4) then
        dx = 1
        dy = 1
    end if
    if (dval == 5) then
        dx = 0
        dy = 1
    end if
    if (dval == 6) then
        dx = -1
        dy = 1
    end if
    if (dval == 7) then
        dx = -1
        dy = 0
    end if
    if (dval == 8) then
        dx = -1
        dy = -1
    end if
    return
    end subroutine next_D8
    !*****************************************************************
    function hubeny_real(lat1, lon1, lat2, lon2)
    implicit none
    !-- for input -----------
    real                                  lat1, lon1, lat2, lon2
    !-- for output-----------
    real                                  hubeny_real  ! (m)

    !-- for calc ------------
    real,parameter                     :: pi = atan(1.0)*4.0
    real,parameter                     :: a  = 6378137
    real,parameter                     :: b  = 6356752.314140
    real,parameter                     :: e2 = 0.00669438002301188
    real,parameter                     :: a_1_e2 = 6335439.32708317
    real                                  M, N, W
    real                                  latrad1, latrad2, lonrad1, lonrad2
    real                                  latave, dlat, dlon
    real                                  dlondeg
    !------------------------
    latrad1   = lat1 * pi / 180.0
    latrad2   = lat2 * pi / 180.0
    lonrad1   = lon1 * pi / 180.0
    lonrad2   = lon2 * pi / 180.0
    !
    latave    = (latrad1 + latrad2)/2.0
    dlat      = latrad2 - latrad1
    dlon      = lonrad2 - lonrad1
    !
    dlondeg   = lon2 - lon1
    if ( abs(dlondeg) .gt. 180.0) then
        dlondeg = 180.0 - mod(abs(dlondeg), 180.0)
        dlon    = dlondeg * pi / 180.0
    end if
    !-------
    W  = sqrt(1.0 - e2 * sin(latave)**2.0 )
    M  =  a_1_e2 / (W**3.0)
    N  =  a / W
    hubeny_real  = sqrt( (dlat * M)**2.0 + (dlon * N * cos(latave))**2.0 )
    return
    end functioN hubeny_real
    !*****************************************************************
    function flow_dist(ix,iy,jx,jy,west,south,csize,flwdir,visual,nx,ny,hiresmap)
    implicit none
    ! calculate the along river distance 
    real                            :: flow_dist
    real                            :: west, south, csize
    integer                         :: nx, ny
    integer*1,dimension(nx,ny)      :: flwdir, visual
    integer                         :: ix, iy, jx, jy, dval, dx, dy
    character*128                   :: hiresmap
    !------------
    integer                         :: iix, iiy, ixx, iyy, iix0, iiy0, flag
    real                            :: west0, south0, north0, tval
    real                            :: lon1, lat1, lon2, lat2
    integer*1,dimension(nx,ny)      :: flwdir0, visual0
    real                            :: hubeny_real
    !--------------------
    flwdir0=flwdir
    visual0=visual
    west0=west+csize/2.0
    south0=south+csize/2.0
    north0=south+10.0+csize/2.0
    !-----------
    flag= 1
    !-----------
    flow_dist = 0.0
    iix = ix 
    iiy = iy
    ixx = jx
    iyy = jy
    lon1=west0+real(ix)*csize
    lat1=north0-real(iy)*csize
    do while ( iix /= ixx .or. iiy /= iyy )
        if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
            ! call got_to_next_tile(iix,iiy,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
            ! west0=west+csize/2.0
            ! south0=south+csize/2.0
            ! north0=south+10.0+csize/2.0
            ! print*, "go to next tile", west0,south0
            flag=-9
            exit
        end if
        if (flwdir0(iix,iiy) == -9 ) then
            ! print*, "River mouth", visual(iix,iiy)
            flag=-1
            exit
        end if
        if ( visual0(iix,iiy) == 20 .or. visual0(iix,iiy) == 25) then
            flag=-1
            exit
        endif
        ! if ( visual(iix,iiy) /= 10 ) then
        !     print*, "not river channel", visual(iix,iiy), flwdir(iix,iiy)
        !     flag=-1
        !     exit
        ! endif
        dval=flwdir0(iix,iiy)
        call next_D8(dval,dx,dy)
        iix = iix + dx 
        iiy = iiy + dy
        if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
            iix0 = iix
            iiy0 = iiy
            ! call got_to_next_tile(iix,iiy,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
            ! west0=west+csize/2.0
            ! south0=south+csize/2.0
            ! north0=south+10.0+csize/2.0
            ! print*, "go to next tile", west0,south0
            flag=-9
            exit
        end if
        lon2=lon1+real(dx)*csize 
        lat2=lat1+real(dy)*csize
        flow_dist=flow_dist+hubeny_real(lat1, lon1, lat2, lon2)
        ! print*, flag, lon1, lat1, flow_dist ,visual(iix,iiy), dval
        lon1=lon2
        lat1=lat2
    end do
    !-------------------
    if (flag == -1) then
        flow_dist = 0.0
        iix = jx 
        iiy = jy
        ixx = ix
        iyy = iy
        lon1=west0+real(ix)*csize
        lat1=north0-real(iy)*csize
        do while ( iix /= jy .or. iiy /= jy )
            if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
                ! call got_to_next_tile(iix0,iiy0,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
                ! west0=west+csize/2.0
                ! south0=south+csize/2.0
                ! north0=south+10.0+csize/2.0
                ! print*, "go to next tile", west0,south0
                flag=-9
                exit
            end if
            if (flwdir0(iix,iiy) == -9 ) then
                ! print*, "River mouth", visual(iix,iiy)
                flag=-9
                exit
            end if
            if ( visual0(iix,iiy) == 20 .or. visual0(iix,iiy) == 25) then
                flag=-9
                exit
            endif
            ! if ( visual(iix,iiy) /= 10 ) then
            !     print*, "not river channel", visual(iix,iiy)
            !     flag=-9
            !     exit
            ! endif
            dval=flwdir0(iix,iiy)
            call next_D8(dval,dx,dy)
            iix = iix + dx 
            iiy = iiy + dy
            if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
                iix0 = iix
                iiy0 = iiy
                ! call got_to_next_tile(iix0,iiy0,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
                ! west0=west+csize/2.0
                ! south0=south+csize/2.0
                ! north0=south+10.0+csize/2.0
                ! print*, "go to next tile", west0,south0
                flag=-9
                exit
            end if
            lon2=lon1+real(dx)*csize 
            lat2=lat1+real(dy)*csize
            flow_dist=flow_dist+hubeny_real(lat1, lon1, lat2, lon2)
            ! print*, flag, lon1, lat1, flow_dist ,visual(iix,iiy), dval
            lon1=lon2
            lat1=lat2
        end do
    end if
    if ( flag == -9 ) flow_dist=-9999.0
    return
    end function flow_dist
    !*****************************************************************
    subroutine next_P2D8(dval,dx,dy)
    implicit none
    integer                         :: dval
    integer                         :: dx, dy
    real                            :: tval
    !-------------------------------------------------------
    ! determine the perpendicular direction to D8 directions
    !-------------------------------------------------------
    ! -----------|
    !  D 8 graph
    !|-----------|
    !| 8 | 1 | 2 |
    !|-----------|
    !| 7 | 0 | 3 |
    !|-----------|
    !| 6 | 5 | 4 |
    !|-----------|
    if (dval == 1) then
        dx = 0
        dy = -1
    end if
    if (dval == 2) then
        dx = 1
        dy = -1
    end if
    if (dval == 3) then
        dx = 1
        dy = 0
    end if
    if (dval == 4) then
        dx = 1
        dy = 1
    end if
    if (dval == 5) then
        dx = 0
        dy = 1
    end if
    if (dval == 6) then
        dx = -1
        dy = 1
    end if
    if (dval == 7) then
        dx = -1
        dy = 0
    end if
    if (dval == 8) then
        dx = -1
        dy = -1
    end if
    return
    end subroutine next_P2D8
    !*****************************************************************
    subroutine loc_pepndD8(ix,iy,nx,ny,flwdir,visual,uparea,west0,south0,hiresmap,oxx,oyy)
    ! river location perpendicular to the flowing direction
    implicit none
    integer                      :: ix, iy, nx, ny
    integer*1,dimension(nx,ny)   :: flwdir, visual
    real,dimension(nx,ny)        :: uparea
    integer                      :: oxx, oyy
    integer,allocatable          :: xlist(:), ylist(:)
    character*128                :: hiresmap
    real                         :: west0, south0
    integer                      :: k

    integer                      :: iix, iiy, dx, dy, j0x, j0y, jx, jy, i, j
    real                         :: tval 
    integer                      :: dval, D8 ! d8 numbering
    real                         :: upa, upn
    !==============================================
    ! -----------|
    !  D 8 graph
    !|-----------|
    !| 8 | 1 | 2 |
    !|-----------|
    !| 7 | 0 | 3 |
    !|-----------|
    !| 6 | 5 | 4 |
    !|-----------|
    !-------------
    dval=flwdir(ix,iy)
    k=30*6
    allocate (xlist(2*k+1),ylist(2*k+1))
    xlist=-9
    ylist=-9
    if (dval==1 .or. dval==5) then
        j=1
        do i=-k,k
            iix=ix+i  
            iiy=iy
            if (iix > nx) cycle
            if (iix < 1) cycle
            if (iiy > ny) cycle
            if (iiy < 1) cycle
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
    elseif (dval==3 .or. dval==7) then
        j=1
        do i=-k,k
            iix=ix 
            iiy=iy+i
            if (iix > nx) cycle
            if (iix < 1) cycle
            if (iiy > ny) cycle
            if (iiy < 1) cycle
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
    elseif (dval==4 .or. dval==8) then
        j=1
        do i=-k,k
            iix=ix+i 
            iiy=iy-i
            if (iix > nx) cycle
            if (iix < 1) cycle
            if (iiy > ny) cycle
            if (iiy < 1) cycle
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
    elseif (dval==2 .or. dval==6) then
        j=1
        do i=-k,k
            iix=ix+i 
            iiy=iy+i
            if (iix > nx) cycle
            if (iix < 1) cycle
            if (iiy > ny) cycle
            if (iiy < 1) cycle
            xlist(j)=iix 
            ylist(j)=iiy
            j=j+1
        end do
    end if
    !-----------------------------
    upa=uparea(ix,iy)
    upn=uparea(ix,iy)
    !---
    oxx=-9
    oyy=-9
    call next_out(ix,iy,flwdir,visual,nx,ny,west0,south0,hiresmap,j0x,j0y)
    if ( j0x/=-9 .or. j0y/=-9 ) then
        do i=1,k 
            iix=xlist(i)
            iiy=ylist(i)
            if (iix == -9 .or. iiy == -9 ) cycle
            if (iix <= 0 .or. iiy <= 0 ) cycle
            ! print*, "L1180", iix, iiy
            call next_out(iix,iiy,flwdir,visual,nx,ny,west0,south0,hiresmap,jx,jy)
            if (jx == -9 .or. jy == -9 ) cycle
            if ( j0x==jx .and. j0y==jy ) then
                upn=uparea(iix,iiy)
                if (upa < upn) then
                    oxx=iix 
                    oyy=iiy 
                    upa=upn 
                end if
            end if
        end do 
    end if
    deallocate (xlist,ylist)
    return
    end subroutine loc_pepndD8
    !*****************************************************************
    subroutine next_out(ix,iy,flwdir,visual,nx,ny,west,south,hiresmap,jx,jy)
    implicit none
    ! calculate the along river distance 
    ! real                            :: flow_dist
    ! real                            :: west, south, csize
    integer                         :: nx, ny
    integer*1,dimension(nx,ny)      :: flwdir, visual
    integer                         :: ix, iy, jx, jy, dval, dx, dy
    character*128                   :: hiresmap
    real                            :: west, south
    !------------
    integer                         :: iix, iiy, ixx, iyy, flag
    integer*1,dimension(nx,ny)      :: flwdir0, visual0
    real                            :: west0, south0
    !
    flwdir0=flwdir
    visual0=visual
    west0=west
    south0=south
    !
    flag=0
    iix=ix
    iiy=iy
    !
    do while (flag<2)
        if ( iix < 1 .or. iiy < 1 .or. iix > nx .or. iiy > ny ) then
            ! call got_to_next_tile(iix,iiy,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
            ! west0=west
            ! south0=south
            ! print*, "go to next tile", west0,south0
            flag=-9
            exit
        end if
        if (flwdir0(iix,iiy) == -9 ) then
            ! print*, "River mouth", visual(iix,iiy)
            flag=10
            exit
        end if
        if ( visual0(iix,iiy) == 1 .or. visual0(iix,iiy) == 2 .or. visual0(iix,iiy) == 3) then
            flag=9
        endif
        if ( visual0(iix,iiy) == 20 ) then
            flag=flag+1
        endif
        if ( visual0(iix,iiy) == 25) then
            flag=2
        endif
        dval=flwdir(iix,iiy)
        call next_D8(dval,dx,dy)
        iix = iix + dx 
        iiy = iiy + dy
        ! print*, "next out: ", flag, iix, iiy
    end do
    if (flag==2) then
        jx=iix
        jy=iiy 
    else
        jx=-9
        jy=-9
    end if
    return
    end subroutine next_out
    !*****************************************************************
    subroutine got_to_next_tile(ix0,iy0,nx,ny,west0,south0,hiresmap,flwdir,visual,west,south,ix,iy)
    implicit none
    integer                         :: nx, ny
    integer*1,dimension(nx,ny)      :: flwdir, visual
    integer                         :: ix0, iy0
    real                            :: west0, south0
    character*128                   :: hiresmap
    integer                         :: ix, iy
    !--
    integer                         :: dval, ios
    real                            :: dlon, dlat
    character*128                   :: rfile, cname
    real                            :: west, south
    !----------------------------------------
    call find_next_tile(ix0,iy0,nx,ny,dval,ix,iy,dlon,dlat)
    west=west0+dlon 
    south=south0+dlat
    call set_name(west,south,cname)

    ! print*, "go to next tile", west, south, trim(cname)
    ! open files
    rfile=trim(hiresmap)//trim(cname)//'.visual.bin'
    open(21,file=rfile,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) visual
    close(21)
    endif

    rfile=trim(hiresmap)//trim(cname)//'.flwdir.bin'
    open(21,file=rfile,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
    read(21,rec=1) flwdir
    close(21)
    endif
    end subroutine got_to_next_tile
    !*****************************************************************
    subroutine find_next_tile(ix0,iy0,nx,ny,dval,ix,iy,dlon,dlat)
    implicit none
    integer                         :: nx, ny
    integer                         :: ix0, iy0
    integer                         :: dval
    real                            :: dlon, dlat
    integer                         :: ix, iy
    real                            :: dlonlat
    !---------------------------------------------
    dlonlat=10.0
    dlon=0.0
    dlat=0.0
    ! find the next tile according to ix,iy values
    if ( ix0 > 1 .and. ix0 <= nx .and. iy0 < 1 ) then
        dval=1
        dlon=dlonlat*0
        dlat=dlonlat*1
        ix=ix0
        iy=iy0+ny
    else if ( ix0 > nx .and. iy0 < 1 ) then
        dval=2
        dlon=dlonlat*1
        dlat=dlonlat*1
        ix=ix0-nx
        iy=iy0+ny
    else if ( ix0 > nx .and. iy0 > 1 .and. iy0 <= ny ) then
        dval=3
        dlon=dlonlat*1
        dlat=dlonlat*0
        ix=ix0-nx
        iy=iy0
    else if ( ix0 > nx .and. iy0 > ny ) then
        dval=4
        dlon=dlonlat*1
        dlat=dlonlat*(-1)
        ix=ix0-nx
        iy=iy0-ny
    else if ( ix0 > 1 .and. ix0 <= nx .and. iy0 > ny ) then
        dval=5
        dlon=dlonlat*0
        dlat=dlonlat*(-1)
        ix=ix0
        iy=iy0-ny
    else if ( ix0 < 1 .and. iy0 > ny ) then
        dval=6
        dlon=dlonlat*(-1)
        dlat=dlonlat*(-1)
        ix=ix0+nx
        iy=iy0-ny
    else if ( ix0 < 1 .and. iy0 > 1 .and. iy0 <= ny ) then
        dval=7
        dlon=dlonlat*(-1)
        dlat=dlonlat*0
        ix=ix0+nx
        iy=iy0
    else if ( ix0 < 1 .and. iy0 < 1 ) then
        dval=7
        dlon=dlonlat*(-1)
        dlat=dlonlat*1 
        ix=ix0+nx
        iy=iy0+ny   
    end if
    end subroutine find_next_tile  
    !*****************************************************************