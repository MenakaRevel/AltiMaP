      ! f2py -c -m river_function river_function.F90 --fcompiler=gnu95
      ! python -m numpy.f2py -c -m river_function river_function.F90
!***************************************************************************
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
!***************************************************************************
    subroutine next_D8(dval,dx,dy)
    implicit none
    integer                         :: dval
    integer                         :: dx, dy
    !real                            :: tval
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
    real                            :: west0, south0, north0 !, tval
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
    subroutine upstream(i,j,nx,ny,flwdir,uparea,visual,x,y)
    implicit none 
    ! find the upstream pixel with closest upstream area to the i,j
    !--
    integer,intent(IN)                        :: i,j,nx,ny
    integer*1,dimension(nx,ny),intent(IN)     :: flwdir, visual !,rivseq
    real,dimension(nx,ny),intent(IN)          :: uparea
    integer,intent(OUT)                       :: x,y
    !--
    real                                      :: dA ! area differnce nextdst,
    integer                                   :: ix,iy,tx,ty,d !,iix,iiy,ud
    !real                                      :: length,rl
    !--
    x=-9999
    y=-9999
    d=3 ! look at 3*3 box
    dA=1.0e20 ! area differnce
    !--
    !write(*,*)i,j
    do tx=i-d,i+d
      do ty=j-d,j+d
        ! print*, "L192: ",tx, ty
        if (tx<=0 .or. ty<=0) cycle
        if (tx>nx .or. ty>ny) cycle
        if (tx==i .and. ty==j) cycle
        if (visual(tx,ty) ==10 .or. visual(tx,ty) ==20 ) then
        !write(*,*)tx,ty
        ! call ixy2iixy(tx,ty, nx, ny, ix, iy)
        !write(*,*)nextX(ix,iy),nextY(ix,iy),ix,iy,uparea(ix,iy),rivseq(ix,iy)
            call next_pixel(tx,ty,flwdir,nx,ny,ix,iy)
            if (ix == i .and. iy == j) then
                ! print*, tx, ty
                if (abs(uparea(i,j)-uparea(tx,ty)) < dA) then
                    dA=abs(uparea(i,j)-uparea(tx,ty))
                    x=tx
                    y=ty
                    ! print*, x,y
                end if
                ! end if
            end if
        end if
      end do
    end do
    ! print*, "L211: upstream ->",x,y
    return
    end subroutine upstream
!*****************************************************************
    subroutine next_pixel(tx,ty,flwdir,nx,ny,ix,iy)
    implicit none
    integer                                   :: tx,ty,ix,iy !iix,iiy,ud,d
    integer                                   :: nx, ny
    integer*1,dimension(nx,ny)                :: flwdir
    integer                                   :: dval, dx, dy
    !----
    dval=flwdir(tx,ty)
    call next_D8(dval,dx,dy)
    ix=tx+dx 
    iy=ty+dy 
    return
    end subroutine next_pixel
!*****************************************************************
    subroutine up_until_mouth(ix,iy,flwdir,uparea,visual,nx,ny,x,y)
    implicit none
    ! find the upstream unit catchment mouth 
    !--
    integer,intent(IN)                        :: ix,iy,nx,ny
    integer*1,dimension(nx,ny),intent(IN)     :: flwdir, visual !,rivseq
    real,dimension(nx,ny),intent(IN)          :: uparea
    integer,intent(OUT)                       :: x,y
    !--
    !real                                      :: dA ! area differnce nextdst,
    integer                                   :: iix,iiy,x0,y0 !,tx,ty,ud,d
    !real                                      :: length,rl
    !------------
    x=-9999
    y=-9999
    iix=ix
    iiy=iy 
    ! print*, iix, iiy, visual(iix,iiy)
    do while (visual(iix,iiy)==10)
        ! print*, "upstream", iix, iiy, visual(iix,iiy)
        call upstream(iix,iiy,nx,ny,flwdir,uparea,visual,x0,y0)
        iix=x0
        iiy=y0
        ! print*, "L252: ",iix, iiy
        if (iix == -9999 .or. iiy == -9999) then
            x=-9999
            y=-9999
            exit
        end if
        ! print*, "at upstream:",visual(iix,iiy)
        if (visual(iix,iiy) == 20) then ! upstream unit-catchment mouth found
            x=x0
            y=y0
            exit
        end if
    end do
    return
    end subroutine up_until_mouth
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
    subroutine river_profile(ix,iy,west,south,csize,nx,ny,hiresmap,len,ele,k)
    implicit none
    ! output the river profile along the river channel from upstream unit-catchment mouth
    ! to unit-catchment mouth of this catchment
    character*128,intent(IN)                  :: hiresmap
    integer,intent(IN)                        :: ix,iy,nx,ny
    real,intent(IN)                           :: west,south,csize
    real,dimension(1000),intent(OUT)          :: len, ele
    integer,intent(OUT)                       :: k
    !- for clculations
    character*128                             :: rfile1
    character*7                               :: cname
    integer                                   :: ios
    ! integer,parameter                         :: nx=1200, ny=1200
    integer*1,dimension(nx,ny)                :: flwdir, visual !,rivseq
    real,dimension(nx,ny)                     :: uparea, elevtn
    integer                                   :: iix,iiy,x0,y0,x,y !,tx,ty,ud,d
    integer                                   :: dval, dx,dy
    real                                      :: down_dist
    real                                      :: west0, south0, north0 !, tval
    real                                      :: lon1, lat1, lon2, lat2
    integer*1,dimension(nx,ny)                :: flwdir0, visual0
    real,dimension(nx,ny)                     :: uparea0
    real                                      :: hubeny_real
    integer                                   :: flag
    !========
    call set_name(west,south,cname)
    
    ! print*, "open high-resolution maps  ",ix, iy, " ", trim(cname)!, "  ",trim(hiresmap)
    ! open high-resolution maps
    rfile1=trim(hiresmap)//trim(cname)//'.elevtn.bin'
    ! print*, "open: ",trim(rfile1)
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
        read(21,rec=1) elevtn
        close(21)
    else
        print*, "no file: ", rfile1
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.uparea.bin'
    ! print*, "open: ",trim(rfile1)
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=4*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
        read(21,rec=1) uparea
        close(21)
    else
        print*, "no file: ", rfile1
    endif

    rfile1=trim(hiresmap)//trim(cname)//'.visual.bin'
    ! print*, "open: ",trim(rfile1)
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
        read(21,rec=1) visual
        close(21)
    else
        print*, "no file: ", rfile1
    endif
    ! print*, "visual:" , visual(ix,iy), visual(iy,ix), ios
    rfile1=trim(hiresmap)//trim(cname)//'.flwdir.bin'
    ! print*, "open: ",trim(rfile1)
    open(21,file=rfile1,form='unformatted',access='direct' , action='READ',recl=1*nx*ny,status='old',iostat=ios)
    if( ios==0 )then
        read(21,rec=1) flwdir
        close(21)
    else
        print*, "no file: ", rfile1
    endif

    
    !--
    ! print*, "intialize"
    flwdir0=flwdir
    visual0=visual
    uparea0=uparea
    ! print*, "coordinates ",west, south
    west0=west+csize/2.0
    ! print*, "L484"
    south0=south+csize/2.0
    north0=south+10.0+csize/2.0
    ! print*, "L487"
    down_dist = 0.0
    ! print*, "L489"
    ! iix=ix
    ! iiy=iy
    ! print*, "start calculation............"!, iix, iiy, visual0(iix,iiy), visual(iix,iiy)
    ! intilize len, ele
    len=-9999.0
    ele=-9999.0

    ! find the upstream unit-catchment mouth
    flag=1
    ! print*, "up_until_mouth", visual0(ix,iy)
    call up_until_mouth(ix,iy,flwdir0,uparea0,visual0,nx,ny,x0,y0)
    ! print* , x0, y0
    if (x0 == -9999 .or. y0 == -9999) then
        flag=-9999
        ! print*, "no upstream location, ", flag
    end if

    if (flag == 1) then
        iix=x0
        iiy=y0
        !==============================
        ! upstream unit-catchment mouth
        len(1)=0.0
        ele(1)=elevtn(iix,iiy)
        lon1=west0+real(iix)*csize
        lat1=north0-real(iiy)*csize
        ! print*, "next_pixel"
        !==============================
        ! go to the next pixel
        call next_pixel(x0,y0,flwdir0,nx,ny,x,y)
        iix=x
        iiy=y
        lon2=west0+real(iix)*csize
        lat2=north0-real(iiy)*csize
        down_dist = hubeny_real(lat1, lon1, lat2, lon2)
        len(2)=down_dist
        ele(2)=elevtn(iix,iiy)
        k=3
        !------------------------------
        lon1=lon2
        lat1=lat2
        ! print*, "Start do-loop:" ,iix,iiy, visual0(iix,iiy)
        do while (visual0(iix,iiy)==10)
            ! print* ,iix,iiy, visual0(iix,iiy)
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
                ! iix0 = iix
                ! iiy0 = iiy 
                ! call got_to_next_tile(iix0,iiy0,nx,ny,west,south,hiresmap,flwdir0,visual0,west0,south0,iix,iiy)
                ! replace west0, south0, flwdir0, visual0
                ! west0=west0+csize/2.0
                ! south0=south0+csize/2.0
                ! north0=south0+10.0+csize/2.0
                ! print*, west, south, "to ",west0, south0
                ! print*, iix0, iiy0, "to ", iix, iiy
                exit
            end if
            lon2=lon1+real(dx)*csize 
            lat2=lat1+real(dy)*csize
            down_dist=down_dist+hubeny_real(lat1, lon1, lat2, lon2)
            len(k)=down_dist
            ele(k)=elevtn(iix,iiy)
            ! print*,k , iix, iiy, down_dist, elevtn(iix,iiy)
            ! print*, iix, iiy, lon1, lat1, down_dist ,visual(iix,iiy), dval
            lon1=lon2
            lat1=lat2
            k=k+1
        end do
        k=k-1 
    else 
        k=0    
    end if
    return
    end subroutine river_profile
    !*****************************************************************