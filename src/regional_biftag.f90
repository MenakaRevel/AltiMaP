program cut_domain
    ! ==================================================
        implicit none
!
        character*256          ::  region_info
        ! parameter                 (region_info='./region_info.txt')
        character*256          ::  global_dir
        character*256          ::  CaMa_dir
        character*256          ::  regional_map
        character*256          ::  org_map
        character*256          ::  buf
        real                   ::  west, east, north, south
        real                   ::  west2, east2, north2, south2

        character*256          ::  global_param, cut_param
        real                   ::  lon_ori, lat_ori, lon_end, lat_end
        real                   ::  glon, glat, d1, d2

        character*256          ::  region_dir, region_param
        ! parameter                 (region_dir='../')
        ! parameter                 (region_param='../params.txt')
!
        integer                ::  ix, iy, jx, jy, kx, ky
        integer                ::  nx, ny, mx, my, dx, dy
        real                   ::  gsize
        integer                ::  iflp, nflp                  !! number of floodplain layer
!
        character*256          ::  fnextxy0, fdownxy0, felevtn0, ffldhgt0, fctmare0, fgrdare0, fuparea0
        character*256          ::  flonlat0, fnxtdst0, frivlen0, fwidth0, fbiftag0
        character*256          ::  fnextxy,  fdownxy,  felevtn,  ffldhgt,  fctmare,  fgrdare,  fuparea
        character*256          ::  flonlat,  fnxtdst,  frivlen,  fwidth, fbiftag
        character*256          ::  flsmask

        integer,allocatable    ::  nextx0(:,:),  nexty0(:,:)      !!  global maps
        integer,allocatable    ::  downx0(:,:),  downy0(:,:)      !!  global maps
        integer,allocatable    ::  biftag0(:,:)
        real,allocatable       ::  elevtn0(:,:), fldhgt0(:,:,:)
        real,allocatable       ::  ctmare0(:,:), uparea0(:,:), grdare0(:,:)
        real,allocatable       ::  lon0(:,:)   , lat0(:,:)
        real,allocatable       ::  nxtdst0(:,:), rivlen0(:,:)
        real,allocatable       ::  width0(:,:)

        integer,allocatable    ::  nextx(:,:),  nexty(:,:)        !!  regional maps
        integer,allocatable    ::  downx(:,:),  downy(:,:)        !!  regional maps
        integer,allocatable    ::  biftag(:,:)
        real,allocatable       ::  elevtn(:,:), fldhgt(:,:,:)
        real,allocatable       ::  ctmare(:,:), uparea(:,:),   grdare(:,:)
        real,allocatable       ::  lon(:,:),    lat(:,:)
        real,allocatable       ::  nxtdst(:,:), rivlen(:,:)
        real,allocatable       ::  width(:,:)
        integer,allocatable    ::  lsmask(:,:)

        real,allocatable       ::  tmp0(:,:), tmp(:,:)            !! for fldhgt I/O

        integer,allocatable    ::  check(:,:)
        integer                ::  icheck

! ==================================================
        call getarg(1,buf)
        read(buf,"(A)") regional_map
        write(*,*) regional_map

        call getarg(2,buf)
        read(buf,"(A)") org_map
        
        call getarg(3,buf)
        read(buf,"(A)") CaMa_dir
        write(*,*) CaMa_dir

        global_dir=trim(CaMa_dir)//'/map/'//trim(org_map)
        region_dir=trim(CaMa_dir)//'/map/'//trim(regional_map)

        region_param=trim(region_dir)//"/params.txt"

        print*, region_param
        open(11,file=region_param,form='formatted')
        read(11,*) mx                  !! regional map nx*ny
        read(11,*) my
        read(11,*) 
        read(11,*) 
        read(11,*) west               !! regional map west
        read(11,*) east
        read(11,*) south
        read(11,*) north             !! regional map north
        close(11)

        
        global_param=trim(global_dir)//'/params.txt'

        open(11,file=global_param,form='formatted')
        read(11,*) nx                  !! global map nx*ny
        read(11,*) ny
        read(11,*) nflp
        read(11,*) gsize
        read(11,*) lon_ori             !! global map west
        read(11,*) lon_end
        read(11,*) lat_end
        read(11,*) lat_ori             !! global map north
        close(11)

        allocate(nextx0(nx,ny),nexty0(nx,ny),downx0(nx,ny),downy0(nx,ny))
        allocate(elevtn0(nx,ny),fldhgt0(nx,ny,nflp),ctmare0(nx,ny),grdare0(nx,ny),uparea0(nx,ny))
        allocate(lon0(nx,ny),  lat0(nx,ny),  nxtdst0(nx,ny),rivlen0(nx,ny), width0(nx,ny))
        allocate(biftag0(nx,ny))

        print *, 'input:  ', west, east, north, south
        ! west2=west
        ! east2=east
        ! north2=north
        ! south2=south

        !! find optimum west, east, north, south, if the regionalized grid is offset from the original grid
        ! d1=1.e20
        ! d2=1.e20
        ! do ix=1, nx
        ! glon=lon_ori+real(ix-1)*gsize
        ! if( west2>=glon .and. west2< glon+gsize .and. abs(west2-glon)<d1 )then
        !     west=glon
        !     d1=abs(west2-glon)
        ! endif
        ! if( east2> glon .and. east2<=glon+gsize .and. abs(east2-glon-gsize)<d2 )then
        !     east=glon+gsize
        !     d2=abs(east2-glon-gsize)
        ! endif
        ! end do

        ! d1=1.e20
        ! d2=1.e20
        ! do iy=1, ny
        ! glat=lat_ori-real(iy-1)*gsize
        ! if( north2> glat-gsize .and. north2<=glat .and. abs(north2-glat)<d1 )then
        !     north=glat
        !     d1=abs(north2-glat)
        ! endif
        ! if( south2>=glat-gsize .and. south2< glat .and. abs(south2-glat+gsize)<d2 )then
        !     south=glat-gsize
        !     d2=abs(south2-glat+gsize)
        ! endif
        ! end do

        ! !! if differemce ios very small, use the original value
        ! if( abs(west- west2 )<0.001 ) west =west2
        ! if( abs(east- east2 )<0.001 ) east =east2
        ! if( abs(north-north2)<0.001 ) north=north2
        ! if( abs(south-south2)<0.001 ) south=south2

        ! print *, 'output: ', west, east, north, south

        ! mx=nint( dble(east-west)    /dble(gsize) +0.001 )    !!  add 0.001 to avoid rounding error
        ! my=nint( dble(north-south)  /dble(gsize) +0.001 )
        dx=nint( dble(west-lon_ori) /dble(gsize) +0.001 )
        dy=nint( dble(lat_ori-north)/dble(gsize) +0.001 )

        allocate(nextx(mx,my),nexty(mx,my),downx(mx,my),downy(mx,my))
        allocate(elevtn(mx,my),fldhgt(mx,my,nflp),ctmare(mx,my),grdare(mx,my),uparea(mx,my))
        allocate(lon(mx,my),  lat(mx,my),  nxtdst(mx,my),rivlen(mx,my), width(mx,my))
        allocate(lsmask(mx,my))
        allocate(biftag(mx,my))

        allocate(tmp0(nx,ny), tmp(mx,my))

        
        fbiftag0=trim(global_dir)//'/biftag.bin'

    print *, 'read global maps'
        print *, trim(fbiftag0)
        open(11,file=fbiftag0,form='unformatted',access='direct',recl=4*nx*ny)
        read(11,rec=1) biftag0
        close(11)
    
    print *, 'cut domain'

        biftag(:,:)=-9999

        do iy=1, my
            do ix=1, mx
              jx=ix+dx
              jy=iy+dy

              if( nextx0(jx,jy)>0 )then
                nextx(ix,iy)=nextx0(jx,jy)-dx
                nexty(ix,iy)=nexty0(jx,jy)-dy
                downx(ix,iy)=downx0(jx,jy)
                downy(ix,iy)=downy0(jx,jy)
                kx=nextx(ix,iy)
                ky=nexty(ix,iy)
                if( kx<1 .or. kx>mx .or. ky<1 .or. ky>my )then  !! if downstream is outside the domain
                  nextx(ix,iy)=-10
                  nexty(ix,iy)=-10
                  downx(ix,iy)=-1000
                  downy(ix,iy)=-1000
                endif
              elseif( nextx0(jx,jy)/=-9999 )then     !! if river mouth
                nextx(ix,iy)=nextx0(jx,jy)
                nexty(ix,iy)=nexty0(jx,jy)
                downx(ix,iy)=downx0(jx,jy)
                downy(ix,iy)=downy0(jx,jy)
              endif

              if( nextx(ix,iy)/=-9999 )then
                biftag(ix,iy)  =biftag0(jx,jy)
              endif
            end do
        end do

    print *, 'write reagional maps'
        
        fbiftag=trim(region_dir)//'/biftag.bin'

        print *, trim(fbiftag)
        open(21,file=fbiftag,form='unformatted',access='direct',recl=4*mx*my)
        write(21,rec=1) biftag
        close(21)

! ==================================================
    end program cut_domain