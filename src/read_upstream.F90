      ! f2py -c -m read_patchMS read_upstream.F90 --fcompiler=gnu95
      ! python -m numpy.f2py -c -m read_upstream read_upstream.F90
    !***************************************************************************
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
      integer,intent(IN)                        :: i,j,nx,ny
      integer,dimension(nx,ny),intent(IN)       :: nextX,nextY !,rivseq
      real,dimension(nx,ny),intent(IN)          :: uparea
      integer,intent(OUT)                       ::x,y
      !--
      real                                      :: dA ! area differnce nextdst,
      integer                                   :: ix,iy,iix,iiy,tx,ty,ud,d
      real                                      :: length,rl
      !--
      x=-9999
      y=-9999
      d=100 ! look at 100*100 box
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
              if (abs(uparea(i,j)-uparea(ix,iy)) < dA) then
                  dA=abs(uparea(i,j)-uparea(ix,iy))
                  x=ix
                  y=iy
                  !write(*,*)x,y
              end if
          end if
        end do
      end do
      return
      end subroutine upstream
      !******************************************************************