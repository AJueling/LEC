program geometry
implicit none

!===============================================================================
!
!  this script calculates the geometrical fields used in the calculation of the 
!  Lorenz Energy Cycle in the high resolution POP model (0.1 deg)
! DXT, DXU
! DYT, DYU
! DZT, DZU
! TAREA, UAREA
!  
!===============================================================================

!===============================================================================
!  variables
!===============================================================================

character*120 :: &
  input_folder,grid_file,kmt_file,in_depths,pbc_file,geometry1_file,geometry2_file

integer :: imt,jmt,km,rec_length,i,j,k
integer, dimension(:,:), allocatable :: kmT

double precision :: volume
double precision, dimension(:),   allocatable :: z1, z2,tdz
double precision, dimension(:,:), allocatable ::                            &
HTN, HTE, WORK, DXT, DYT, TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK2,    &
WORK3  
double precision, parameter ::                                              &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

real, dimension(:),     allocatable :: area,dz,tdepth,p_z,vol
real, dimension(:,:,:), allocatable :: DZT, DZU


imt            = 3600
jmt            = 2400
km             =   42

input_folder   = '/home/dijkbio/andre/LEC/input/'
grid_file      = trim(input_folder)//'grid.3600x2400.fob.da'
kmt_file       = trim(input_folder)//'kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'
in_depths      = trim(input_folder)//'in_depths.42.dat'
pbc_file       = trim(input_folder)//'dzbc_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'

geometry1_file  = trim(input_folder)//'geometry1'
geometry2_file  = trim(input_folder)//'geometry2'

write (*,*) ''
write (*,*) '--- GEOMETRY ---'
write (*,*) ''

!===============================================================================
!  read and create horizontal grid spacing, define TAREA
!===============================================================================

allocate(                                                                      &
  HTN(imt,jmt), HTE(imt,jmt), HUW(imt,jmt), HUS(imt,jmt), WORK(imt,jmt),       &
  WORK2(imt,jmt), WORK3(imt,jmt), DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt),  &
  DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt) )

inquire (iolength = rec_length) HTN
!write   (*,*) 'HTN  rec_length', rec_length, imt,jmt
open    (1,file=grid_file,access='direct',form='unformatted',recl=rec_length,  &
           status='unknown')
  read  (1,rec=3) HTN ! [cm]
  read  (1,rec=4) HTE ! [cm]
close   (1)

where (HTN <= c0) HTN = c1
where (HTE <= c0) HTE = c1

call s_rshift(WORK,HTN,imt,jmt)
DXT   = p5*(HTN + WORK)
call w_rshift(WORK,HTE,imt,jmt)
DYT   = p5*(HTE + WORK)
TAREA = DXT*DYT ! [cm^2]

call e_rshift(WORK,HTN,imt,jmt)
DXU   = p5*(HTN + WORK)
call n_rshift(WORK,HTE,imt,jmt)
DYU   = p5*(HTE + WORK)
UAREA = DXU*DYU

! create HUW and HUS (as straight averages of surrrounding HTE/HTN according
! to Reference Manual)
call n_rshift(WORK,HTE,imt,jmt)
call w_rshift(WORK2,WORK,imt,jmt)
call w_rshift(WORK3,HTE,imt,jmt)
HUW   = p5*p5 * ( HTE + WORK + WORK2 + WORK3 ) ! [cm]
call s_rshift(WORK,HTN,imt,jmt)
call e_rshift(WORK2,WORK,imt,jmt)
call e_rshift(WORK3,HTN,imt,jmt)
HUS   = p5*p5 * ( HTN + WORK + WORK2 + WORK3 ) ! [cm]

!   deallocate( WORK, WORK2, WORK3 )

!===============================================================================
!  read bathymetry and depth level variables, define DZT and DZU
!===============================================================================

allocate( dz(km), kmT(imt,jmt), DZBC(imt,jmt), DZT(imt,jmt,km), DZU(imt,jmt,km),&
          tdz(km), tdepth(km), area(km), p_z(km), vol(km) )

! kmT
inquire (iolength=rec_length) kmT
!write   (*,*) 'kmT  rec_length', rec_length, imt,jmt
open    (1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,   &
           status='unknown')
  read  (1,rec=1) kmT
close   (1)
!write   (*,*) 'kmT file: ', kmt_file

! dz
open     (1,file=in_depths,status='old')
  do k = 1, km
    read (1,*) dz(k) ! [cm]
  enddo  !k
close    (1)
!write   (*,*) 'dz(k) file: ', in_depths

! partial bottom cell depths
inquire (iolength=rec_length) DZBC
!write   (*,*) 'DZBC rec_length', rec_length, imt,jmt
open    (1,file=pbc_file,access='direct',form='unformatted',recl=rec_length,   &
           status='unknown')
  read  (1,rec=1) DZBC ! [cm]
close   (1)
!write   (*,*) 'PBC file: ', pbc_file

! create DZT-field, also used as 3D TTT-mask
do k=1,km
  where(kmT>k)
    DZT(:,:,k) = dz(k)
  elsewhere(kmT==k)
    DZT(:,:,k) = DZBC
  elsewhere
    DZT(:,:,k) = c0
  endwhere
enddo !k
  
! DZU field
do k=1,km
  do j=1,jmt-1   ! eastern edge
      DZU(imt,j,k) = min(DZT(imt,j,k),DZT(1,j,k),DZT(imt,j+1,k),DZT(1,j+1,k))
    do i=1,imt-1  ! rest of array
      if ( k.le.kmT(i,j) ) then
        DZU(i,j,k) = min(DZT(i,j,k),DZT(i+1,j,k),DZT(i,j+1,k),DZT(i+1,j+1,k))
      endif
    enddo !j
  enddo !j
  ! northern edge
  do i=1,imt-1
    if ( i.lt.imt/2 ) then      ! first half of array, before wrapping around
      DZU(i,j,k) = min( DZT(i,jmt,k),      &
                        DZT(i+1,jmt,k),    &
                        DZT(imt-i,jmt,k),  &
                        DZT(imt-i+1,jmt,k) )   
    elseif ( i.gt.imt/2 ) then  ! second half, after wrap
      DZU(i,jmt,k) = DZU(imt-i,jmt,k)
    endif
  enddo !i
enddo !k

!  create tdepth 1D array containing vertical extent of T-cells
open(1,file=in_depths,status='old')
do k = 1, km
  read (1,*) tdz(k)           ! [cm]
  dz(k) = real(tdz(k))/1.0E02 ! to [m]
  !write(*,*) dz(k)
enddo  !k
close(1)

! depths of T-cells
tdepth(1) = dz(1)/2.0
do k = 2,km
  tdepth(k) = tdepth(k-1)+p5*(dz(k-1)+dz(k)) ! [m]
  !write(*,*) tdepth(k)
enddo !k

! (T)area per level [m^2]
do k = 1, km
  area(k) = sum(TAREA/1.0E04,DZT(:,:,k).ne.0) ! [m^2]
  p_z(k) = pressure(tdepth(k))         ! [bar]
  !write(*,*) area(k), p_z(k)
enddo

! volume per level [m^3] including partial bottom cells
do k=1,km
  vol(k) = sum(DZT(:,:,k)*TAREA(:,:),DZT(:,:,k).ne.0.0)/1.0E06
enddo

!===============================================================================
!  OUTPUT
!===============================================================================

! DXT, DXU
! DYT, DYU
! TAREA, UAREA
! DZT, DZU


open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,   &
           status='unknown')

write(1,rec=1) real(DXT/1.0E02)   ! [m]
write(1,rec=2) real(DXU/1.0E02)   ! [m]
write(1,rec=3) real(DYT/1.0E02)   ! [m]
write(1,rec=4) real(DYU/1.0E02)   ! [m]
write(1,rec=5) real(TAREA/1.0E04) ! [m^2]
write(1,rec=6) real(UAREA/1.0E04) ! [m^2]

do k=1,km
  write(1,rec=6+   k) DZT(:,:,k)/1.0E02 ! [m]
  write(1,rec=6+km+k) DZU(:,:,k)/1.0E02 ! [m]
enddo

close(1)

! k, dz, tdepth, area, p, vol
open(1,file=geometry2_file,access='direct',form='unformatted',recl=6,         &
       status='unknown')
do k=1,km
  write(1,rec=k) real(k),dz(k),tdepth(k),area(k),p_z(k),vol(k)
enddo
close(1)

!===============================================================================
!===============================================================================
contains

!===============================================================================
subroutine load_3D_field(imt,jmt,km,nfile,nrec,FIELD)
!
!  loads 3D field (starting with record number nrec) from file with number nfile
!
implicit none

integer, intent(in)                         :: imt, jmt, km, nfile,  nrec
real,    dimension(imt,jmt)                 :: WORK
real,    dimension(imt,jmt,km), intent(out) :: FIELD

do k = 1,km
  read (nfile,rec=nrec+k-1) WORK
  FIELD(:,:,k) = WORK
enddo

end subroutine load_3d_field
!===============================================================================

!===============================================================================
subroutine load_2D_field(imt,jmt,nfile,nrec,FIELD)
!
!  loads 2D field (starting with record number nrec) from file with number nfile
!
implicit none

integer, intent(in)                      :: imt, jmt, nfile,  nrec
real,    dimension(imt,jmt), intent(out) :: FIELD

read (nfile,rec=nrec) FIELD

end subroutine load_2D_field
!===============================================================================

!===============================================================================
subroutine surf_int(imin,jmin,imax,jmax,FIELD,TAREA,MASK,INTEGRAL)
!
!     calculates surface integral (*[m^2]) for TT-surface, using Kahan Summation
!
implicit none

! input/output variables
integer,                                 intent(in)  :: imin,jmin,imax,jmax
real,                dimension(imt,jmt), intent(in)  :: FIELD, MASK!=DZT(:,:,1)
double precision,    dimension(imt,jmt), intent(in)  :: TAREA 
real,                                    intent(out) :: INTEGRAL
real                                                 :: c, y, t

! INTEGRAL = 0.0
! c = 0.0
! do j = jmin,jmax
!  do i = imin,imax
!   if ( MASK(i,j).ne.0.0 ) then
!    y = TAREA(i,j) * FIELD(i,j) - c
!    t = INTEGRAL + y
!    c = ( t - INTEGRAL ) - y
!    INTEGRAL = t
!   endif
!  enddo !i
! enddo !j
! 
! INTEGRAL = INTEGRAL/1.0E04

!      INTEGRAL = sum(TAREA*FIELD, MASK.ne.0.0)/1.0E04
INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)*FIELD(imin:imax,jmin:jmax),          &
               MASK(imin:imax,jmin:jmax).ne.0.0)/1.0E04 

end subroutine surf_int
!===============================================================================

!===============================================================================
subroutine vol_int(imin,jmin,kmin,imax,jmax,kmax,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field, using Kahan Summation
!
implicit none

! input/output variables
integer,             intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,                dimension(imt,jmt,km), intent(in)  :: FIELD
double precision,    dimension(imt,jmt   ), intent(in)  :: TAREA
real,                dimension(imt,jmt,km), intent(in)  :: DZT
real,                                       intent(out) :: INTEGRAL
real                                                    :: c, y, t

INTEGRAL = 0.0
c = 0.0
do k = kmin,kmax
  do j = jmin,jmax
    do i = imin,imax
      if ( DZT(i,j,k).ne.0.0 ) then
        y = TAREA(i,j) * DZT(i,j,k) * FIELD(i,j,k) - c
        t = INTEGRAL + y
        c = ( t - INTEGRAL ) - y
        INTEGRAL = t
      endif
    enddo !i
  enddo !j
enddo !k

INTEGRAL = INTEGRAL/1.0E06

!INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)                                      &
!         * sum(DZT(imin:imax,jmin:jmax,kmin:kmax)                              &
!         * FIELD(imin:imax,jmin:jmax,kmin:kmax),3))/1.0E06

end subroutine vol_int
!===============================================================================

!===============================================================================
subroutine e_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from east (double precision arrays)
!-----------------------------------------------------------------------
!     inputs
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
!     outputs
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result
!     local variables
integer :: i,j,n,ip

do j=1,jmt
  do i=1,imt
    ip = i + 1
    if(i == imt) ip = 1
    XOUT(i,j) = X(ip,j)
  enddo
enddo

end subroutine e_rshift
!===============================================================================

!===============================================================================
subroutine s_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from south (real arrays)
!-----------------------------------------------------------------------
!     inputs
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
!     outputs
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result
!     local variables
integer :: i,j,n,jm

do j=1,jmt
  jm = j - 1
  if(j == 1) jm = 1
  do i=1,imt
    XOUT(i,j) = X(i,jm)
  enddo
enddo

end subroutine s_rshift
!===============================================================================

!===============================================================================
subroutine w_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from west (real arrays)
!-----------------------------------------------------------------------
!     inputs
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
!     outputs
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result
!     local variables
integer :: i,j,n,im

do j=1,jmt
  do i=1,imt
    im = i - 1
    if(i == 1) im = imt
    XOUT(i,j) = X(im,j)
  end do
end do

end subroutine w_rshift
!===============================================================================

!===============================================================================
subroutine n_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from north (double precision arrays)
!-----------------------------------------------------------------------
!     inputs
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
!     outputs
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result
!     local variables
integer :: i,j,n,jp

do j=1,jmt
  jp = j+1
  if(j == jmt) jp = jmt
  do i=1,imt
    XOUT(i,j) = X(i,jp)
  end do
end do

end subroutine n_rshift
!===============================================================================

!===============================================================================
function pressure(depth)

! !DESCRIPTION:
!  This function computes pressure in bars from depth in meters
!  using a mean density derived from depth-dependent global 
!  average temperatures and salinities from Levitus 1994, and 
!  integrating using hydrostatic balance.
!
!  References:
!
!     Levitus, S., R. Burgett, and T.P. Boyer, World Ocean Atlas 
!          1994, Volume 3: Salinity, NOAA Atlas NESDIS 3, US Dept. of 
!          Commerce, 1994.
!
!     Levitus, S. and T.P. Boyer, World Ocean Atlas 1994, 
!          Volume 4: Temperature, NOAA Atlas NESDIS 4, US Dept. of 
!          Commerce, 1994.
!
!     Dukowicz, J. K., 2000: Reduction of Pressure and Pressure
!          Gradient Errors in Ocean Simulations, J. Phys. Oceanogr.,
!          submitted.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

real, intent(in) :: depth    ! depth in meters

! !OUTPUT PARAMETERS:

real :: pressure   ! pressure [bars]

! !LOCAL PARAMETERS:
!  adjustment from original code which was written as 0.1_r8
!  which the compiler disliked so constants are now pre-defined here
double precision, parameter ::              &
   pc1 = 0.059808,   &
   pc2 = 0.025,      &
   pc3 = 0.100766,   &
   pc4 = 2.28405e-7, &
   c1  = 1.0
!EOP
!BOC
!-----------------------------------------------------------------------
!  convert depth in meters to pressure in bars
!-----------------------------------------------------------------------

pressure = pc1*(exp(-pc2*depth) - c1) + pc3*depth + pc4*depth**2

end function pressure
!===============================================================================

!===============================================================================
!===============================================================================
end program geometry
