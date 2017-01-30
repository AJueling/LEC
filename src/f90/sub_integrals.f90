!===============================================================================
subroutine vol_int(imin,jmin,kmin,imax,jmax,kmax,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field
!
implicit none

! input/output variables
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,    dimension(imt,jmt,km), intent(in)  :: FIELD
real,    dimension(imt,jmt   ), intent(in)  :: TAREA
real,    dimension(imt,jmt,km), intent(in)  :: DZT
real,                           intent(out) :: INTEGRAL

INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)&
         * sum(DZT(imin:imax,jmin:jmax,kmin:kmax)&
         * FIELD(imin:imax,jmin:jmax,kmin:kmax),3,DZT.ne.0.0))

end subroutine vol_int
!===============================================================================

!===============================================================================
subroutine vol_int2(imin,jmin,kmin,imax,jmax,kmax,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field via Kahan summation
!
implicit none

! input/output variables
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,    dimension(imt,jmt,km), intent(in)  :: FIELD
real,    dimension(imt,jmt   ), intent(in)  :: TAREA
real,    dimension(imt,jmt,km), intent(in)  :: DZT
real,                           intent(out) :: INTEGRAL
real                                        :: c, y, t

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

INTEGRAL = INTEGRAL

end subroutine vol_int2
!===============================================================================

!===============================================================================
subroutine surf_int(imin,jmin,imax,jmax,FIELD,TAREA,MASK,INTEGRAL)
!
!     calculates surface integral (*[m^2]) for TT-surface, using Kahan Summation
!
implicit none

integer, parameter                       :: imt=3600,jmt=2400,km=42
integer,                     intent(in)  :: imin,jmin,imax,jmax
real,    dimension(imt,jmt), intent(in)  :: FIELD, MASK!=DZT(:,:,1)
real,    dimension(imt,jmt), intent(in)  :: TAREA
real,                        intent(out) :: INTEGRAL

INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)*FIELD(imin:imax,jmin:jmax),          &
               MASK(imin:imax,jmin:jmax).ne.0.0)/1.0E04

end subroutine surf_int
!===============================================================================

!===============================================================================
subroutine area_avg(FIELD,DZT,AREA,area_k,avg)
! creates unweighted area average

real, dimension(:,:,:), intent(in)  :: FIELD, DZT
real, dimension(:,:),   intent(in)  :: AREA
real, dimension(:),     intent(in)  :: area_k
real, dimension(:),     intent(out) :: avg

do k=1,km
  avg(k) = sum(FIELD(:,:,k)*AREA,DZT(:,:,k).ne.0.0)/area_k(k)
enddo

end subroutine area_avg
!===============================================================================

!===============================================================================
subroutine area_avg_weighted(FIELD,DZT,AREA,vol,w_avg)
! creates cell depth weighted area average

real, dimension(:,:,:), intent(in)  :: FIELD, DZT
real, dimension(:,:),   intent(in)  :: AREA
real, dimension(:),     intent(in)  :: vol
real, dimension(:),     intent(out) :: w_avg

do k=1,km
  w_avg(k) = sum(FIELD(:,:,k)*AREA(:,:)*DZT(:,:,k),DZT(:,:,k).ne.0.0)/vol(k)
enddo

end subroutine area_avg_weighted
!===============================================================================

!===============================================================================
subroutine masked_avg(imt,jmt,DXT,DYT,MASK,mask_area,FIELD,value)
implicit none
!
!  reads yearly data that was previoously calculated and written into a binary
!

integer, intent(in) :: imt,jmt
double precision, dimension(imt,jmt), intent(in) :: DXT,DYT
real, intent(in) :: mask_area
real, dimension(imt,jmt), intent(in) :: FIELD
integer, dimension(imt,jmt), intent(in) :: MASK
real, intent(out) :: value

value = sum(FIELD*DXT*DYT*MASK)/mask_area

end subroutine masked_avg
!===============================================================================

!===============================================================================
subroutine vert_int(FIELD,DZT,INTEGRAL)
!
!     calculates vertical integral (*[m]) for TT-surface, using Kahan Summation
!
implicit none

! input/output variables
real,    dimension(imt,jmt,km), intent(in)  :: FIELD
real,    dimension(imt,jmt,km), intent(in)  :: DZT
real,    dimension(imt,jmt),    intent(out) :: INTEGRAL

INTEGRAL = sum(DZT*FIELD,3)/1.0e02
! DZT in [cm], factor 1e-2 used to achieve [m]

end subroutine vert_int
!===============================================================================

!===============================================================================
subroutine zonal_int(imt,jmt,DXT,FIELD,FIELD_zint)
!
!  "zonal" integral, true for Southern Hemisphere, Northern one is distorted
!
implicit none

integer,                                intent(in)  :: imt,jmt
double precision,   dimension(imt,jmt), intent(in)  :: DXT
real,               dimension(imt,jmt), intent(in)  :: FIELD
real,               dimension(jmt),     intent(out) :: FIELD_zint

FIELD_zint(:) = sum(FIELD(:,:)*DYT(:,:),1)

end subroutine
!===============================================================================

