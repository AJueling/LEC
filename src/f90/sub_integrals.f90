!===============================================================================
subroutine vol_int_part(imin,jmin,kmin,imax,jmax,kmax,FIELD,AREA,dz,DZ_3D,INTEGRAL)
!
!     calculates volume integral (*[m^3]) with partial bottom cells
!
implicit none

! input/output variables
integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,    dimension(km),         intent(in)  :: dz
real,    dimension(imt,jmt,km), intent(in)  :: DZ_3D,FIELD
real,    dimension(imt,jmt),    intent(in)  :: AREA
real,                           intent(out) :: INTEGRAL

INTEGRAL = sum( AREA(imin:imax,jmin:jmax)                                      &
                * sum( DZ_3D(imin:imax,jmin:jmax,kmin:kmax)                    &
                       * FIELD(imin:imax,jmin:jmax,kmin:kmax),                 &
               3,DZ_3D(imin:imax,jmin:jmax,kmin:kmax).ne.0.0 ) )

end subroutine vol_int_part
!===============================================================================

!===============================================================================
subroutine vol_int_full(imin,jmin,kmin,imax,jmax,kmax,FIELD,AREA,dz,DZ_3D,INTEGRAL)
!
!     calculates volume integral (*[m^3]) with full bottom cells
!
implicit none

! input/output variables
integer                                     :: k
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,    dimension(km),         intent(in)  :: dz
real,    dimension(imt,jmt,km), intent(in)  :: DZ_3D,FIELD
real,    dimension(imt,jmt),    intent(in)  :: AREA
real,                           intent(out) :: INTEGRAL

INTEGRAL = 0.0
do k = kmin,kmax
  INTEGRAL = INTEGRAL                                                          &
           + sum( AREA(imin:imax,jmin:jmax)*FIELD(imin:imax,jmin:jmax,k),      &
                  DZ_3D(imin:imax,jmin:jmax,k).ne.0.0 ) * dz(k)
enddo !k

end subroutine vol_int_full
!===============================================================================
!===============================================================================
subroutine vol_int2(imin,jmin,kmin,imax,jmax,kmax,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field via Kahan summation
!
implicit none

! input/output variables
integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,    dimension(imt,jmt,km), intent(in)  :: DZT,FIELD
real,    dimension(imt,jmt),    intent(in)  :: TAREA
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
subroutine surf_int_3D(imin,jmin,imax,jmax,FIELD,TAREA,DZT_k,INTEGRAL)
!  calculates surface integrals *[m^2] weighting contribution by (pbc)cell-depth
implicit none

integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                     intent(in)  :: imin,jmin,imax,jmax
real,    dimension(imt,jmt), intent(in)  :: FIELD,TAREA,DZT_k
real,                        intent(out) :: INTEGRAL
real                                     :: avg_dz

avg_dz   = sum( TAREA(imin:imax,jmin:jmax) * DZT_k(imin:imax,jmin:jmax),       &
                DZT_k(imin:imax,jmin:jmax).ne.0.0  ) /    &
           sum( TAREA(imin:imax,jmin:jmax),  DZT_k(imin:imax,jmin:jmax).ne.0.0 )
INTEGRAL = sum( TAREA(imin:imax,jmin:jmax) * FIELD(imin:imax,jmin:jmax)        &
              * DZT_k(imin:imax,jmin:jmax),  DZT_k(imin:imax,jmin:jmax).ne.0.0)
INTEGRAL = INTEGRAL/avg_dz

end subroutine surf_int_3D
!===============================================================================

!===============================================================================
subroutine surf_int_2D(imin,jmin,imax,jmax,FIELD,TAREA,MASK,INTEGRAL)
!
!     calculates surface integral (*[m^2]) for TT-surface, using Kahan Summation
!
implicit none

integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                     intent(in)  :: imin,jmin,imax,jmax
real,    dimension(imt,jmt), intent(in)  :: FIELD,TAREA,MASK!=DZT(:,:,1)
real,                        intent(out) :: INTEGRAL

INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)*FIELD(imin:imax,jmin:jmax),          &
               MASK(imin:imax,jmin:jmax).ne.0.0)

end subroutine surf_int_2D
!===============================================================================

!===============================================================================
subroutine area_avg(FIELD,DZT,AREA,area_k,avg)
! creates unweighted area average

integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt,km), intent(in)  :: FIELD, DZT
real,    dimension(imt,jmt),    intent(in)  :: AREA
real,    dimension(km),         intent(in)  :: area_k
real,    dimension(km),         intent(out) :: avg

do k=1,km
  avg(k) = sum(FIELD(:,:,k)*AREA(:,:),DZT(:,:,k).ne.0.0)/area_k(k)
enddo

end subroutine area_avg
!===============================================================================

!===============================================================================
subroutine area_avg_weighted(FIELD,DZT,AREA,vol,w_avg)
! creates cell depth weighted area average

integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt,km), intent(in)  :: FIELD, DZT
real,    dimension(imt,jmt),    intent(in)  :: AREA
real,    dimension(km),         intent(in)  :: vol
real,    dimension(km),         intent(out) :: w_avg

do k=1,km
  w_avg(k) = sum(FIELD(:,:,k)*AREA(imt,jmt)*DZT(:,:,k),DZT(:,:,k).ne.0.0)/vol(k)
enddo

end subroutine area_avg_weighted
!===============================================================================

!===============================================================================
subroutine masked_avg(DXT,DYT,MASK,mask_area,FIELD,value)
implicit none
!
!  reads yearly data that was previoously calculated and written into a binary
!

integer                                  :: i,j,k
integer, parameter                       :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt), intent(in)  :: DXT,DYT,FIELD
real,                        intent(in)  :: mask_area
integer, dimension(imt,jmt), intent(in)  :: MASK
real,                        intent(out) :: value

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
integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt,km), intent(in)  :: DZT, FIELD
real,    dimension(imt,jmt),    intent(out) :: INTEGRAL

INTEGRAL = sum(DZT*FIELD,3)
! DZT in [cm], factor 1e-2 used to achieve [m]

end subroutine vert_int
!===============================================================================

!===============================================================================
subroutine zonal_int(DXT,FIELD,FIELD_zint)
!
!  "zonal" integral, true for Southern Hemisphere, Northern one is distorted
!
implicit none

integer                                  :: i,j,k
integer, parameter                       :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt), intent(in)  :: DXT, FIELD
real,    dimension(jmt),     intent(out) :: FIELD_zint

FIELD_zint(:) = sum(FIELD(:,:)*DXT(:,:),1)

end subroutine
!===============================================================================

