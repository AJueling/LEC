program extra_analysis
implicit none

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  this script calculates various quantities from the original POP output files
!  and the LEC files created with LEC.f90
!
!  calculated quantities:
!  1. vertical integrals of cPKm/cPKe
!  2. Eulerian meridional overturning stream function
!
!  generated output
!  1. two binary files with verticaly integrated cPKm/cPKe for each year 
!     + an average of said years
!  2. one binary file with Eul. stream function as an output
!
!  for curvilinear grid:
!  longitude (i) 0=250E, 1100=0E, 1500=40E
!  latitude  (j) 0=78.5S, 518=55S, 677=45S, 866=30S, 1181=0S
!  depth     (k) 1=5m, 9=96.9m, 10=112m, 14=198m, 17=318m, 35=4125, 42=5875m
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!===============================================================================
!  variables
!===============================================================================

! user input 
integer ::                                                                     &
  imt,jmt,km,ntavg,                                                            &
  y, start_year, end_year, nt,                                                 &
  nrec_rPm, nrec_rPe, nrec_rKm, nrec_rKe,                                      &
  nrec_gPm, nrec_gPe, nrec_gKm, nrec_gKe,                                      &
  nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                                  &
  nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                                  &
  nrec_TEMP, nrec_VVEL, nrec_WVEL

character*120 :: grid_file,kmt_file,in_depths,pbc_file
character*38  :: LEC_file
character*58  :: bin_file
character*3   :: year
character*62  :: filename1
character*42  :: filename2

! internal model variables 
integer                                         :: rec_length,i,j,k
integer,          dimension(:,:),   allocatable :: kmT

double precision                                :: volume 
double precision, dimension(:),     allocatable :: dz,area
double precision, dimension(:,:),   allocatable ::                             &
  HTN, HTE, DXT, DYT, TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK, WORK2, WORK3  

real,             dimension(:),     allocatable :: tdepth
real,             dimension(:,:,:), allocatable ::                             &
  DZT, DZU, cPKm, cPKe, VVEL, WVEL
real,             dimension(:,:),   allocatable ::                             &
  cPKm_vint, cPKe_vint, cPKm_avg, cPKe_avg

! Parameters
double precision, parameter ::                                                 &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,                          &
S0ref = 0.035, rho0 = 4.1/3.996*1000, c = 3996,                                &
g = 9.806, salinity_factor = -34.7*1.0E-04, RE = 6.37E06

write (*,*) ''
write (*,*) '--- EXTRA ANALYSIS ---'
write (*,*) ''

!===============================================================================
!  INPUT
!===============================================================================

read  (*,*)        imt, jmt, km
read  (*,'(a120)') grid_file
read  (*,'(a120)') kmt_file
read  (*,'(a120)') in_depths
read  (*,'(a120)') pbc_file
read  (*,'(a58)')  bin_file
read  (*,'(a38)')  LEC_file
!read  (*,'(a3)')   year
read  (*,*)        ntavg
read  (*,*)        nrec_rPm,  nrec_rPe,  nrec_rKm,  nrec_rKe,                  &
                   nrec_gPm,  nrec_gPe,  nrec_gKm,  nrec_gKe,                  &
                   nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                 &
                   nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                 &
                   nrec_TEMP, nrec_VVEL, nrec_WVEL

!===============================================================================
!  GEOMETRY 
!===============================================================================


allocate(                                                                      &
HTN(imt,jmt), HTE(imt,jmt), HUW(imt,jmt), HUS(imt,jmt), WORK(imt,jmt),         &
WORK2(imt,jmt), WORK3(imt,jmt), DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt),    &
DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), dz(km), tdepth(km), kmT(imt,jmt),  &
DZBC(imt,jmt), DZT(imt,jmt,km), DZU(imt,jmt,km) )

inquire (iolength = rec_length) HTN
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

! kmT
inquire (iolength=rec_length) kmT
open    (1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,   &
           status='unknown')
  read  (1,rec=1) kmT
close   (1)

! dz
open    (1,file=in_depths,status='old')
  do k = 1, km
    read  (1,*) dz(k) ! [cm]
  enddo !k
close   (1)

! partial bottom cell depths
inquire (iolength=rec_length) DZBC
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
      DZU(i,j,k) = min( DZT(i    ,jmt,k), DZT(i+1    ,jmt,k),                  &
                        DZT(imt-i,jmt,k), DZT(imt-i+1,jmt,k) )   
    elseif ( i.gt.imt/2 ) then  ! second half, after wrap
      DZU(i,jmt,k) = DZU(imt-i,jmt,k)
    endif
  enddo !i
enddo !k

!  create tdepth 1D array containing depths of TTT-points in [cm]
tdepth(1) = dz(1)/2.0
do k = 2,km
  tdepth(k) = tdepth(k-1)+p5*(dz(k-1)+dz(k))
enddo !k

deallocate( WORK, WORK2, WORK3, DZBC )

!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'

allocate( area(km) )

volume  = sum(TAREA*sum(DZT,3))*1.0E-06 ! volume between bottom and z=0 [m^3]
!write (*,"(A10, ES13.4E2, A4)") 'volume  =',  volume, 'm^3'

do k = 1, km
  area(k) = sum(TAREA,DZT(:,:,k).ne.0)*1.0E-04 ! [m^2]
!  write (*,*) area(k)
enddo !k

! years
if ( ntavg==1 ) then
  start_year = 276
  end_year   = 326
elseif ( ntavg==5 ) then
  start_year = 278
  end_year   = 322
elseif ( ntavg==11 ) then
  start_year = 281
  end_year   = 315 ! check this one again
endif
nt = end_year-start_year+1


allocate( cPKm(imt,jmt,km),   cPKe(imt,jmt,km),                                &
          cPKm_avg(imt,jmt),  cPKe_avg(imt,jmt),                               &
          cPKm_vint(imt,jmt), cPKe_vint(imt,jmt),                              &
          VVEL(imt,jmt,km) )

! naming the output files
if ( ntavg==1 ) then
  open(3,file='/projects/0/samoc/jan/Andree/cPKm_1', &
       form='unformatted',access='direct',recl=imt*jmt)
  open(4,file='/projects/0/samoc/jan/Andree/cPKe_1', &
       form='unformatted',access='direct',recl=imt*jmt)
else if ( ntavg==5 ) then
  open(3,file='/projects/0/samoc/jan/Andree/cPKm_5', &
      form='unformatted',access='direct',recl=imt*jmt)
  open(4,file='/projects/0/samoc/jan/Andree/cPKe_5', &
      form='unformatted',access='direct',recl=imt*jmt)
else if ( ntavg==11 ) then
  open(3,file='/projects/0/samoc/jan/Andree/cPKm_11',&
       form='unformatted',access='direct',recl=imt*jmt)
  open(4,file='/projects/0/samoc/jan/Andree/cPKe_11',&
       form='unformatted',access='direct',recl=imt*jmt)
endif

! yearly vertical integrals
cPKm_avg(:,:) = 0.0
cPKe_avg(:,:) = 0.0

do y = start_year, end_year

  write(year,"(I3)") y
  filename1 = bin_file//year
  filename2 = LEC_file//'_'//year
  write (*,*) filename1
  write (*,*) filename2

  ! open file
  open (1,file=filename1,access='direct',form='unformatted',recl=imt*jmt,      &
        status='unknown')
  open (2,file=filename2,access='direct',form='unformatted',recl=imt*jmt,      &
        status='unknown')

  !=============================================================================
  !  1. cPKm/cPKe
  !=============================================================================
  ! load fields
!  call load_3D_field(imt,jmt,km,2,nrec_cPKm,cPKm)
!  call load_3D_field(imt,jmt,km,2,nrec_cPKe,cPKe)

  !vertical integral
!  call vert_int(cPKm,DZT,cPKm_vint)
!  call vert_int(cPKe,DZT,cPKe_vint)

!  cPKm_avg(:,:) = cPKm_avg(:,:) + 1.0/nt*cPKm_vint(:,:)
!  cPKe_avg(:,:) = cPKe_avg(:,:) + 1.0/nt*cPKe_vint(:,:)

  ! write into interal fields cPKm_vint
!  write (3,rec=y-start_year+1) cPKm_vint(:,:)
!  write (4,rec=y-start_year+1) cPKe_vint(:,:)

  !=============================================================================
  !  2. Eulerian Oversturning Stream Function
  !=============================================================================

  ! load fields
!  call load_3D_field(imt,jmt,km,1,nrec_VVEL,VVEL)




  close(1)
  close(2)

enddo ! y

! averages
!write (3,rec=y-start_year+1) cPKm_avg(:,:)
!write (4,rec=y-start_year+1) cPKe_avg(:,:)










!===============================================================================
!  OUTPUT
!===============================================================================

!write (*,*) ''
!write (*,*) 'OUTPUT'

!===============================================================================
!  SUBROUTINES
!===============================================================================

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine osf(imt,jmt,km,DXU,DZU,VVEL,WVEL,FIELD1,FIELD2)
!
!  overturning stream function calculated both via VVEL and WVEL
!
implicit none

integer, intent(in)                        :: imt, jmt, km
integer                                    :: l
real,    dimension(imt,jmt,km)             :: DZU, VVEL, WVEL
real,    dimension(imt,jmt),   intent(in)  :: DXU
real,    dimension(jmt,km),    intent(out) :: FIELD1,FIELD2


do k = 1,km
  l = km-k+1
  if     ( k==1 ) then
    FIELD1(:,l) = sum(VVEL(:,:,l)*DXU(:,:)/1.0E-04,1)
  elseif ( k>1 ) then
    FIELD1(:,l) = FIELD1(:,l-1) + sum(VVEL(:,:,l)*DXU(:,:)/1.0E-04,1)
  endif
enddo

do j = 1,jmt
  if     ( j==1 ) then
    FIELD2(j,:) = sum(WVEL(:,j,:)*DZU(:,j,:)/1.0E-04,1)
  elseif ( j>1 ) then
    FIELD2(j,:) = FIELD1(j-1,:) + sum(WVEL(:,:,l)*DXU(:,:)/1.0E-04,1)
  endif
enddo



end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine load_2D_field(imt,jmt,nfile,nrec,FIELD)
!
!  loads 2D field (starting with record number nrec) from file with number nfile
!
implicit none

integer, intent(in)                      :: imt, jmt, nfile,  nrec
real,    dimension(imt,jmt), intent(out) :: FIELD

read (nfile,rec=nrec) FIELD

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

INTEGRAL = 0.0
c = 0.0
do j = jmin,jmax
 do i = imin,imax
  if ( MASK(i,j).ne.0.0 ) then
   y = TAREA(i,j) * FIELD(i,j) - c
   t = INTEGRAL + y
   c = ( t - INTEGRAL ) - y
   INTEGRAL = t
  endif
 enddo !i
enddo !j

INTEGRAL = INTEGRAL/1.0E04

!INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)*FIELD(imin:imax,jmin:jmax),          &
!               MASK(imin:imax,jmin:jmax).ne.0.0)/1.0E04
write (*,*) '2D: ', INTEGRAL
end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
subroutine vol_int(imin,jmin,kmin,imax,jmax,kmax,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field, using Kahan Summation
!
implicit none

! input/output variables
integer,             intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,                dimension(imt,jmt,km), intent(in)  :: DZT, FIELD
double precision,    dimension(imt,jmt   ), intent(in)  :: TAREA
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
! write (*,*) '3D: ', INTEGRAL

!INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)                                      &
!         * sum(DZT(imin:imax,jmin:jmax,kmin:kmax)                            &
!         * FIELD(imin:imax,jmin:jmax,kmin:kmax),3))/1.0E06
write (*,*) '3D: ', INTEGRAL

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
subroutine wtt2ttt(imt,jmt,km,W_FIELD,DZT,T_FIELD)
!
!     interpolated real field from WTT-grid to TTT-grid
!     interpolating between two W-levels
!     T-point is located exactly in the middle of two W-points
!     T_FIELD(i,j,k) = W_FIELD(i,j,k) + ( W_FIELD(i,j,k+1) - W_FIELD(i,j,k) ) &
!                                       / dz(k) * dz(k) / 2
!
implicit none

! input/output variables
integer,                         intent(in)  :: imt,jmt,km
real,     dimension(imt,jmt,km), intent(in)  :: DZT
real,     dimension(imt,jmt,km), intent(in)  :: W_FIELD ! WTT
real,     dimension(imt,jmt,km), intent(out) :: T_FIELD ! TTT
! internal variables
double precision,  parameter                 :: p5=0.5, c0=0.0, c3=3.0

do k = 1,km-1
 do j = 1,jmt
  do i = 1,imt
   if ( DZT(i,j,k).ne.c0 ) then
    T_FIELD(i,j,k) = p5 * ( W_FIELD(i,j,k) + W_FIELD(i,j,k+1) )
   endif
  enddo !i
 enddo !j
enddo !k

! bottom level
do j = 1,jmt
  do i = 1,imt
    if ( DZT(i,j,km).ne.c0 ) then
      ! extrapolating to last level with gradient between last two W-levels
      ! T_FIELD(i,j,km) = W_FIELD(i,j,km) 
      !                 + (W_FIELD(i,j,km)-W_FIELD(i,j,km-1)) / dz(k) * dz(k)/2
      T_FIELD(i,j,km) = p5 * (c3*W_FIELD(i,j,km) - W_FIELD(i,j,km-1) )
    endif
  enddo !i
enddo !j

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine uu2tt(imt,jmt,k,DZT_k,DZU_k,TAREA,UAREA,FIELD_k,NEW_FIELD_k)
!
!     interpolates real 2D field from UU-grid to TT-grid on T-depth levels
!
implicit none

! input/output variables
integer,                              intent(in)  :: imt,jmt,k
real,             dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
double precision, dimension(imt,jmt), intent(in)  :: TAREA,UAREA
real,             dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
real,             dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT
real,             dimension(imt,jmt)              :: NEW_FIELD_k_old ! (T)TT
real                                              :: A, B 
! surface integral, used for comparing old and new method od interpolating 

! southernmost gridpoints (j=1) on TTT-grid must be 0, as there would not be 
! UU gridpoints south of it otherwise

do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    NEW_FIELD_k_old(1,j) = 0.25 * ( FIELD_k(1  ,j  )                       &
                                  + FIELD_k(imt,j  )                       &
                                  + FIELD_k(1  ,j-1)                       &
                                  + FIELD_k(imt,j-1) )
    NEW_FIELD_k(1,j) = 0.25                                                & 
                  * ( FIELD_k(1  ,j  ) * UAREA(1  ,j  ) * DZU_k(1  ,j  )   &
                    + FIELD_k(imt,j  ) * UAREA(imt,j  ) * DZU_k(imt,j  )   &
                    + FIELD_k(1  ,j-1) * UAREA(1  ,j-1) * DZU_k(1  ,j-1)   &
                    + FIELD_k(imt,j-1) * UAREA(imt,j-1) * DZU_k(imt,j-1) ) &
                  / TAREA(1,j) / DZT_k(1,j)
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      NEW_FIELD_k_old(i,j) = 0.25 * ( FIELD_k(i  ,j  )                      &
                                    + FIELD_k(i-1,j  )                      &
                                    + FIELD_k(i  ,j-1)                      &
                                    + FIELD_k(i-1,j-1) )
      NEW_FIELD_k(i,j) = 0.25                                               &
                   * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  ) * DZU_k(i  ,j  )   &
                     + FIELD_k(i-1,j  ) * UAREA(i-1,j  ) * DZU_k(i-1,j  )   &
                     + FIELD_k(i  ,j-1) * UAREA(i  ,j-1) * DZU_k(i  ,j-1)   &
                     + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) * DZU_k(i-1,j-1) ) &
                   / TAREA(i,j) / DZT_k(i,j)
    endif
  enddo !i
enddo !j

call surf_int(1,1,imt,jmt,NEW_FIELD_k_old,TAREA,DZT_k,A)
call surf_int(1,1,imt,jmt,NEW_FIELD_k    ,TAREA,DZT_k,B)
write (*,*) 'uu2tt method error: ', k, (A-B)/B

end subroutine uu2tt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine uu2tt_scalar(imt,jmt,k,DZT_k,DZU_k,TAREA,UAREA,FIELD_k,NEW_FIELD_k)
!
!     interpolates real scalar 2D field (as opposed to velocity fields which  
!     can be related to volume transports) from UU-grid to TT-grid on T-depth 
!     levels
!
implicit none

! input/output variables
integer,                              intent(in)  :: imt,jmt,k
real,             dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
double precision, dimension(imt,jmt), intent(in)  :: TAREA,UAREA
real,             dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
real,             dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT
real :: A, B ! surface integral, used for comparing old/new interp. method

! southernmost gridpoints (j=1) on TTT-grid must be 0, as there would not
! be UU gridpoints south of it otherwise
do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    NEW_FIELD_k(1,j) = 0.25 * ( FIELD_k(1  ,j  ) * UAREA(1  ,j  )          &
                              + FIELD_k(imt,j  ) * UAREA(imt,j  )          &
                              + FIELD_k(1  ,j-1) * UAREA(1  ,j-1)          &
                              + FIELD_k(imt,j-1) * UAREA(imt,j-1) )        &
                            / TAREA(1,j)
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      NEW_FIELD_k(i,j) = 0.25 * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  )         &
                                + FIELD_k(i-1,j  ) * UAREA(i-1,j  )         &
                                + FIELD_k(i  ,j-1) * UAREA(i  ,j-1)         &
                                + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) )       &
                              / TAREA(i,j)
    endif
  enddo !i
enddo !j

! results in very small correction O(-4)

end subroutine uu2tt_scalar

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end program
