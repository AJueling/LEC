program analyze_LEC
implicit none

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  this script calculates the integrals of the LEC output file created with
!  LEC.f90
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
  nrec_rPm, nrec_rPe, nrec_rKm, nrec_rKe,                                      &
  nrec_gPm, nrec_gPe, nrec_gKm, nrec_gKe,                                      &
  nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                                  &
  nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                                  &
  nrec_TEMP

character*120 :: grid_file,kmt_file,in_depths,pbc_file,LEC_file
character*3   :: year


! internal model variables 
integer                                         :: rec_length,i,j,k
integer,          dimension(:,:),   allocatable :: kmT

double precision                                :: volume 
double precision, dimension(:),     allocatable :: dz,area
double precision, dimension(:,:),   allocatable ::                             &
  HTN,HTE,DXT,DYT,TAREA,DXU,DYU,UAREA,DZBC,HUW,HUS,WORK,WORK2,WORK3  


real                                            ::                             &
! global integrals
  rPm_int,rPe_int,rKm_int,rKe_int,                                             &
  gPm_int,gPe_int,gKm_int,gKe_int,                                             &
  gPmh_int,gPeh_int,gPms_int,gPes_int,                                         &
  cPem_int,cKem_int,cPKm_int,cPKe_int,                                         &
! Southern Ocean integrals [90S,30S]
! full depth
  rPm_int_SO30,rPe_int_SO30,rKm_int_SO30,rKe_int_SO30,                         &
  gPm_int_SO30,gPe_int_SO30,gKm_int_SO30,gKe_int_SO30,                         &
  gPmh_int_SO30,gPeh_int_SO30,gPms_int_SO30,gPes_int_SO30,                     &
  cPem_int_SO30,cKem_int_SO30,cPKm_int_SO30,cPKe_int_SO30,                     &
! top [0m,300m]
  rPm_int_SO30_top,rPe_int_SO30_top,rKm_int_SO30_top,rKe_int_SO30_top,         &
  cPem_int_SO30_top,cKem_int_SO30_top,cPKm_int_SO30_top,cPKe_int_SO30_top,     &
! bottom [4000m,6000m]
  rPm_int_SO30_bot,rPe_int_SO30_bot,rKm_int_SO30_bot,rKe_int_SO30_bot,         &
  cPem_int_SO30_bot,cKem_int_SO30_bot,cPKm_int_SO30_bot,cPKe_int_SO30_bot,     &
! Southern Ocean Atlantic sector integrals [90S,30S]x[60W,60E]
! full depth
  rPm_int_SA30,rPe_int_SA30,rKm_int_SA30,rKe_int_SA30,                         &
  gPm_int_SA30,gPe_int_SA30,gKm_int_SA30,gKe_int_SA30,                         &
  gPmh_int_SA30,gPeh_int_SA30,gPms_int_SA30,gPes_int_SA30,                     &
  cPem_int_SA30,cKem_int_SA30,cPKm_int_SA30,cPKe_int_SA30,                     &
! top [0m,300m]
  rPm_int_SA30_top,rPe_int_SA30_top,rKm_int_SA30_top,rKe_int_SA30_top,         &
  cPem_int_SA30_top,cKem_int_SA30_top,cPKm_int_SA30_top,cPKe_int_SA30_top,     &
! bottom [4000m,6000m]
  rPm_int_SA30_bot,rPe_int_SA30_bot,rKm_int_SA30_bot,rKe_int_SA30_bot,         &
  cPem_int_SA30_bot,cKem_int_SA30_bot,cPKm_int_SA30_bot,cPKe_int_SA30_bot,     &
! Southern Ocean integrals [90S,45S]
! full depth
  rPm_int_SO45,rPe_int_SO45,rKm_int_SO45,rKe_int_SO45,                         &
  gPm_int_SO45,gPe_int_SO45,gKm_int_SO45,gKe_int_SO45,                         &
  gPmh_int_SO45,gPeh_int_SO45,gPms_int_SO45,gPes_int_SO45,                     &
  cPem_int_SO45,cKem_int_SO45,cPKm_int_SO45,cPKe_int_SO45,                     &
! top [0m,300m]
  rPm_int_SO45_top,rPe_int_SO45_top,rKm_int_SO45_top,rKe_int_SO45_top,         &
  cPem_int_SO45_top,cKem_int_SO45_top,cPKm_int_SO45_top,cPKe_int_SO45_top,     &
! bottom [4000m,6000m]
  rPm_int_SO45_bot,rPe_int_SO45_bot,rKm_int_SO45_bot,rKe_int_SO45_bot,         &
  cPem_int_SO45_bot,cKem_int_SO45_bot,cPKm_int_SO45_bot,cPKe_int_SO45_bot,     &
! Southern Ocean integrals [90S,55S]
! full depth
  rPm_int_SO55,rPe_int_SO55,rKm_int_SO55,rKe_int_SO55,                         &
  gPm_int_SO55,gPe_int_SO55,gKm_int_SO55,gKe_int_SO55,                         &
  gPmh_int_SO55,gPeh_int_SO55,gPms_int_SO55,gPes_int_SO55,                     &
  cPem_int_SO55,cKem_int_SO55,cPKm_int_SO55,cPKe_int_SO55,                     &
! top [0m,300m]
  rPm_int_SO55_top,rPe_int_SO55_top,rKm_int_SO55_top,rKe_int_SO55_top,         &
  cPem_int_SO55_top,cKem_int_SO55_top,cPKm_int_SO55_top,cPKe_int_SO55_top,     &
! bottom [4000m,6000m]
  rPm_int_SO55_bot,rPe_int_SO55_bot,rKm_int_SO55_bot,rKe_int_SO55_bot,         &
  cPem_int_SO55_bot,cKem_int_SO55_bot,cPKm_int_SO55_bot,cPKe_int_SO55_bot,     &
! OHC
  ohc_int, ohc_int_SO30, ohc_int_SA30, ohc_int_SO45, ohc_int_SO55,             &
  ohc_int_SO30_top, ohc_int_SA30_top, ohc_int_SO45_top, ohc_int_SO55_top,      &
  ohc_int_SO30_bot, ohc_int_SA30_bot, ohc_int_SO45_bot, ohc_int_SO55_bot

real,             dimension(:),     allocatable ::                             &
  tdepth
real,             dimension(:,:),   allocatable ::                             &
  gPm,gPe,gKm,gKe,gPmh,gPeh,gPms,gPes  
real,             dimension(:,:,:), allocatable ::                             &
  rPm,rPe,rKm,rKe,cPem,cKem,cPKm,cPKe,TEMP,DZT,DZU

! Parameters
double precision, parameter ::                                                 &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,                          &
S0ref = 0.035, rho0 = 4.1/3.996*1000, c = 3996,                                &
g = 9.806, salinity_factor = -34.7*1.0E-04, RE = 6.37E06

write (*,*) ''
write (*,*) '--- ANALYZING ENERGY BUDGET ---'
write (*,*) ''

!===============================================================================
!  INPUT
!===============================================================================

read  (*,*)        imt, jmt, km
read  (*,'(a120)') grid_file
read  (*,'(a120)') kmt_file
read  (*,'(a120)') in_depths
read  (*,'(a120)') pbc_file
read  (*,'(a120)') LEC_file
read  (*,'(a3)')   year
read  (*,*)        ntavg
read  (*,*)        nrec_rPm,  nrec_rPe,  nrec_rKm,  nrec_rKe,                  &
                   nrec_gPm,  nrec_gPe,  nrec_gKm,  nrec_gKe,                  &
                   nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                 &
                   nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                 &
                   nrec_TEMP

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
!  LOAD LEC FIELDS
!===============================================================================

allocate(                                                                      &
rPm(imt,jmt,km),  rPe(imt,jmt,km),  rKm(imt,jmt,km),  rKe(imt,jmt,km),         &
gPm(imt,jmt),     gPe(imt,jmt),     gKm(imt,jmt),     gKe(imt,jmt),            &
gPmh(imt,jmt),    gPms(imt,jmt),    gPeh(imt,jmt),    gPes(imt,jmt),           &
cPem(imt,jmt,km), cKem(imt,jmt,km), cPKm(imt,jmt,km), cPKe(imt,jmt,km),        & 
TEMP(imt,jmt,km) )


!open file
open (1,file=LEC_file,access='direct',form='unformatted',recl=imt*jmt,  &
        status='unknown')
write (*,*) 'file: ', LEC_file

call load_3D_field(imt,jmt,km,1,nrec_rPm ,rPm )
call load_3D_field(imt,jmt,km,1,nrec_rPe ,rPe )
call load_3D_field(imt,jmt,km,1,nrec_rKm ,rKm )
call load_3D_field(imt,jmt,km,1,nrec_rKe ,rKe )

write(*,*) 'after 3d loading'

call load_2D_field(imt,jmt,   1,nrec_gPm ,gPm )
call load_2D_field(imt,jmt,   1,nrec_gPe ,gPe )
call load_2D_field(imt,jmt,   1,nrec_gKm ,gKm )
call load_2D_field(imt,jmt,   1,nrec_gKe ,gKe )
call load_2D_field(imt,jmt,   1,nrec_gPmh,gPmh)
call load_2D_field(imt,jmt,   1,nrec_gPms,gPms)
call load_2D_field(imt,jmt,   1,nrec_gPeh,gPeh)
call load_2D_field(imt,jmt,   1,nrec_gPes,gPes)
write(*,*) 'after 2d loading'

call load_3D_field(imt,jmt,km,1,nrec_cPem,cPem)
call load_3D_field(imt,jmt,km,1,nrec_cKem,cKem)
call load_3D_field(imt,jmt,km,1,nrec_cPKm,cPKm)
call load_3D_field(imt,jmt,km,1,nrec_cPKe,cPKe)

call load_3D_field(imt,jmt,km,1,nrec_TEMP,TEMP)

close(1)

!===============================================================================
!  CALCULATIONS
!===============================================================================
write (*,*) ''
write (*,*) 'CALCULATIONS'

allocate( area(km) )

volume  = sum(TAREA*sum(DZT,3))*1.0E-06 ! volume between bottom and z=0 [m^3]
write (*,"(A10, ES13.4E2, A4)") 'volume  =',  volume, 'm^3'

do k = 1, km
  area(k) = sum(TAREA,DZT(:,:,k).ne.0)*1.0E-04 ! [m^2]
  write (*,*) area(k)
enddo !k

! global integrals
call vol_int(1,1,1,imt,jmt,km,rPm,TAREA,DZT,rPm_int)
call vol_int(1,1,1,imt,jmt,km,rPe,TAREA,DZT,rPe_int)
call vol_int(1,1,1,imt,jmt,km,rKm,TAREA,DZT,rKm_int)
call vol_int(1,1,1,imt,jmt,km,rKe,TAREA,DZT,rKe_int)

call surf_int(1,1,imt,jmt,gPm,TAREA,DZT(:,:,1),gPm_int)
call surf_int(1,1,imt,jmt,gPe,TAREA,DZT(:,:,1),gPe_int)
call surf_int(1,1,imt,jmt,gKm,TAREA,DZT(:,:,1),gKm_int)
call surf_int(1,1,imt,jmt,gKe,TAREA,DZT(:,:,1),gKe_int)

call surf_int(1,1,imt,jmt,gPmh,TAREA,DZT(:,:,1),gPmh_int)
call surf_int(1,1,imt,jmt,gPms,TAREA,DZT(:,:,1),gPms_int)
call surf_int(1,1,imt,jmt,gPeh,TAREA,DZT(:,:,1),gPeh_int)
call surf_int(1,1,imt,jmt,gPes,TAREA,DZT(:,:,1),gPes_int)

call vol_int(1,1,1,imt,jmt,km,cPem,TAREA,DZT,cPem_int)
call vol_int(1,1,1,imt,jmt,km,cKem,TAREA,DZT,cKem_int)
call vol_int(1,1,1,imt,jmt,km,cPKm,TAREA,DZT,cPKm_int)
call vol_int(1,1,1,imt,jmt,km,cPKe,TAREA,DZT,cPKe_int)


! Southern Ocean integrals [90S,30S]
! full depth
call vol_int(1,1,1,imt,866,km,rPm,TAREA,DZT,rPm_int_SO30)
call vol_int(1,1,1,imt,866,km,rPe,TAREA,DZT,rPe_int_SO30)
call vol_int(1,1,1,imt,866,km,rKm,TAREA,DZT,rKm_int_SO30)
call vol_int(1,1,1,imt,866,km,rKe,TAREA,DZT,rKe_int_SO30)

call surf_int(1,1,imt,866,gPm,TAREA,DZT(:,:,1),gPm_int_SO30)
call surf_int(1,1,imt,866,gPe,TAREA,DZT(:,:,1),gPe_int_SO30)
call surf_int(1,1,imt,866,gKm,TAREA,DZT(:,:,1),gKm_int_SO30)
call surf_int(1,1,imt,866,gKe,TAREA,DZT(:,:,1),gKe_int_SO30)

call surf_int(1,1,imt,866,gPmh,TAREA,DZT(:,:,1),gPmh_int_SO30)
call surf_int(1,1,imt,866,gPms,TAREA,DZT(:,:,1),gPms_int_SO30)
call surf_int(1,1,imt,866,gPeh,TAREA,DZT(:,:,1),gPeh_int_SO30)
call surf_int(1,1,imt,866,gPes,TAREA,DZT(:,:,1),gPes_int_SO30)

call vol_int(1,1,1,imt,866,km,cPem,TAREA,DZT,cPem_int_SO30)
call vol_int(1,1,1,imt,866,km,cKem,TAREA,DZT,cKem_int_SO30)
call vol_int(1,1,1,imt,866,km,cPKm,TAREA,DZT,cPKm_int_SO30)
call vol_int(1,1,1,imt,866,km,cPKe,TAREA,DZT,cPKe_int_SO30)

! top
call vol_int(1,1,1,imt,866,17,rPm,TAREA,DZT,rPm_int_SO30_top)
call vol_int(1,1,1,imt,866,17,rPe,TAREA,DZT,rPe_int_SO30_top)
call vol_int(1,1,1,imt,866,17,rKm,TAREA,DZT,rKm_int_SO30_top)
call vol_int(1,1,1,imt,866,17,rKe,TAREA,DZT,rKe_int_SO30_top)

call vol_int(1,1,1,imt,866,17,cPem,TAREA,DZT,cPem_int_SO30_top)
call vol_int(1,1,1,imt,866,17,cKem,TAREA,DZT,cKem_int_SO30_top)
call vol_int(1,1,1,imt,866,17,cPKm,TAREA,DZT,cPKm_int_SO30_top)
call vol_int(1,1,1,imt,866,17,cPKe,TAREA,DZT,cPKe_int_SO30_top)

! bottom
call vol_int(1,1,35,imt,866,km,rPm,TAREA,DZT,rPm_int_SO30_bot)
call vol_int(1,1,35,imt,866,km,rPe,TAREA,DZT,rPe_int_SO30_bot)
call vol_int(1,1,35,imt,866,km,rKm,TAREA,DZT,rKm_int_SO30_bot)
call vol_int(1,1,35,imt,866,km,rKe,TAREA,DZT,rKe_int_SO30_bot)

call vol_int(1,1,35,imt,866,km,cPem,TAREA,DZT,cPem_int_SO30_bot)
call vol_int(1,1,35,imt,866,km,cKem,TAREA,DZT,cKem_int_SO30_bot)
call vol_int(1,1,35,imt,866,km,cPKm,TAREA,DZT,cPKm_int_SO30_bot)
call vol_int(1,1,35,imt,866,km,cPKe,TAREA,DZT,cPKe_int_SO30_bot)

! Southern Ocean integrals [90S,30S]x[60W,60E]
! full depth
call vol_int(500,1,1,1700,866,km,rPm,TAREA,DZT,rPm_int_SA30)
call vol_int(500,1,1,1700,866,km,rPe,TAREA,DZT,rPe_int_SA30)
call vol_int(500,1,1,1700,866,km,rKm,TAREA,DZT,rKm_int_SA30)
call vol_int(500,1,1,1700,866,km,rKe,TAREA,DZT,rKe_int_SA30)

call surf_int(500,1,1700,866,gPm,TAREA,DZT(:,:,1),gPm_int_SA30)
call surf_int(500,1,1700,866,gPe,TAREA,DZT(:,:,1),gPe_int_SA30)
call surf_int(500,1,1700,866,gKm,TAREA,DZT(:,:,1),gKm_int_SA30)
call surf_int(500,1,1700,866,gKe,TAREA,DZT(:,:,1),gKe_int_SA30)

call surf_int(500,1,1700,866,gPmh,TAREA,DZT(:,:,1),gPmh_int_SA30)
call surf_int(500,1,1700,866,gPms,TAREA,DZT(:,:,1),gPms_int_SA30)
call surf_int(500,1,1700,866,gPeh,TAREA,DZT(:,:,1),gPeh_int_SA30)
call surf_int(500,1,1700,866,gPes,TAREA,DZT(:,:,1),gPes_int_SA30)

call vol_int(500,1,1,1700,866,km,cPem,TAREA,DZT,cPem_int_SA30)
call vol_int(500,1,1,1700,866,km,cKem,TAREA,DZT,cKem_int_SA30)
call vol_int(500,1,1,1700,866,km,cPKm,TAREA,DZT,cPKm_int_SA30)
call vol_int(500,1,1,1700,866,km,cPKe,TAREA,DZT,cPKe_int_SA30)

! top
call vol_int(500,1,1,1700,866,17,rPm,TAREA,DZT,rPm_int_SA30_top)
call vol_int(500,1,1,1700,866,17,rPe,TAREA,DZT,rPe_int_SA30_top)
call vol_int(500,1,1,1700,866,17,rKm,TAREA,DZT,rKm_int_SA30_top)
call vol_int(500,1,1,1700,866,17,rKe,TAREA,DZT,rKe_int_SA30_top)

call vol_int(500,1,1,1700,866,17,cPem,TAREA,DZT,cPem_int_SA30_top)
call vol_int(500,1,1,1700,866,17,cKem,TAREA,DZT,cKem_int_SA30_top)
call vol_int(500,1,1,1700,866,17,cPKm,TAREA,DZT,cPKm_int_SA30_top)
call vol_int(500,1,1,1700,866,17,cPKe,TAREA,DZT,cPKe_int_SA30_top)

! bottom
call vol_int(500,1,35,1700,866,km,rPm,TAREA,DZT,rPm_int_SA30_bot)
call vol_int(500,1,35,1700,866,km,rPe,TAREA,DZT,rPe_int_SA30_bot)
call vol_int(500,1,35,1700,866,km,rKm,TAREA,DZT,rKm_int_SA30_bot)
call vol_int(500,1,35,1700,866,km,rKe,TAREA,DZT,rKe_int_SA30_bot)

call vol_int(500,1,35,1700,866,km,cPem,TAREA,DZT,cPem_int_SA30_bot)
call vol_int(500,1,35,1700,866,km,cKem,TAREA,DZT,cKem_int_SA30_bot)
call vol_int(500,1,35,1700,866,km,cPKm,TAREA,DZT,cPKm_int_SA30_bot)
call vol_int(500,1,35,1700,866,km,cPKe,TAREA,DZT,cPKe_int_SA30_bot)

!  ADD SO45 and SO55 quantities
! Southern Ocean integrals [90S,45S]
! full depth
call vol_int(1,1,1,imt,677,km,rPm,TAREA,DZT,rPm_int_SO45)
call vol_int(1,1,1,imt,677,km,rPe,TAREA,DZT,rPe_int_SO45)
call vol_int(1,1,1,imt,677,km,rKm,TAREA,DZT,rKm_int_SO45)
call vol_int(1,1,1,imt,677,km,rKe,TAREA,DZT,rKe_int_SO45)

call surf_int(1,1,imt,677,gPm,TAREA,DZT(:,:,1),gPm_int_SO45)
call surf_int(1,1,imt,677,gPe,TAREA,DZT(:,:,1),gPe_int_SO45)
call surf_int(1,1,imt,677,gKm,TAREA,DZT(:,:,1),gKm_int_SO45)
call surf_int(1,1,imt,677,gKe,TAREA,DZT(:,:,1),gKe_int_SO45)

call surf_int(1,1,imt,677,gPmh,TAREA,DZT(:,:,1),gPmh_int_SO45)
call surf_int(1,1,imt,677,gPms,TAREA,DZT(:,:,1),gPms_int_SO45)
call surf_int(1,1,imt,677,gPeh,TAREA,DZT(:,:,1),gPeh_int_SO45)
call surf_int(1,1,imt,677,gPes,TAREA,DZT(:,:,1),gPes_int_SO45)

call vol_int(1,1,1,imt,677,km,cPem,TAREA,DZT,cPem_int_SO45)
call vol_int(1,1,1,imt,677,km,cKem,TAREA,DZT,cKem_int_SO45)
call vol_int(1,1,1,imt,677,km,cPKm,TAREA,DZT,cPKm_int_SO45)
call vol_int(1,1,1,imt,677,km,cPKe,TAREA,DZT,cPKe_int_SO45)

! top
call vol_int(1,1,1,imt,677,17,rPm,TAREA,DZT,rPm_int_SO45_top)
call vol_int(1,1,1,imt,677,17,rPe,TAREA,DZT,rPe_int_SO45_top)
call vol_int(1,1,1,imt,677,17,rKm,TAREA,DZT,rKm_int_SO45_top)
call vol_int(1,1,1,imt,677,17,rKe,TAREA,DZT,rKe_int_SO45_top)

call vol_int(1,1,1,imt,677,17,cPem,TAREA,DZT,cPem_int_SO45_top)
call vol_int(1,1,1,imt,677,17,cKem,TAREA,DZT,cKem_int_SO45_top)
call vol_int(1,1,1,imt,677,17,cPKm,TAREA,DZT,cPKm_int_SO45_top)
call vol_int(1,1,1,imt,677,17,cPKe,TAREA,DZT,cPKe_int_SO45_top)

! bottom
call vol_int(1,1,35,imt,677,km,rPm,TAREA,DZT,rPm_int_SO45_bot)
call vol_int(1,1,35,imt,677,km,rPe,TAREA,DZT,rPe_int_SO45_bot)
call vol_int(1,1,35,imt,677,km,rKm,TAREA,DZT,rKm_int_SO45_bot)
call vol_int(1,1,35,imt,677,km,rKe,TAREA,DZT,rKe_int_SO45_bot)

call vol_int(1,1,35,imt,677,km,cPem,TAREA,DZT,cPem_int_SO45_bot)
call vol_int(1,1,35,imt,677,km,cKem,TAREA,DZT,cKem_int_SO45_bot)
call vol_int(1,1,35,imt,677,km,cPKm,TAREA,DZT,cPKm_int_SO45_bot)
call vol_int(1,1,35,imt,677,km,cPKe,TAREA,DZT,cPKe_int_SO45_bot)


! Southern Ocean integrals [90S,55S]
! full depth
call vol_int(1,1,1,imt,518,km,rPm,TAREA,DZT,rPm_int_SO55)
call vol_int(1,1,1,imt,518,km,rPe,TAREA,DZT,rPe_int_SO55)
call vol_int(1,1,1,imt,518,km,rKm,TAREA,DZT,rKm_int_SO55)
call vol_int(1,1,1,imt,518,km,rKe,TAREA,DZT,rKe_int_SO55)

call surf_int(1,1,imt,518,gPm,TAREA,DZT(:,:,1),gPm_int_SO55)
call surf_int(1,1,imt,518,gPe,TAREA,DZT(:,:,1),gPe_int_SO55)
call surf_int(1,1,imt,518,gKm,TAREA,DZT(:,:,1),gKm_int_SO55)
call surf_int(1,1,imt,518,gKe,TAREA,DZT(:,:,1),gKe_int_SO55)

call surf_int(1,1,imt,518,gPmh,TAREA,DZT(:,:,1),gPmh_int_SO55)
call surf_int(1,1,imt,518,gPms,TAREA,DZT(:,:,1),gPms_int_SO55)
call surf_int(1,1,imt,518,gPeh,TAREA,DZT(:,:,1),gPeh_int_SO55)
call surf_int(1,1,imt,518,gPes,TAREA,DZT(:,:,1),gPes_int_SO55)

call vol_int(1,1,1,imt,518,km,cPem,TAREA,DZT,cPem_int_SO55)
call vol_int(1,1,1,imt,518,km,cKem,TAREA,DZT,cKem_int_SO55)
call vol_int(1,1,1,imt,518,km,cPKm,TAREA,DZT,cPKm_int_SO55)
call vol_int(1,1,1,imt,518,km,cPKe,TAREA,DZT,cPKe_int_SO55)

! top
call vol_int(1,1,1,imt,518,17,rPm,TAREA,DZT,rPm_int_SO55_top)
call vol_int(1,1,1,imt,518,17,rPe,TAREA,DZT,rPe_int_SO55_top)
call vol_int(1,1,1,imt,518,17,rKm,TAREA,DZT,rKm_int_SO55_top)
call vol_int(1,1,1,imt,518,17,rKe,TAREA,DZT,rKe_int_SO55_top)

call vol_int(1,1,1,imt,518,17,cPem,TAREA,DZT,cPem_int_SO55_top)
call vol_int(1,1,1,imt,518,17,cKem,TAREA,DZT,cKem_int_SO55_top)
call vol_int(1,1,1,imt,518,17,cPKm,TAREA,DZT,cPKm_int_SO55_top)
call vol_int(1,1,1,imt,518,17,cPKe,TAREA,DZT,cPKe_int_SO55_top)

! bottom
call vol_int(1,1,35,imt,518,km,rPm,TAREA,DZT,rPm_int_SO55_bot)
call vol_int(1,1,35,imt,518,km,rPe,TAREA,DZT,rPe_int_SO55_bot)
call vol_int(1,1,35,imt,518,km,rKm,TAREA,DZT,rKm_int_SO55_bot)
call vol_int(1,1,35,imt,518,km,rKe,TAREA,DZT,rKe_int_SO55_bot)

call vol_int(1,1,35,imt,518,km,cPem,TAREA,DZT,cPem_int_SO55_bot)
call vol_int(1,1,35,imt,518,km,cKem,TAREA,DZT,cKem_int_SO55_bot)
call vol_int(1,1,35,imt,518,km,cPKm,TAREA,DZT,cPKm_int_SO55_bot)
call vol_int(1,1,35,imt,518,km,cPKe,TAREA,DZT,cPKe_int_SO55_bot)


! OHC
call vol_int(1,1,1,imt,jmt,km,TEMP,TAREA,DZT,ohc_int)
call vol_int(1,1,1,imt,866,km,TEMP,TAREA,DZT,ohc_int_SO30)
call vol_int(500,1,1,1700,866,km,TEMP,TAREA,DZT,ohc_int_SA30)
call vol_int(1,1,1,imt,677,km,TEMP,TAREA,DZT,ohc_int_SO45)
call vol_int(1,1,1,imt,518,km,TEMP,TAREA,DZT,ohc_int_SO55)
call vol_int(1,1,1,imt,866,17,TEMP,TAREA,DZT,ohc_int_SO30_top)
call vol_int(500,1,1,1700,866,17,TEMP,TAREA,DZT,ohc_int_SA30_top)
call vol_int(1,1,1,imt,677,17,TEMP,TAREA,DZT,ohc_int_SO45_top)
call vol_int(1,1,1,imt,518,17,TEMP,TAREA,DZT,ohc_int_SO55_top)
call vol_int(1,1,35,imt,866,km,TEMP,TAREA,DZT,ohc_int_SO30_bot)
call vol_int(500,1,35,1700,866,km,TEMP,TAREA,DZT,ohc_int_SA30_bot)
call vol_int(1,1,35,imt,677,km,TEMP,TAREA,DZT,ohc_int_SO45_bot)
call vol_int(1,1,35,imt,518,km,TEMP,TAREA,DZT,ohc_int_SO55_bot)


!===============================================================================
!  OUTPUT
!===============================================================================

write (*,*) ''
write (*,*) 'OUTPUT'


! open  (2,file='output/level_LEC_'//(year//'.out'),form='formatted')
! write (2,200) 'k','depth_k','area_k','rho_ref_k','pd_ref_k','n0_k',            &
!               'rPm_k','rPe_k','rKm_k','rKe_k',                                 &
!               'cPem_k','cKem_k','cPKm_k','cPKe_k'
! do k = 1,km
!   !surface integrals of energy budget terms
!   write (2,201) real(k), tdepth(k), area(k), rho_ref(k), pd_ref(k), n0(k),     &
!                 rPm_sint(k),  rPe_sint(k),  rKm_sint(k),  rKe_sint(k),         &
!                 cPem_sint(k), cKem_sint(k), cPKm_sint(k), cPKe_sint(k)
! enddo
! close (2)

200 FORMAT (156(A,","),A)
201 FORMAT (156(E15.7,","),E15.7)
if ( ntavg==1 ) then
  open  (2,file='output/analysis_LEC_1_'//(year//'.out'),form='formatted')
else if ( ntavg==5 ) then
  open  (2,file='output/analysis_LEC_5_'//(year//'.out'),form='formatted')
else if ( ntavg==11 ) then
  open  (2,file='output/analysis_LEC_11_'//(year//'.out'),form='formatted')
endif
write (2,200)                                                                  &
! global integrals
  'rPm','rPe','rKm','rKe',                                                     &
  'gPm','gPe','gKm','gKe',                                                     &
  'gPmh','gPeh','gPms','gPes',                                                 &
  'cPem','cKem','cPKm','cPKe',                                                 &
! Southern Ocean integrals [90S,30S]
! full depth
  'rPm_SO30','rPe_SO30','rKm_SO30','rKe_SO30',                                 &
  'gPm_SO30','gPe_SO30','gKm_SO30','gKe_SO30',                                 &
  'gPmh_SO30','gPeh_SO30','gPms_SO30','gPes_SO30',                             &
  'cPem_SO30','cKem_SO30','cPKm_SO30','cPKe_SO30',                             &
! top [0m,300m]
  'rPm_SO30_top','rPe_SO30_top','rKm_SO30_top','rKe_SO30_top',                 &
  'cPem_SO30_top','cKem_SO30_top',                                             &
  'cPKm_SO30_top','cPKe_SO30_top',                                             &
! bottom [4000m,6000m]
  'rPm_SO30_bot','rPe_SO30_bot','rKm_SO30_bot','rKe_SO30_bot',                 &
  'cPem_SO30_bot','cKem_SO30_bot',                                             &
  'cPKm_SO30_bot','cPKe_SO30_bot',                                             &
! Southern Ocean integrals [90S,30S]
! full depth
  'rPm_SA30','rPe_SA30','rKm_SA30','rKe_SA30',                                 &
  'gPm_SA30','gPe_SA30','gKm_SA30','gKe_SA30',                                 &
  'gPmh_SA30','gPeh_SA30','gPms_SA30','gPes_SA30',                             &
  'cPem_SA30','cKem_SA30','cPKm_SA30','cPKe_SA30',                             &
! top [0m,300m]
  'rPm_SA30_top','rPe_SA30_top','rKm_SA30_top','rKe_SA30_top',                 &
  'cPem_SA30_top','cKem_SA30_top',                                             &
  'cPKm_SA30_top','cPKe_SA30_top',                                             &
! bottom [4000m,6000m]
  'rPm_SA30_bot','rPe_SA30_bot','rKm_SA30_bot','rKe_SA30_bot',                 &
  'cPem_SA30_bot','cKem_SA30_bot',                                             &
  'cPKm_SA30_bot','cPKe_SA30_bot',                                             &
! Southern Ocean integrals [90S,45S]
! full depth
  'rPm_SO45','rPe_SO45','rKm_SO45','rKe_SO45',                                 &
  'gPm_SO45','gPe_SO45','gKm_SO45','gKe_SO45',                                 &
  'gPmh_SO45','gPeh_SO45','gPms_SO45','gPes_SO45',                             &
  'cPem_SO45','cKem_SO45','cPKm_SO45','cPKe_SO45',                             &
! top [0m,300m]
  'rPm_SO45_top','rPe_SO45_top','rKm_SO45_top','rKe_SO45_top',                 &
  'cPem_SO45_top','cKem_SO45_top',                                             &
  'cPKm_SO45_top','cPKe_SO45_top',                                             &
! bottom [4000m,6000m]
  'rPm_SO45_bot','rPe_SO45_bot','rKm_SO45_bot','rKe_SO45_bot',                 &
  'cPem_SO45_bot','cKem_SO45_bot',                                             &
  'cPKm_SO45_bot','cPKe_SO45_bot',                                             &
! Southern Ocean integrals [90S,55S]
! full depth
  'rPm_SO55','rPe_SO55','rKm_SO55','rKe_SO55',                                 &
  'gPm_SO55','gPe_SO55','gKm_SO55','gKe_SO55',                                 &
  'gPmh_SO55','gPeh_SO55','gPms_SO55','gPes_SO55',                             &
  'cPem_SO55','cKem_SO55','cPKm_SO55','cPKe_SO55',                             &
! top [0m,300m]
  'rPm_SO55_top','rPe_SO55_top','rKm_SO55_top','rKe_SO55_top',                 &
  'cPem_SO55_top','cKem_SO55_top',                                             &
  'cPKm_SO55_top','cPKe_SO55_top',                                             &
! bottom [4000m,6000m]
  'rPm_SO55_bot','rPe_SO55_bot','rKm_SO55_bot','rKe_SO55_bot',                 &
  'cPem_SO55_bot','cKem_SO55_bot',                                             &
  'cPKm_SO55_bot','cPKe_SO55_bot',                                             &
! OHC
  'ohc_int', 'ohc_int_SO30', 'ohc_int_SA30', 'ohc_int_SO45', 'ohc_int_SO55',   &
  'ohc_int_SO30_top', 'ohc_int_SA30_top', 'ohc_int_SO45_top', 'ohc_int_SO55_top', &
  'ohc_int_SO30_bot', 'ohc_int_SA30_top', 'ohc_int_SO45_bot', 'ohc_int_SO55_bot'

write (2,201)                                                                  &
! global integrals
  rPm_int,rPe_int,rKm_int,rKe_int,                                             &
  gPm_int,gPe_int,gKm_int,gKe_int,                                             &
  gPmh_int,gPeh_int,gPms_int,gPes_int,                                         &
  cPem_int,cKem_int,cPKm_int,cPKe_int,                                         &
! Southern Ocean integrals [90S,30S]
! full depth
  rPm_int_SO30,rPe_int_SO30,rKm_int_SO30,rKe_int_SO30,                         &
  gPm_int_SO30,gPe_int_SO30,gKm_int_SO30,gKe_int_SO30,                         &
  gPmh_int_SO30,gPeh_int_SO30,gPms_int_SO30,gPes_int_SO30,                     &
  cPem_int_SO30,cKem_int_SO30,cPKm_int_SO30,cPKe_int_SO30,                     &
! top [0m,300m]
  rPm_int_SO30_top,rPe_int_SO30_top,rKm_int_SO30_top,rKe_int_SO30_top,         &
  cPem_int_SO30_top,cKem_int_SO30_top,cPKm_int_SO30_top,cPKe_int_SO30_top,     &
! bottom [4000m,6000m]
  rPm_int_SO30_bot,rPe_int_SO30_bot,rKm_int_SO30_bot,rKe_int_SO30_bot,         &
  cPem_int_SO30_bot,cKem_int_SO30_bot,cPKm_int_SO30_bot,cPKe_int_SO30_bot,     &
! Southern Ocean integrals [90S,30S]
! full depth
  rPm_int_SA30,rPe_int_SA30,rKm_int_SA30,rKe_int_SA30,                         &
  gPm_int_SA30,gPe_int_SA30,gKm_int_SA30,gKe_int_SA30,                         &
  gPmh_int_SA30,gPeh_int_SA30,gPms_int_SA30,gPes_int_SA30,                     &
  cPem_int_SA30,cKem_int_SA30,cPKm_int_SA30,cPKe_int_SA30,                     &
! top [0m,300m]
  rPm_int_SA30_top,rPe_int_SA30_top,rKm_int_SA30_top,rKe_int_SA30_top,         &
  cPem_int_SA30_top,cKem_int_SA30_top,cPKm_int_SA30_top,cPKe_int_SA30_top,     &
! bottom [4000m,6000m]
  rPm_int_SA30_bot,rPe_int_SA30_bot,rKm_int_SA30_bot,rKe_int_SA30_bot,         &
  cPem_int_SA30_bot,cKem_int_SA30_bot,cPKm_int_SA30_bot,cPKe_int_SA30_bot,     &
! Southern Ocean integrals [90S,45S]
! full depth
  rPm_int_SO45,rPe_int_SO45,rKm_int_SO45,rKe_int_SO45,                         &
  gPm_int_SO45,gPe_int_SO45,gKm_int_SO45,gKe_int_SO45,                         &
  gPmh_int_SO45,gPeh_int_SO45,gPms_int_SO45,gPes_int_SO45,                     &
  cPem_int_SO45,cKem_int_SO45,cPKm_int_SO45,cPKe_int_SO45,                     &
! top [0m,300m]
  rPm_int_SO45_top,rPe_int_SO45_top,rKm_int_SO45_top,rKe_int_SO45_top,         &
  cPem_int_SO45_top,cKem_int_SO45_top,cPKm_int_SO45_top,cPKe_int_SO45_top,     &
! bottom [4000m,6000m]
  rPm_int_SO45_bot,rPe_int_SO45_bot,rKm_int_SO45_bot,rKe_int_SO45_bot,         &
  cPem_int_SO45_bot,cKem_int_SO45_bot,cPKm_int_SO45_bot,cPKe_int_SO45_bot,     &
! Southern Ocean integrals [90S,55S]
! full depth
  rPm_int_SO55,rPe_int_SO55,rKm_int_SO55,rKe_int_SO55,                         &
  gPm_int_SO55,gPe_int_SO55,gKm_int_SO55,gKe_int_SO55,                         &
  gPmh_int_SO55,gPeh_int_SO55,gPms_int_SO55,gPes_int_SO55,                     &
  cPem_int_SO55,cKem_int_SO55,cPKm_int_SO55,cPKe_int_SO55,                     &
! top [0m,300m]
  rPm_int_SO55_top,rPe_int_SO55_top,rKm_int_SO55_top,rKe_int_SO55_top,         &
  cPem_int_SO55_top,cKem_int_SO55_top,cPKm_int_SO55_top,cPKe_int_SO55_top,     &
! bottom [4000m,6000m]
  rPm_int_SO55_bot,rPe_int_SO55_bot,rKm_int_SO55_bot,rKe_int_SO55_bot,         &
  cPem_int_SO55_bot,cKem_int_SO55_bot,cPKm_int_SO55_bot,cPKe_int_SO55_bot,     &
! OHC
  ohc_int, ohc_int_SO30, ohc_int_SA30, ohc_int_SO45, ohc_int_SO55,             &
  ohc_int_SO30_top, ohc_int_SA30_top, ohc_int_SO45_top, ohc_int_SO55_top,      &
  ohc_int_SO30_bot, ohc_int_SA30_bot, ohc_int_SO45_bot, ohc_int_SO55_bot
write (*,*) '>>> done writing output'


!===============================================================================
!===============================================================================
end program
