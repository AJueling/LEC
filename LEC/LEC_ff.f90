program LEC
implicit none

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  this script calculates the Lorenz Energy Cycle terms for POP model output
!  it should be called with a script that will compile this file and supply all
!  the input parameters
!  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!===============================================================================
!  variables
!===============================================================================

! user input 
integer :: imt,jmt,km,ntavg

integer ::                                                                  &
nrec_UVEL, nrec_VVEL, nrec_KE, nrec_RHO, nrec_RHO2, nrec_Q, nrec_ALPHA01,   & 
nrec_SHF, nrec_BETA01, nrec_SFWF, nrec_STFRHO, nrec_SSFRHO, nrec_TAUX,      &
nrec_TAUY, nrec_UTAUX, nrec_VTAUY, nrec_RHOU, nrec_RHOV, nrec_WVEL,         &
nrec_UVEL2, nrec_VVEL2, nrec_UV, nrec_UW, nrec_VW, nrec_RHOW, nrec_PD,      &
nrec_PD2, nrec_PDU, nrec_PDV,                                               &
nrec_EVAPRHO, nrec_PRECIPRHO, nrec_SWRESTRHO, nrec_SSRESTRHO,               &
nrec_RUNOFFRHO, nrec_TEMP, nrec_SALT, nrec_SSH, nrec_SWNET, nrec_LWNET,     &
nrec_LATENT, nrec_SENSIBLE, nrec_TSREST, nrec_EVAP, nrec_PRECIP,            &
nrec_SWREST, nrec_SSREST, nrec_RUNOFF, nrec_UEU, nrec_VNU, nrec_UEV,        &
nrec_VNV, nrec_WTU, nrec_WTV, nrec_UET, nrec_VNT, nrec_WTT, nrec_UES,       &
nrec_VNS, nrec_WTS, nrec_ADVT, nrec_ADVS, nrec_URHO, nrec_VRHO, nrec_WRHO,  &
nrec_PDW, nrec_UPD, nrec_VPD, nrec_WPD
! obsolete: RHO2

character*120 :: grid_file, kmt_file, in_depths, pbc_file, tavg_file
character*3   :: year

! program output
real :: rho_avg, pd_avg,                                                    &
rKm_int, rKe_int, rPm_int, rPe_int,                                         &
gKm_int, gPmh_int, gPms_int, gPm_int,                                       &
gKe_int, gPeh_int, gPes_int, gPe_int,                                       &
cPem_int, cKem_int, cPKm_int, cPKe_int,                                     &
dPm, dPe, dKm, dKe,                                                         &
V2U2, V2U2_vol

double precision :: volume, RHO_new, PD_top, PD_bottom

real, dimension(:), allocatable ::                                          &
rKm_sint,     rKe_sint,     rPm_sint,     rPe_sint,                         &
cPem_sint,    cKem_sint,    cPKm_sint,    cPKe_sint

! internal model variables 
integer :: rec_length, i, j, k

integer, dimension(:,:), allocatable :: kmT

real, dimension(:,:), allocatable ::                                        &
ALPHA01, SHF, BETA01, SFWF, STFRHO, SSFRHO, TAUX, TAUY, UTAUX, VTAUY,       &
RHOW_k, RHOU_k, RHOV_k,    &
gPmh, gPms, gPm, gPeh, gPes, gPe, gKm, gKe, TT_gKe, TT_gKm,    &
EVAP, PRECIP, SWREST, SSREST, RUNOFF, EVAPRHO, PRECIPRHO, SWRESTRHO,        &
SSRESTRHO, RUNOFFRHO, SWNET, LWNET, LATENT, SENSIBLE, TSREST,&
gPt, gKt, TT_gKt, rKt_k,                                                    &
UEU_k, UEV_k, VNU_k, VNV_k, WTU_k, WTV_k, UUE_k, VUN_k, UTE_k, VTN_k,       &
UEU_eddy_k, UEV_eddy_k, VNU_eddy_k, VNV_eddy_k,                             &
TTT_UEU_k, TTT_UEV_k, TTT_VNU_k, TTT_VNV_k, TTT_WTU_k, TTT_WTV_k,           &
TTT_UEU_eddy_k, TTT_UEV_eddy_k, TTT_VNU_eddy_k, TTT_VNV_eddy_k,             &
UPD_k, VPD_k, TTT_UPD_k, TTT_VPD_k,                                         &
UPD_eddy_k, VPD_eddy_k, TTT_UPD_eddy_k, TTT_VPD_eddy_k,                     &
WUU_WTU_eddy_k,  WUU_WTV_eddy_k, WTT_WTU_eddy_k, WTT_WTV_eddy_k


real, dimension(:,:,:), allocatable ::                                      &
UVEL, VVEL, KE, RHO, Q, RHOU, RHOV, WVEL, UVEL2, VVEL2, UV, UW, VW,   &
RHOW, PD, PD2, PDU, PDV, rKm, rKe, rPm, rPe, cPem, cKem, cPKm, cPKe, SALT,  &
TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV, TTT_UW, TTT_VW,           &
DRHODX, DRHODY, DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DPDDX, DPDDY,           &
TTT_WVEL, TTT_RHOW, TTT_PDW,  TTT_rKm, TTT_rKe, TEMP, PDW,                  &
WUU_WVEL, WTT_WTU_eddy, WTT_WTV_eddy, TTT_WTU_eddy, TTT_WTV_eddy

double precision, dimension(:), allocatable ::                              &
dz,z1, z2, area, rho_ref, n0, n0_inv, pd_ref, tdepth, salt_ref, temp_ref,   &
pressz, rho_new_ref, pd_new_ref, n0_new

double precision, dimension(:,:), allocatable ::                            &
HTN, HTE, WORK, DXT, DYT, TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK2,    &
WORK3  

real, dimension(:,:,:), allocatable :: DZT, DZU

! Testing
real :: WVEL_int, RHO_int

! WTT-calculations
real, dimension(:,:,:), allocatable ::                                      &
WTT_RHO, WTT_UVEL, WTT_VVEL, WTT_UVEL2, WTT_VVEL2, WTT_UV, WTT_DUDX,        &
WTT_DUDY, WTT_DUDZ, WTT_DVDX, WTT_DVDY, WTT_DVDZ, cKem_WTT, cPKm_WTT,       &
cPKe_WTT

real :: cKem_WTT_int, cPKm_WTT_int, cPKe_WTT_int

! Parameters
double precision, parameter ::                                              &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,                       &
S0ref = 0.035, rho0 = 4.1/3.996*1000, c = 3996,                             &
g = 9.806, salinity_factor = -34.7*1.0E-04, RE = 6.37E06
! salinity factor from http://www.cesm.ucar.edu/models/cesm1.0/cesm/
! cesmBbrowser/html_code/pop/constants.F90.html "flux (kg/m^2/s) to salt flux
! (msu*cm/s)"


write (*,*) ''
write (*,*) '--- ENERGY BUDGET ---'
write (*,*) ''

!===============================================================================
!  user input
!===============================================================================

read  (*,*)        imt, jmt, km
read  (*,'(a120)') grid_file
read  (*,'(a120)') kmt_file
read  (*,'(a120)') in_depths
read  (*,'(a120)') pbc_file
read  (*,'(a120)') tavg_file
read  (*,'(a3)')   year
read  (*,*)        ntavg
read  (*,*) nrec_UVEL, nrec_VVEL, nrec_UVEL2, nrec_VVEL2, nrec_STFRHO,         &
  nrec_SSFRHO, nrec_EVAPRHO, nrec_PRECIPRHO, nrec_SWRESTRHO, nrec_SSRESTRHO,   &
  nrec_RUNOFFRHO, nrec_UTAUX, nrec_VTAUY, nrec_ALPHA01, nrec_BETA01, nrec_KE,  &
  nrec_TEMP, nrec_SALT, nrec_RHO, nrec_UV, nrec_SSH, nrec_SWNET, nrec_LWNET,   &
  nrec_LATENT, nrec_SENSIBLE, nrec_TSREST, nrec_EVAP, nrec_PRECIP,             &
  nrec_SWREST, nrec_SSREST, nrec_RUNOFF, nrec_SHF, nrec_SFWF,                  &
  nrec_TAUX, nrec_TAUY, nrec_WVEL, nrec_UEU, nrec_VNU, nrec_UEV, nrec_VNV,     &
  nrec_WTU, nrec_WTV, nrec_UET, nrec_VNT, nrec_WTT, nrec_UES, nrec_VNS,        &
  nrec_WTS, nrec_ADVT, nrec_ADVS, nrec_Q, nrec_PD, nrec_PD2, nrec_RHO2,        &
  nrec_UW, nrec_VW, nrec_RHOU, nrec_RHOV, nrec_RHOW, nrec_URHO, nrec_VRHO,     &
  nrec_WRHO, nrec_PDW, nrec_PDU, nrec_PDV, nrec_UPD, nrec_VPD, nrec_WPD,       &
  nrec_UEU, nrec_UEV, nrec_VNU, nrec_VNV, nrec_WTU, nrec_WTV

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

allocate( dz(km), tdepth(km), kmT(imt,jmt), DZBC(imt,jmt) )
allocate( DZT(imt,jmt,km), DZU(imt,jmt,km) )

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
deallocate( DZBC )
  
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

!  create tdepth 1D array containing depths of TTT-points in [cm]
tdepth(1) = dz(1)/2.0
do k = 2,km
  tdepth(k) = tdepth(k-1)+p5*(dz(k-1)+dz(k))
enddo !k

!===============================================================================
!  read fields
!===============================================================================


! 3D
allocate(                                                                      &
  UVEL(imt,jmt,km), VVEL(imt,jmt,km), RHO(imt,jmt,km),                         &
  Q(imt,jmt,km),    PD(imt,jmt,km),   SALT(imt,jmt,km), TEMP(imt,jmt,km) )

!open file
inquire (iolength=rec_length) kmT
!write (*,*) 'UVEL rec_length', rec_length, imt,jmt, imt*jmt
open (1,file=tavg_file,access='direct',form='unformatted',recl=rec_length,     &
        status='unknown')

!read 3-D fields
call load_3D_field(imt,jmt,km,1,nrec_RHO, RHO ) ! [g/cm^3]
call load_3D_field(imt,jmt,km,1,nrec_PD,  PD  ) ! [g/cm^3]
call load_3D_field(imt,jmt,km,1,nrec_Q,   Q   ) ! [g/cm^4]
call load_3D_field(imt,jmt,km,1,nrec_UVEL,UVEL) ! [cm/s]
call load_3D_field(imt,jmt,km,1,nrec_VVEL,VVEL) ! [cm/s]
call load_3D_field(imt,jmt,km,1,nrec_SALT,SALT) ! [g/kg]
call load_3D_field(imt,jmt,km,1,nrec_TEMP,TEMP) ! [degC]

! to show how Q is calculated incorrectly at surface (and bottom)
! write (*,*) 'Q  ', Q(500,500,1)
! write (*,*) 'PD', PD(500,500,1)
! write (*,*) 'PD', PD(500,500,2)
! write (*,*) 'dPD/dz', (PD(500,500,1)-PD(500,500,2))/1000.0
! write (*,*) 'Q  ', Q(600,600,1)
! write (*,*) 'PD', PD(600,600,1)
! write (*,*) 'PD', PD(600,600,2)
! write (*,*) 'dPD/dz', (PD(600,600,1)-PD(600,600,2))/1000.0

write (*,*) 'file: ', tavg_file

!===============================================================================
!  calculate averaged quantities
!===============================================================================

allocate( area(km), rho_ref(km), pd_ref(km), n0(km), n0_inv(km),            &
          salt_ref(km), temp_ref(km), pressz(km), rho_new_ref(km),          &
          pd_new_ref(km), n0_new(km) )

volume  = sum(TAREA*sum(DZT,3))*1.0E-06 ! volume between bottom and z=0 [m^3]
call vol_int(1,1,1,imt,jmt,km,RHO,TAREA,DZT,rho_avg)
rho_avg = rho_avg/volume
call vol_int(1,1,1,imt,jmt,km, PD,TAREA,DZT, pd_avg)
pd_avg  =  pd_avg/volume
write (*,"(A10, ES13.4E2, A7)") 'rho_avg =', rho_avg, 'g/cm^3'
write (*,"(A10, ES13.4E2, A7)") 'pd_avg  =',  pd_avg, 'g/cm^3'
write (*,"(A10, ES13.4E2, A7)") 'rho0    =',    rho0, 'kg/m^3'
write (*,"(A10, ES13.4E2, A4)") 'volume  =',  volume, 'm^3'

do k = 1, km
  ! (T)area per level [m^2]
  area(k)     = sum(TAREA,             DZT(:,:,k).ne.0)*1.0E-04
  ! rho_ref: average in situ density at level k [kg/m^3]
  rho_ref(k)  = sum(TAREA*RHO(:,:,k),  DZT(:,:,k).ne.0)/area(k)*1.0E-01 
  ! pd_ref: average potential density at level k [kg/m^3]
  pd_ref(k)   = sum(TAREA*PD(:,:,k),   DZT(:,:,k).ne.0)/area(k)*1.0E-01
  ! salt_ref: average salinity at level k [g/kg]
  salt_ref(k) = sum(TAREA*SALT(:,:,k), DZT(:,:,k).ne.0)/area(k)*1.0E-04
  ! temp_ref: average potential temperature at level k [degC]
  temp_ref(k) = sum(TAREA*TEMP(:,:,k), DZT(:,:,k).ne.0)/area(k)*1.0E-04
  ! n0: average Q at level k [kg/m^4]
  n0(k)       = sum(TAREA*Q(:,:,k),    DZT(:,:,k).ne.0)/area(k)*1.0E01
  if (k==1)  n0(k) = n0(k)*2.0
  if (k==km) n0(k) = n0(k)*2.0
  ! inverse n0 [m^4/kg]
  n0_inv(k)   = c1/n0(k)
  ! pressure [bar], pressure function takes depth in [m]
  pressz(k) = pressure(tdepth(k)*1.0E-02)
  ! vertical gradient of local potential density
  call state(salt_ref(k), temp_ref(k), pressz(k), RHO_new)
  rho_new_ref(k) = RHO_new
  call state(salt_ref(k), temp_ref(k), c0, RHO_new)
  pd_new_ref(k)  = RHO_new
enddo !k

write (*,*) '-- n0 at surface and bottom corrected by factor 2 --'


do k = 1, km
  if ( k==1 ) then
    call state(salt_ref(k  ), temp_ref(k  ), pressz(k), PD_top)
    call state(salt_ref(k+1), temp_ref(k+1), pressz(k), PD_bottom)
    n0_new(k) = ( PD_top - PD_bottom ) / dz(k) * 1.0E02
  else if ( k==km ) then
    call state(salt_ref(k-1), temp_ref(k-1), pressz(k), PD_top)
    call state(salt_ref(k  ), temp_ref(k  ), pressz(k), PD_bottom)
    n0_new(k) = ( PD_top - PD_bottom ) / dz(k) * 1.0E02
  else
    call state(salt_ref(k-1), temp_ref(k-1), pressz(k), PD_top)
    call state(salt_ref(k+1), temp_ref(k+1), pressz(k), PD_bottom)
    n0_new(k) = ( PD_top - PD_bottom ) / 2.0 / dz(k) * 1.0E02
  endif
!    n0_inv(k) = c1/n0_new(k)
enddo !k

!  to test writing out RHO, its deallocation is ommitted out here
!   deallocate( RHO, SALT, TEMP )
   deallocate( SALT )

write (*,*) ''
write (*,*) 'k      dz [m]  tdepth [m]  area [m^2]     p [bar]    salt [1]  temp [deC]'
do k = 1, km
  write (*,"( I2.1, 6ES12.3E2)") &
  k, dz(k)*1.0E-02, tdepth(k)*1.0E-02, area(k), pressz(k), salt_ref(k), temp_ref(k)
enddo !k

write (*,*) ''
write (*,*) 'k     rho_ref rho_new_ref      pd_ref  pd_new_ref          n0      n0_new'
do k = 1, km
  write (*,"( I2.1, 6ES12.3E1)") &
  k, rho_ref(k), rho_new_ref(k)-rho0, pd_ref(k), pd_new_ref(k)-rho0, &
  n0(k), n0_new(k)
enddo !k

write (*,*) ''
write (*,*) 'n0_new/n0'
do k = 1, km
  write (*,*) k, n0_new(k)/n0(k)
enddo !k

write (*,*) ''

!  test value : rho = 1.033213242 for
!  S = 35.0 PSU, theta = 20.0, pressz = 200.0
!
!  from paper:
!  [ (PSU, degC, dbar) -> kg m^-3 ]
!  (35, 25, 2000) -> 1031.654 229
!  ( 0, 20, 1000) -> 1017.726 743
!  (40, 12, 8000) -> 1062.928 258


!   call state(S, T, P, RHO_new)
!   write (*,*) '[ (PSU, degC, dbar) -> kg m^-3 ]'
!   write (*,*) '(35, 20, 2000) -> 1033.312 242'
!   call state(0.035, 20.0, 200.0, RHO_new)
!   write (*,*) RHO_new
!
!   write (*,*) '(35, 25, 2000) -> 1031.654 229'
!   call state(0.035, 25.0, 200.0, RHO_new)
!   write (*,*) RHO_new
!
!   write (*,*) '( 20, 20, 1000) -> 1017.726 743'
!   call state(0.02, 20.0, 100.0, RHO_new)
!   write (*,*) RHO_new
!
!   write (*,*) '(40, 12, 8000) -> 1062.928 258'
!   call state(0.04, 12.0, 800.0, RHO_new)
!   write (*,*) RHO_new





!===============================================================================
!  ENERGY CALCULATIONS
!  ( from now everything will be converted into SI units )
!===============================================================================

!===============================================================================
!  1. RESERVOIRS
!===============================================================================

write (*,*) ''
write (*,*) '1. RESERVOIRS'

!  all PE terms are on TTT-grid, all KE terms on TUU grid
allocate( PD2(imt,jmt,km), KE(imt,jmt,km) )
call load_3D_field(imt,jmt,km,1,nrec_PD2, PD2) ! [g^2/cm^6]
call load_3D_field(imt,jmt,km,1,nrec_KE,  KE ) ! [cm^2/s^2]

allocate(                                                                      &
rKm(imt,jmt,km), rKe(imt,jmt,km), rPm(imt,jmt,km), rPe(imt,jmt,km),            &
rKm_sint(km),    rKe_sint(km),    rPm_sint(km),    rPe_sint(km) )

do k = 1, km
  rPm(:,:,k) = -p5*g*n0_inv(k) * ( PD(:,:,k)*1.0E03  -  pd_ref(k) )**2 
  rPe(:,:,k) = -p5*g*n0_inv(k) * (   PD2(:,:,k)      -  PD(:,:,k)**2 )* 1.0E06
  rKm(:,:,k) =  p5*rho0 * ( UVEL(:,:,k)**2 + VVEL(:,:,k)**2 ) * 1.0E-04 
  rKe(:,:,k) =  rho0 * KE(:,:,k) * 1.0E-04 - rKm(:,:,k)
enddo !k

! interpolate KE terms onto TTT-grid
allocate( TTT_rKm(imt,jmt,km), TTT_rKe(imt,jmt,km) )

call uu2tt_scalar_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,rKm,TTT_rKm)
call uu2tt_scalar_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,rKe,TTT_rKe)

! integrate
call vol_int(1,1,1,imt,jmt,km,    rPm,TAREA,DZT,rPm_int)
call vol_int(1,1,1,imt,jmt,km,    rPe,TAREA,DZT,rPe_int)
call vol_int(1,1,1,imt,jmt,km,TTT_rKm,TAREA,DZT,rKm_int)
call vol_int(1,1,1,imt,jmt,km,TTT_rKe,TAREA,DZT,rKe_int)

! surface integrals
do k = 1,km
  call surf_int(1,1,imt,jmt,    rPm(:,:,k),TAREA,DZT(:,:,k),rPm_sint(k))
  call surf_int(1,1,imt,jmt,    rPe(:,:,k),TAREA,DZT(:,:,k),rPe_sint(k))
  call surf_int(1,1,imt,jmt,TTT_rKm(:,:,k),TAREA,DZT(:,:,k),rKm_sint(k))
  call surf_int(1,1,imt,jmt,TTT_rKe(:,:,k),TAREA,DZT(:,:,k),rKe_sint(k))
enddo !k

! output
300 FORMAT(A10, ES14.6E2, A2, ES10.2E2)
301 FORMAT(A10, ES14.6E2, A2, ES10.2E2, A15, ES14.6E2, A2, ES10.2E2)
write (*,300) 'rPm     =',    rPm_int, 'J',    rPm_int/1.66E23
write (*,300) 'rPe     =',    rPe_int, 'J',    rPe_int/6.38E18
write (*,300) 'rKm     =',    rKm_int, 'J',    rKm_int/1.27E18
write (*,300) 'rKe     =',    rKe_int, 'J',    rKe_int/3.55E18

if ( ntavg==1 ) then
  open(3,file='/projects/0/samoc/jan/Andree/LEC_bin_1_'//year, &
       form='unformatted',access='direct',recl=rec_length)
else if ( ntavg==5 ) then
  open(3,file='/projects/0/samoc/jan/Andree/LEC_bin_5_'//year, &
      form='unformatted',access='direct',recl=rec_length)
else if ( ntavg==11 ) then
  open(3,file='/projects/0/samoc/jan/Andree/LEC_bin_11_'//year,&
       form='unformatted',access='direct',recl=rec_length)
endif

do i = 1,4
  do k = 1,km
    if ( i==1 ) write (3,rec=(i-1)*km+k) rPm(:,:,k) 
    if ( i==2 ) write (3,rec=(i-1)*km+k) rPe(:,:,k) 
    if ( i==3 ) write (3,rec=(i-1)*km+k) rKm(:,:,k) 
    if ( i==4 ) write (3,rec=(i-1)*km+k) rKe(:,:,k) 
  enddo !k
enddo !i

! rKm not deallocated
deallocate( rKm, rKe, rPm, rPe, PD2, KE, TTT_rKm, TTT_rKe)

!===============================================================================
!  2. GENERATION
!===============================================================================

write (*,*) ''
write (*,*) '2. GENERATION'

! 2D
allocate(                                                                      &
  ALPHA01(imt,jmt), SHF(imt,jmt), BETA01(imt,jmt), SFWF(imt,jmt),              &
  STFRHO(imt,jmt), SSFRHO(imt,jmt), TAUX(imt,jmt), TAUY(imt,jmt),              &
  UTAUX(imt,jmt), VTAUY(imt,jmt) )

!read 2-D fields first
read (1,rec=nrec_ALPHA01)     ALPHA01   ! [g/cm^3/K]
read (1,rec=nrec_BETA01)      BETA01    ! [g/cm^3/msu]
read (1,rec=nrec_SHF)         SHF       ! [W/m^2]
read (1,rec=nrec_SFWF)        SFWF      ! [kg/m^2/s]
read (1,rec=nrec_STFRHO)      STFRHO    ! [C*g/cm^2/s]
read (1,rec=nrec_SSFRHO)      SSFRHO    ! [msu*g/cm^2/s]
read (1,rec=nrec_TAUX)        TAUX      ! [dyne/cm^2]
read (1,rec=nrec_TAUY)        TAUY      ! [dyne/cm^2]
read (1,rec=nrec_UTAUX)       UTAUX     ! [g/s^3]
read (1,rec=nrec_VTAUY)       VTAUY     ! [g/^3]

!
!  all PE terms quantities are on TT-grid, all KE terms on UU-grid
!

allocate(                                                                      &
  gPmh(imt,jmt), gPms(imt,jmt), gPm(imt,jmt), gPeh(imt,jmt),                   &
  gPes(imt,jmt), gPe(imt,jmt),  gKm(imt,jmt), gKe(imt,jmt),                    &
  TT_gKm(imt,jmt), TT_gKe(imt,jmt),                                            &
  gPt(imt,jmt), gKt(imt,jmt), TT_gKt(imt,jmt) )

!
!  2.1/2.2 Mean/Eddy Potential Energy Generation
!

gPmh = -g*ALPHA01*1.0E03*n0_inv(1) *                                           &
        SHF /c/rho0 * ( PD(:,:,1)*1.0E03 - pd_ref(1) )  ! []
gPeh = -g*ALPHA01*1.0E03*n0_inv(1) *                                           &
       ( STFRHO*1.0E01 - SHF /c/rho0 * pd_ref(1) ) - gPmh

gPms = -g* BETA01*1.0E03*n0_inv(1) *                                           &
        SFWF *salinity_factor*1.0E-02 * ( PD(:,:,1)*1.0E03 - pd_ref(1) )
gPes = -g* BETA01*1.0E03*n0_inv(1) *                                           &
        ( SSFRHO*1.0E01 - SFWF *salinity_factor*1.0E-02 * pd_ref(1) ) - gPms

gPm  = gPmh + gPms
gPe  = gPeh + gPes

!
!  2.3/2.4 Mean/Eddy Kinetic Energy Generation
!

gKm  = ( UVEL(:,:,1)*TAUX + VVEL(:,:,1)*TAUY )*1.0E-03
gKe  = ( UTAUX + VTAUY )*1.0E-03 - gKm

call uu2tt_scalar(imt,jmt,k,DZT(:,:,1),DZU(:,:,1),TAREA,UAREA,gKm,TT_gKm)
call uu2tt_scalar(imt,jmt,k,DZT(:,:,1),DZU(:,:,1),TAREA,UAREA,gKe,TT_gKe)

call surf_int(1,1,imt,jmt,  gPmh,TAREA,DZT(:,:,1),gPmh_int)
call surf_int(1,1,imt,jmt,  gPms,TAREA,DZT(:,:,1),gPms_int)
call surf_int(1,1,imt,jmt,   gPm,TAREA,DZT(:,:,1), gPm_int)
call surf_int(1,1,imt,jmt,  gPeh,TAREA,DZT(:,:,1),gPeh_int)
call surf_int(1,1,imt,jmt,  gPes,TAREA,DZT(:,:,1),gPes_int)
call surf_int(1,1,imt,jmt,   gPe,TAREA,DZT(:,:,1), gPe_int)
call surf_int(1,1,imt,jmt,TT_gKm,TAREA,DZT(:,:,1), gKm_int)
call surf_int(1,1,imt,jmt,TT_gKe,TAREA,DZT(:,:,1), gKe_int)

write (3,rec=169)    gPm(:,:) 
write (3,rec=170)    gPe(:,:) 
write (3,rec=171) TT_gKm(:,:) 
write (3,rec=172) TT_gKe(:,:) 
write (3,rec=173)   gPmh(:,:) 
write (3,rec=174)   gPms(:,:) 
write (3,rec=175)   gPeh(:,:) 
write (3,rec=176)   gPes(:,:) 


deallocate( gPmh, gPms, gPm, gPeh, gPes, gPe, gKm, gKe, TT_gKm, TT_gKe, TT_gKt,&
ALPHA01, SHF, BETA01, SFWF, STFRHO, SSFRHO, TAUX, TAUY, UTAUX, VTAUY, gPt, gKt )

write (*,300) '|gPmh   =', gPmh_int, 'W', gPmh_int/2.00E12
write (*,300) '|gPms   =', gPms_int, 'W', gPms_int/2.00E12
write (*,300) 'gPm     =',  gPm_int, 'W',  gPm_int/2.00E12

write (*,300) '|gPeh   =', gPeh_int, 'W', gPeh_int/ 5.8E11
write (*,300) '|gPes   =', gPes_int, 'W', gPes_int/ 5.8E11
write (*,300) 'gPe     =',  gPe_int, 'W',  gPe_int/ 5.8E11

write (*,300) 'gKm     =',  gKm_int, 'W',  gKm_int/1.85E12
write (*,300) 'gKe     =',  gKe_int, 'W',  gKe_int/2.19E12


!  =============================================================================
!  3. CONVERSION
!  =============================================================================

write (*,*) ''
write (*,*) '3. CONVERSION'

!  #############################################################################
!  ## Product Formulation ##
!  #############################################################################

allocate(                                                                      &
  UVEL2(imt,jmt,km), VVEL2(imt,jmt,km), UV(imt,jmt,km),                        &
  VW(imt,jmt,km),    UW(imt,jmt,km) ,   WVEL(imt,jmt,km),                      &
  PDU(imt,jmt,km),   PDV(imt,jmt,km) )

call load_3D_field(imt,jmt,km,1,nrec_WVEL, WVEL ) ! [cm/s]
call load_3D_field(imt,jmt,km,1,nrec_UVEL2,UVEL2) ! [cm^2/s^2]
call load_3D_field(imt,jmt,km,1,nrec_VVEL2,VVEL2) ! [cm^2/s^2]
call load_3D_field(imt,jmt,km,1,nrec_UV,   UV   ) ! [cm^2/s^2]
call load_3D_field(imt,jmt,km,1,nrec_UW,   UW   ) ! [cm^2/s^2]
call load_3D_field(imt,jmt,km,1,nrec_VW,   VW   ) ! [cm^2/s^2]
call load_3D_field(imt,jmt,km,1,nrec_PDU,  PDU  ) ! [g/cm^2/s]
call load_3D_field(imt,jmt,km,1,nrec_PDV,  PDV  ) ! [g/cm^2/s]

! calculate non-zonality parameter
! VVEL/UVEL integrated over [0E,40E]x[90S,45S]=[1100,1500]x[0,676]


call vol_int(1,1100,1,676,1500,17,VVEL2/UVEL2,UAREA,DZU,V2U2    )
call vol_int(1,1100,1,676,1500,17,DZU/DZU    ,UAREA,DZU,V2U2_vol)
V2U2 = V2U2/V2U2_vol
write (*,*) 'V2U2:', V2U2
!
!  first create new fields:
!  3.1: (tuu->ttt) uvel,vvel;  (quadratic? central difference) drho/dx, drho/dy,
!       (nabla, horizontally linear central difference) 
!          du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz
!  3.2: (tuu->ttt) uvel2,vvel2,uvelvel
!  3.3: (wtt->ttt) wvel
!  4.4: (wtt->ttt) wvelrho
!

!
! Interpolation
!

!  TEST vertical velocities
!  write (*,*) 'surface integrals of vertical velocities'
!  do k = 1,km
!   call surf_int(1,1,imt,jmt,WVEL(:,:,k),TAREA,DZT(:,:,k),WVEL_int)
!   write (*,*) k, 'WVEL_int =', WVEL_int/area(k)/1.0E02, 'm/s'
!  enddo


!
! > TUU->TTT: TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV
!

allocate(                                                                      &
  TTT_UVEL(imt,jmt,km), TTT_VVEL(imt,jmt,km), TTT_UVEL2(imt,jmt,km),           &
  TTT_VVEL2(imt,jmt,km), TTT_UV(imt,jmt,km) )
call uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,UVEL ,TTT_UVEL )
call uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,VVEL ,TTT_VVEL )
call uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,UVEL2,TTT_UVEL2)
call uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,VVEL2,TTT_VVEL2)
call uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,UV   ,TTT_UV   )
deallocate( UVEL2, VVEL2, UV )

!
! > WTT->TTT: TTT_UW, TTT_VW, TTT_WVEL, TTT_PDW
!

allocate( TTT_UW(imt,jmt,km), TTT_VW(imt,jmt,km), TTT_WVEL(imt,jmt,km) )
write (*,*) 'UW'
!call wtt2ttt(imt,jmt,km,  UW,DZT,  TTT_UW)
TTT_UW = UW
write (*,*) 'VW'
!call wtt2ttt(imt,jmt,km,  VW,DZT,  TTT_VW)
TTT_VW = VW
write (*,*) 'WVEL'
call wtt2ttt(imt,jmt,km,WVEL,DZT,TTT_WVEL)

deallocate( UW, VW )

!
! Derivatives
!

!
! > DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ
!

allocate( DUDX(imt,jmt,km), DUDY(imt,jmt,km), DUDZ(imt,jmt,km),                &
          DVDX(imt,jmt,km), DVDY(imt,jmt,km), DVDZ(imt,jmt,km) )
call nabla_hvel(imt,jmt,km,UVEL,TTT_UVEL,dz,DXT,DYT,DXU,DYU,DZT,DZU,TAREA,     &
                DUDX,DUDY,DUDZ)
call nabla_hvel(imt,jmt,km,VVEL,TTT_VVEL,dz,DXT,DYT,DXU,DYU,DZT,DZU,TAREA,     &
                DVDX,DVDY,DVDZ)
!
! > DRHODX, DRHODY
!

allocate( DPDDX(imt,jmt,km), DPDDY(imt,jmt,km) )
call grad_rho(imt,jmt,km,kmT,DZT,DZU,DXU,DYU,TAREA,UAREA,PD,DPDDX,DPDDY)

!===============================================================================

allocate(                                                                      &
  cPem(imt,jmt,km), cKem(imt,jmt,km), cPKm(imt,jmt,km), cPKe(imt,jmt,km),      &
  cPem_sint(km),    cKem_sint(km),    cPKm_sint(km),    cPKe_sint(km) )

!
!  3.1/3.2 Eddy to Mean Energy
!


do k = 1, km
  cPem(:,:,k) = - g * n0_inv(k) * 1.0E06                                       &
              * ( ( PDU(:,:,k) - PD(:,:,k) * TTT_UVEL(:,:,k)) * DPDDX(:,:,k)   &
              +   ( PDV(:,:,k) - PD(:,:,k) * TTT_VVEL(:,:,k)) * DPDDY(:,:,k) )
  cKem(:,:,k) = rho0 * 1.0E-04                                                 &
              * ( ( TTT_UVEL2(:,:,k) - TTT_UVEL(:,:,k)**2 ) * DUDX(:,:,k)      &
              +   ( TTT_VVEL2(:,:,k) - TTT_VVEL(:,:,k)**2 ) * DVDY(:,:,k)      &
              +   ( TTT_UV(:,:,k)    - TTT_UVEL(:,:,k) * TTT_VVEL(:,:,k) )     &
              *                              ( DUDY(:,:,k) + DVDX(:,:,k) )     &
              +   ( TTT_UW(:,:,k)    - TTT_UVEL(:,:,k) * TTT_WVEL(:,:,k) )     &
              *                                              DUDZ(:,:,k)       &
              +   ( TTT_VW(:,:,k)    - TTT_VVEL(:,:,k) * TTT_WVEL(:,:,k) )     &
              *                                              DVDZ(:,:,k)   )
enddo !k

deallocate( PDU, PDV, TTT_UVEL2, TTT_VVEL2, TTT_UV, TTT_UW, TTT_VW )

!
!  3.3/3.4 Potential to Kinetic Energy
!

allocate( PDW(imt,jmt,km), TTT_PDW(imt,jmt,km) )

call load_3D_field(imt,jmt,km,1,nrec_PDW,PDW)
write (*,*) 'PDW'
call wtt2ttt(imt,jmt,km, PDW,DZT, TTT_PDW)
!TTT_PDW = PDW

do k = 1, km
  cPKm(:,:,k) = -g * PD(:,:,k) * TTT_WVEL(:,:,k) * 1.0E01 
  cPKe(:,:,k) = -g * TTT_PDW(:,:,k)*1.0E01 - cPKm(:,:,k)
enddo !k

deallocate( PDW, TTT_PDW, TTT_WVEL )

!
!  volume integrals
!

call vol_int(1,1,1,imt,jmt,km,   cPem,TAREA,DZT,   cPem_int)
call vol_int(1,1,1,imt,jmt,km,   cKem,TAREA,DZT,   cKem_int)
call vol_int(1,1,1,imt,jmt,km,   cPKm,TAREA,DZT,   cPKm_int)
call vol_int(1,1,1,imt,jmt,km,   cPKe,TAREA,DZT,   cPKe_int)

!
!  surface integrals
!

do k = 1,km
  call surf_int(1,1,imt,jmt,cPem(:,:,k),TAREA,DZT(:,:,k),cPem_sint(k))
  call surf_int(1,1,imt,jmt,cKem(:,:,k),TAREA,DZT(:,:,k),cKem_sint(k))
  call surf_int(1,1,imt,jmt,cPKm(:,:,k),TAREA,DZT(:,:,k),cPKm_sint(k))
  call surf_int(1,1,imt,jmt,cPKe(:,:,k),TAREA,DZT(:,:,k),cPKe_sint(k))
enddo !k

do i = 1,5
  do k = 1,km
    if ( i==1 ) write (3,rec=176+(i-1)*km+k) cPem(:,:,k) 
    if ( i==2 ) write (3,rec=176+(i-1)*km+k) cKem(:,:,k) 
    if ( i==3 ) write (3,rec=176+(i-1)*km+k) cPKm(:,:,k) 
    if ( i==4 ) write (3,rec=176+(i-1)*km+k) cPKe(:,:,k) 
    if ( i==5 ) write (3,rec=176+(i-1)*km+k) TEMP(:,:,k) 
  enddo !k
enddo !i

deallocate( cPKm, cPKe )

write (*,300) 'cPem    =',    cPem_int, 'W',    cPem_int/-8.3E11
write (*,300) 'cKem    =',    cKem_int, 'W',    cKem_int/-1.1E11
write (*,300) 'cPKm    =',    cPKm_int, 'W',    cPKm_int/-4.9E11
write (*,300) 'cPKe    =',    cPKe_int, 'W',    cPKe_int/ 7.3E11

write (*,*) ''
write (*,*) '4. DISSIPATION'

dPm = gPm_int + cPem_int - cPKm_int
dPe = gPe_int - cPem_int - cPKe_int
dKm = gKm_int + cKem_int + cPKm_int
dKe = gKe_int - cKem_int + cPKe_int

write (*,300) 'dPm     =', dPm   , 'W',    dPm/1.66E12
write (*,300) 'dPe     =', dPe   , 'W',    dPe/6.80E11
write (*,300) 'dKm     =', dKm   , 'W',    dKm/1.36E12
write (*,300) 'dKe     =', dKe   , 'W',    dKe/3.03E12

100 FORMAT (21(A,","),A)
101 FORMAT (21(E15.7,","),E15.7)
if ( ntavg==1 ) then
  open  (2,file=('output/total_LEC_1_'//year)//'.out',form='formatted')
else if ( ntavg==5 ) then
  open  (2,file=('output/total_LEC_5_'//year)//'.out',form='formatted')
else if ( ntavg==11 ) then
  open  (2,file=('output/total_LEC_11_'//year)//'.out',form='formatted')
endif
write (2,100) 'pd_avg',                                                        &
              'rPm','rPe','rKm','rKe',                                         &
              'gPmh','gPms','gPm',                                             &
              'gPeh','gPes','gPe',                                             &
              'gKm','gKe',                                                     &
              'cPem','cKem','cPKm','cPKe',                                     &
              'dPm','dPe','dKm','dKe',                                         &
              'V2U2'
write (2,101) pd_avg,                                                          &
              rPm_int,     rPe_int,     rKm_int,     rKe_int,                  &
              gPmh_int,    gPms_int,    gPm_int,                               &
              gPeh_int,    gPes_int,    gPe_int,                               &
              gKm_int,     gKe_int,                                            &
              cPem_int,    cKem_int,    cPKm_int,    cPKe_int,                 &
              dPm,         dPe,         dKm,         dKe,                      &
              V2U2
close (2)

! ##############################################################################
! ## Flux Formulation ##
! ##############################################################################

! allocate(                                                                   &
! WUU_WVEL(imt,jmt,km),                                                       &
! WTU_k(imt,jmt),           WTV_k(imt,jmt),                                   &
! WUU_WTU_eddy_k(imt,jmt),  WUU_WTV_eddy_k(imt,jmt),                          &
! WTT_WTU_eddy_k(imt,jmt),  WTT_WTV_eddy_k(imt,jmt),                          &
! WTT_WTU_eddy(imt,jmt,km), WTT_WTV_eddy(imt,jmt,km),                         &
! TTT_WTU_eddy(imt,jmt,km), TTT_WTV_eddy(imt,jmt,km) )
! 
! do k = 1, km
!   read (1,rec=nrec_WTU+k-1) WTU_k    ! [cm/s^2]
!   read (1,rec=nrec_WTV+k-1) WTV_k    ! [cm/s^2]
! 
!   do j = 2,jmt
!     WUU_WVEL(1,j,k) = p25 * ( WVEL(1  ,j  ,k) * TAREA(1  ,j  )                &
!                             + WVEL(1  ,j-1,k) * TAREA(1  ,j-1)                &
!                             + WVEL(imt,j  ,k) * TAREA(imt,j  )                &
!                             + WVEL(imt,j-1,k) * TAREA(imt,j-1))               &
!                           / UAREA(1,j)
!     do i = 2,imt
!       if ( j.eq.2 ) then
!         WUU_WVEL(i,1,k) = p25 * ( WVEL(i  ,1  ,k) * TAREA(i  ,1  )            &
!                                 + WVEL(i-1,1  ,k) * TAREA(i-1,1  ))           &
!                               / UAREA(i,1)
!       endif
!     WUU_WVEL(i,j,k) = p25 * ( WVEL(i  ,j  ,k) * TAREA(i  ,j  )               &
!                             + WVEL(i  ,j-1,k) * TAREA(i  ,j-1)               &
!                             + WVEL(i-1,j  ,k) * TAREA(i-1,j  )               &
!                             + WVEL(i-1,j-1,k) * TAREA(i-1,j-1))              &
!                           / UAREA(i,j)
!     enddo !i
!   enddo !j
! 
!   do j = 1,jmt
!     do i = 1,imt
!       if ( DZU(i,j,k).ne.c0 ) then
!         if ( k.eq.1 ) then
!           WUU_WTU_eddy_k(i,j) = WTU_k(i,j) - WUU_WVEL(i,j,k)*UVEL(i,j,k)/DZU(i,j,k)
!           WUU_WTV_eddy_k(i,j) = WTV_k(i,j) - WUU_WVEL(i,j,k)*VVEL(i,j,k)/DZU(i,j,k)
!         else
!           WUU_WTU_eddy_k(i,j) = WTU_k(i,j) - p5 * WUU_WVEL(i,j,k)                &
!                                     * (UVEL(i,j,k) + UVEL(i,j,k-1))/DZU(i,j,k)
!           WUU_WTV_eddy_k(i,j) = WTV_k(i,j) - p5 * WUU_WVEL(i,j,k)                &
!                                       * (VVEL(i,j,k) + VVEL(i,j,k-1))/DZU(i,j,k)
!         endif
!       endif
!     enddo !j
!   enddo !i
! 
!   call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,             &
!                       WUU_WTU_eddy_k(:,:),WTT_WTU_eddy_k)
!   call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,             &
!                       WUU_WTV_eddy_k(:,:),WTT_WTV_eddy_k)
!     
!   WTT_WTU_eddy(:,:,k) =  WTT_WTU_eddy_k(:,:)
!   WTT_WTV_eddy(:,:,k) =  WTT_WTV_eddy_k(:,:)
! 
! enddo !k
! 
!    
! write (*,*) 'after first k loop'
! 
! call wtt2ttt(imt,jmt,km,WTT_WTU_eddy,DZT,TTT_WTU_eddy)
! call wtt2ttt(imt,jmt,km,WTT_WTV_eddy,DZT,TTT_WTV_eddy)
! 
! allocate(                                                                   &
! UUE_k(imt,jmt), VUN_k(imt,jmt),                                             &
! UEU_k(imt,jmt), UEV_k(imt,jmt),                                             &
! VNU_k(imt,jmt), VNV_k(imt,jmt),                                             &
! UTE_k(imt,jmt), VTN_k(imt,jmt),                                             &
! UPD_k(imt,jmt), VPD_k(imt,jmt),                                             &
! TTT_UEU_k(imt,jmt), TTT_UEV_k(imt,jmt),                                     &
! TTT_VNU_k(imt,jmt), TTT_VNV_k(imt,jmt),                                     &
! TTT_UPD_k(imt,jmt), TTT_VPD_k(imt,jmt),                                     &
! UEU_eddy_k(imt,jmt), UEV_eddy_k(imt,jmt),                                   &
! VNU_eddy_k(imt,jmt), VNV_eddy_k(imt,jmt),                                   & 
! UPD_eddy_k(imt,jmt), VPD_eddy_k(imt,jmt),                                   & 
! TTT_UEU_eddy_k(imt,jmt), TTT_UEV_eddy_k(imt,jmt),                           &
! TTT_VNU_eddy_k(imt,jmt), TTT_VNV_eddy_k(imt,jmt),                           & 
! TTT_UPD_eddy_k(imt,jmt), TTT_VPD_eddy_k(imt,jmt) )
! 
! write (*,*) 'just before k loop'
! do k = 1, km
!   read (1,rec=nrec_UEU+k-1) UEU_k    ! [cm/s^2]
!   read (1,rec=nrec_VNU+k-1) VNU_k    ! [cm/s^2]
!   read (1,rec=nrec_UEV+k-1) UEV_k    ! [cm/s^2]
!   read (1,rec=nrec_VNV+k-1) VNV_k    ! [cm/s^2]
! 
! 
!   ! interpolate to TT-grid
!   !call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,WTU_k(:,:),TTT_WTU_k)
!   !call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,WTV_k(:,:),TTT_WTV_k)
! 
!   !call interp_mom_fluxes(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),DXT,DXU,DYT,DYU,UEU_k,VNU_k,TTT_UEU_k,TTT_VNU_k)
!   !call interp_mom_fluxes(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),DXT,DXU,DYT,DYU,UEV_k,VNV_k,TTT_UEV_k,TTT_VNV_k)
! 
!  ! edge cases!
!   do j = 2,jmt-1
!     UUE_k(1  ,j) = p25 *(UVEL(2  ,j  ,k)*DYU(2  ,j  )*DZU(2  ,j  ,k) + &
!                          UVEL(1  ,j  ,k)*DYU(1  ,j  )*DZU(1  ,j  ,k))+ &
!                    p125*(UVEL(2  ,j-1,k)*DYU(2  ,j-1)*DZU(2  ,j-1,k) + &
!                          UVEL(1  ,j-1,k)*DYU(1  ,j-1)*DZU(1  ,j-1,k) + &
!                          UVEL(2  ,j+1,k)*DYU(2  ,j+1)*DZU(2  ,j+1,k) + &
!                          UVEL(1  ,j+1,k)*DYU(1  ,j+1)*DZU(1  ,j+1,k))
!     VUN_k(1  ,j) = p25 *(VVEL(1  ,j+1,k)*DXU(1  ,j+1)*DZU(1  ,j+1,k) + &
!                          VVEL(1  ,j  ,k)*DXU(1  ,j  )*DZU(1  ,j  ,k))+ &
!                    p125*(VVEL(imt,j+1,k)*DXU(imt,j+1)*DZU(imt,j+1,k) + &
!                          VVEL(imt,j  ,k)*DXU(imt,j  )*DZU(imt,j  ,k) + &
!                          VVEL(2  ,j+1,k)*DXU(2  ,j+1)*DZU(2  ,j+1,k) + &
!                          VVEL(2  ,j  ,k)*DXU(2  ,j  )*DZU(2  ,j  ,k))
! 
!     UUE_k(imt,j) = p25 *(UVEL(1    ,j  ,k)*DYU(1    ,j  )*DZU(1    ,j  ,k) + &
!                          UVEL(imt  ,j  ,k)*DYU(imt  ,j  )*DZU(imt  ,j  ,k))+ &
!                    p125*(UVEL(1    ,j-1,k)*DYU(1    ,j-1)*DZU(1    ,j-1,k) + &
!                          UVEL(imt  ,j-1,k)*DYU(imt  ,j-1)*DZU(imt  ,j-1,k) + &
!                          UVEL(1    ,j+1,k)*DYU(1    ,j+1)*DZU(1    ,j+1,k) + &
!                          UVEL(imt  ,j+1,k)*DYU(imt  ,j+1)*DZU(imt  ,j+1,k))
!     VUN_k(imt,j) = p25 *(VVEL(imt  ,j+1,k)*DXU(imt  ,j+1)*DZU(imt  ,j+1,k) + &
!                          VVEL(imt  ,j  ,k)*DXU(imt  ,j  )*DZU(imt  ,j  ,k))+ &
!                    p125*(VVEL(imt-1,j+1,k)*DXU(imt-1,j+1)*DZU(imt-1,j+1,k) + &
!                          VVEL(imt-1,j  ,k)*DXU(imt-1,j  )*DZU(imt-1,j  ,k) + &
!                          VVEL(1    ,j+1,k)*DXU(1    ,j+1)*DZU(1    ,j+1,k) + &
!                          VVEL(1    ,j  ,k)*DXU(1    ,j  )*DZU(1    ,j  ,k))
!     do i = 2,imt-1
!       UUE_k(i,j) = p25 *(UVEL(i+1,j  ,k)*DYU(i+1,j  )*DZU(i+1,j  ,k) + &
!                          UVEL(i  ,j  ,k)*DYU(i  ,j  )*DZU(i  ,j  ,k))+ &
!                    p125*(UVEL(i+1,j-1,k)*DYU(i+1,j-1)*DZU(i+1,j-1,k) + &
!                          UVEL(i  ,j-1,k)*DYU(i  ,j-1)*DZU(i  ,j-1,k) + &
!                          UVEL(i+1,j+1,k)*DYU(i+1,j+1)*DZU(i+1,j+1,k) + &
!                          UVEL(i  ,j+1,k)*DYU(i  ,j+1)*DZU(i  ,j+1,k))
!       VUN_k(i,j) = p25 *(VVEL(i  ,j+1,k)*DXU(i  ,j+1)*DZU(i  ,j+1,k) + &
!                          VVEL(i  ,j  ,k)*DXU(i  ,j  )*DZU(i  ,j  ,k))+ &
!                    p125*(VVEL(i-1,j+1,k)*DXU(i-1,j+1)*DZU(i-1,j+1,k) + &
!                          VVEL(i-1,j  ,k)*DXU(i-1,j  )*DZU(i-1,j  ,k) + &
!                          VVEL(i+1,j+1,k)*DXU(i+1,j+1)*DZU(i+1,j+1,k) + &
!                          VVEL(i+1,j  ,k)*DXU(i+1,j  )*DZU(i+1,j  ,k))
!      enddo !i
!   enddo !j
! 
!   do j = 2,jmt-1
!     do i = 1,imt-1
!       if ( DZU(i,j,k).ne.c0 ) then
!         UEU_eddy_k(i,j) = UEU_k(i,j) - p5*UUE_k(i,j) / UAREA(i,j) / DZU(i,j,k)  &
!                           * ( UVEL(i,j,k) + UVEL(i+1,j  ,k) )
!         UEV_eddy_k(i,j) = UEV_k(i,j) - p5*UUE_k(i,j) / UAREA(i,j) / DZU(i,j,k)  &
!                           * ( VVEL(i,j,k) + VVEL(i+1,j  ,k) )
!         VNU_eddy_k(i,j) = VNU_k(i,j) - p5*VUN_k(i,j) / UAREA(i,j) / DZU(i,j,k)  &
!                           * ( UVEL(i,j,k) + UVEL(i  ,j+1,k) )
!         VNV_eddy_k(i,j) = VNV_k(i,j) - p5*VUN_k(i,j) / UAREA(i,j) / DZU(i,j,k)  &
!                           * ( VVEL(i,j,k) + VVEL(i  ,j+1,k) )
!       endif
!     enddo !i
!   enddo !j
! 
!   ! interpolate to TT-grid
!   call interp_mom_fluxes(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),DXT,DXU,DYT,DYU,    &
!                          UEU_eddy_k,VNU_eddy_k,TTT_UEU_eddy_k,TTT_VNU_eddy_k)
!   call interp_mom_fluxes(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),DXT,DXU,DYT,DYU,    &
!                          UEV_eddy_k,VNV_eddy_k,TTT_UEV_eddy_k,TTT_VNV_eddy_k)
! 
!   cKem(:,:,k) = rho0 * 1.0E-04 *                                             &
!                           ( TTT_UEU_eddy_k(:,:) * DXT(:,:)   * DUDX(:,:,k) + &
!                             TTT_UEV_eddy_k(:,:) * DXT(:,:)   * DVDX(:,:,k) + &
!                             TTT_VNU_eddy_k(:,:) * DYT(:,:)   * DUDY(:,:,k) + &
!                             TTT_VNV_eddy_k(:,:) * DYT(:,:)   * DVDY(:,:,k) + &
!                             TTT_WTU_eddy(:,:,k) * DZT(:,:,k) * DUDZ(:,:,k) + &
!                             TTT_WTV_eddy(:,:,k) * DZT(:,:,k) * DVDZ(:,:,k) )
!    
! !===============================================================================
!    ! cPEm
! !===============================================================================
! 
!   read (1,rec=nrec_UPD+k-1) UPD_k    ! [g/cm^3/s]
!   read (1,rec=nrec_VPD+k-1) VPD_k    ! [g/cm^3/s]
! 
!   do j = 2,jmt-1
!     UTE_k(1  ,j) = p5 * ( UVEL(1  ,j  ,k)*DYU(1  ,j  )*DZU(1  ,j  ,k)         &
!                         + UVEL(1  ,j-1,k)*DYU(1  ,j-1)*DZU(1  ,j-1,k) )
!     VTN_k(1  ,j) = p5 * ( VVEL(1  ,j  ,k)*DXU(1  ,j  )*DZU(1  ,j  ,k)         &
!                         + VVEL(imt,j  ,k)*DXU(imt,j  )*DZU(imt,j  ,k) )
! 
!     UTE_k(imt,j) = p5 * ( UVEL(imt  ,j  ,k)*DYU(imt  ,j  )*DZU(imt  ,j  ,k)   &
!                         + UVEL(imt  ,j-1,k)*DYU(imt  ,j-1)*DZU(imt  ,j-1,k) )
!     VTN_k(imt,j) = p5 * ( VVEL(imt  ,j  ,k)*DXU(imt  ,j  )*DZU(imt  ,j  ,k)   &
!                         + VVEL(imt-1,j  ,k)*DXU(imt-1,j  )*DZU(imt-1,j  ,k) )
! 
!     do i = 2,imt-1
!       UTE_k(i,j)  = p5 * ( UVEL(i  ,j  ,k)*DYU(i  ,j  )*DZU(i  ,j  ,k)         &
!                          + UVEL(i  ,j-1,k)*DYU(i  ,j-1)*DZU(i  ,j-1,k) )
!       VTN_k(i,j)  = p5 * ( VVEL(i  ,j  ,k)*DXU(i  ,j  )*DZU(i  ,j  ,k)         &
!                          + VVEL(i-1,j  ,k)*DXU(i-1,j  )*DZU(i-1,j  ,k) )
!     enddo !i
!   enddo !j
! 
!   do j = 2,jmt-1
!     do i = 1,imt-1
!       if ( DZT(i,j,k).ne.c0 ) then
!         UPD_eddy_k(i,j) = UPD_k(i,j) - p5 * UTE_k(i,j) / TAREA(i,j)           &
!                                    / DZT(i,j,k) * ( PD(i,j,k) + PD(i+1,j  ,k) )
!         VPD_eddy_k(i,j) = VPD_k(i,j) - p5 * VTN_k(i,j) / TAREA(i,j)           &
!                                    / DZT(i,j,k) * ( PD(i,j,k) + PD(i  ,j-1,k) )
!       endif
!       if ( ISNAN(UPD_eddy_k(i,j)) .and. k.lt.3 ) then
!         write (*,*) i,j,k, 'UPD_eddy_k is NaN'
!       endif 
!       if ( ISNAN(VPD_eddy_k(i,j)) .and. k.lt.3 ) then
!         write (*,*) i,j,k, 'VPD_eddy_k is NaN'
!       endif 
!     enddo !i
!   enddo !j
! 
!     ! interpolate UPD/VPD to TTT grid
!   do j = 2, jmt
!     do i = 1, imt
!       if ( DZT(i,j,k).ne.c0 ) then
!         if ( i.eq.1 ) then
!           TTT_UPD_eddy_k(1,j) = p5 * ( UPD_eddy_k(1,j) + UPD_eddy_k(imt,j) )
!         else
!           TTT_UPD_eddy_k(i,j) = p5 * ( UPD_eddy_k(i,j) + UPD_eddy_k(i-1,j) )
!         endif
!           TTT_VPD_eddy_k(i,j) = p5 * ( VPD_eddy_k(i,j) + VPD_eddy_k(i,j-1) )
!       endif
!       if ( ISNAN(TTT_UPD_eddy_k(i,j)) .and. k.lt.3 ) then
!         write (*,*) i,j,k, 'TTT_UPD_eddy_k is NaN'
!       endif 
!       if ( ISNAN(TTT_VPD_eddy_k(i,j)) .and. k.lt.3 ) then
!         write (*,*) i,j,k, 'TTT_VPD_eddy_k is NaN'
!       endif 
!     enddo !i
!   enddo !j
! 
!   cPem(:,:,k) = -g * n0_inv(k) * 1.0E06 *                                    &
!                             ( TTT_UPD_eddy_k(:,:) * DXT(:,:) * DPDDX(:,:,k)  &
!                             + TTT_VPD_eddy_k(:,:) * DYT(:,:) * DPDDY(:,:,k) )
! 
! enddo !k
! 
! call vol_int(1,1,1,imt,jmt,km,cKem,TAREA,DZT,cKem_int)
! call vol_int(1,1,1,imt,jmt,km,cPem,TAREA,DZT,cPem_int)
! 
! write (*,*) '------------------------'
! write (*,*) 'FLUX FORMULATION RESULTS'
! write (*,*) ''
! 
! write (*,300) 'cPem    =',    cPem_int, 'W',    cPem_int/-8.3E11
! write (*,300) 'cKem    =',    cKem_int, 'W',    cKem_int/-1.1E11
! write (*,300) 'cPKm    =',    cPKm_int, 'W',    cPKm_int/-4.9E11
! write (*,300) 'cPKe    =',    cPKe_int, 'W',    cPKe_int/ 7.3E11
! 
! !  closing the output file
! close(1)
! 
! !===============================================================================
! !  4. DISSIPATION
! !===============================================================================
! 
! write (*,*) ''
! write (*,*) '4. DISSIPATION'
! 
! dPm    = gPm_int + cPem_int - cPKm_int
! dPe    = gPe_int - cPem_int - cPKe_int
! dKm    = gKm_int + cKem_int + cPKm_int
! dKe    = gKe_int - cKem_int + cPKe_int
! 
! write (*,300) 'dPm     =', dPm   , 'W',    dPm/1.66E12
! write (*,300) 'dPe     =', dPe   , 'W',    dPe/6.80E11
! write (*,300) 'dKm     =', dKm   , 'W',    dKm/1.36E12
! write (*,300) 'dKe     =', dKe   , 'W',    dKe/3.03E12
! 
! !===============================================================================
! !  TESTING
! !===============================================================================

!   write (*,*) ''
!   write (*,*) 'TESTING'

!   write (*,*) 'surface integrals of vertical velocities'
!   do k = 1,km
!   call surf_int(1,1,imt,jmt,WVEL(:,:,k),TAREA,DZT(:,:,k),WVEL_int)
!   write (*,*) k, 'WVEL_int =', WVEL_int/area(k)/1.0E02, 'm/s'
!   enddo
!   >>> ca. 10^-13 m/s, seems correct

!   write (*,*) 'surface integrals of TTT-vertical velocities'
!   do k = 1,km
!   call surf_int(1,1,imt,jmt,TTT_WVEL(:,:,k),TAREA,DZT(:,:,k),WVEL_int)
!   write (*,*) k, 'WVEL_int =', WVEL_int/area(k)/1.0E02, 'm/s'
!   enddo
!   >>> same as above

!   write (*,*) 'volume integral of vertical velocity'
!   call vol_int(1,1,1,imt,jmt,km,TTT_WVEL,TAREA,DZT,WVEL_int)
!   write (*,*) 'WVEL_int =', WVEL_int/volume/1.0E02, 'm/s'
!   >>> 8.6*10^-8 m/s, seems correct
!   >>> about the same for volume integral of WVEL, instead of TTT_WVEL


!===============================================================================
!   FILE OUTPUT
!===============================================================================

! closing the output binary file that includes all fields
close(3)

write (*,*) ''
write (*,*) 'OUTPUT'

! open  (1,file=('output/total_LEC_ff_'//year)//'.out',form='formatted')
! write (1,100) 'pd_avg',                                                        &
!               'rPm','rPe','rKm','rKe',                                         &
!               'gPmh','gPms','gPm',                                             &
!               'gPeh','gPes','gPe',                                             &
!               'gKm','gKe',                                                     &
!               'cPem','cKem','cPKm','cPKe',                                     &
!               'dPm','dPe','dKm','dKe',                                         &
!               'V2U2'
! write (1,101) pd_avg,                                                          &
!               rPm_int,     rPe_int,     rKm_int,     rKe_int,                  &
!               gPmh_int,    gPms_int,    gPm_int,                               &
!               gPeh_int,    gPes_int,    gPe_int,                               &
!               gKm_int,     gKe_int,                                            &
!               cPem_int,    cKem_int,    cPKm_int,    cPKe_int,                 &
!               dPm,          dPe,        dKm,         dKe,                      &
!               V2U2
! close (1)


200 FORMAT (13(A,","),A)
201 FORMAT (13(E15.7,","),E15.7)
if ( ntavg==1 ) then 
  open  (2,file='output/level_LEC_1_'//(year//'.out'),form='formatted')
else if ( ntavg==5 ) then 
  open  (2,file='output/level_LEC_5_'//(year//'.out'),form='formatted')
else if ( ntavg==11 ) then 
  open  (2,file='output/level_LEC_11_'//(year//'.out'),form='formatted')
endif
write (2,200) 'k','depth_k','area_k','rho_ref_k','pd_ref_k','n0_k',            &
              'rPm_k','rPe_k','rKm_k','rKe_k',                                 &
              'cPem_k','cKem_k','cPKm_k','cPKe_k'
do k = 1,km
  !surface integrals of energy budget terms
  write (2,201) real(k), tdepth(k), area(k), rho_ref(k), pd_ref(k), n0(k),     &
                rPm_sint(k),  rPe_sint(k),  rKm_sint(k),  rKe_sint(k),         &
                cPem_sint(k), cKem_sint(k), cPKm_sint(k), cPKe_sint(k)
enddo
close (2)

write (*,*) '>>> done writing output'

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

double precision, intent(in) :: depth    ! depth in meters

! !OUTPUT PARAMETERS:

double precision :: pressure   ! pressure [bars]

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
!
!  convert depth in meters to pressure in bars
!
!-----------------------------------------------------------------------

pressure = pc1*(exp(-pc2*depth) - c1) + pc3*depth + pc4*depth**2

!-----------------------------------------------------------------------
!EOC

end function pressure


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

end subroutine load_3d_field

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine load_2D_field(imt,jmt,nfile,nrec,FIELD)
!
!  loads 2D field (starting with record number nrec) from file with number nfile
!
implicit none

integer, intent(in)                      :: imt, jmt, nfile,  nrec
real,    dimension(imt,jmt), intent(out) :: FIELD

read (nfile,rec=nrec) FIELD

end subroutine load_2D_field

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
subroutine grad_rho(imt,jmt,km,kmT,DZT,DZU,DXU,DYU,TAREA,UAREA,RHO,      &
                          TT_DRHODX,TT_DRHODY)
!   call grad_rho(imt,jmt,km,DZT,DZU,DXU,DYU,TAREA,UAREA,PD,DPDDX,DPDDY)
!
!     calculates gradient components of RHO/PD
!     as described in POP reference manual
!     RHO/PD on 3111 grid, output on 3221 grid
!     units [g/cm^4]
!
implicit none

! input/output variables
integer,                                 intent(in)  :: imt,jmt,km
integer,          dimension(imt,jmt   ), intent(in)  :: kmT
double precision, dimension(imt,jmt   ), intent(in)  :: DXU, DYU ! [cm]
double precision, dimension(imt,jmt   ), intent(in)  :: TAREA, UAREA ! [cm^2]
real,             dimension(imt,jmt,km), intent(in)  :: DZT, DZU, RHO ! [cm], [g/cm^3]
real,             dimension(imt,jmt,km)              :: DRHODX, DRHODY
real,             dimension(imt,jmt,km), intent(out) :: TT_DRHODX, TT_DRHODY
! internal variables
real,             parameter                         :: p5=0.5, c0=0.0


write (*,*) 'grad rho started'

DRHODX = c0
DRHODY = c0

do k = 1,km

  ! northermost line between poles of tripolar grid
  do i = 1,imt-1
    if ( DZU(i,jmt,k).ne.c0 ) then
      if ( i.lt.imt/2 ) then      ! first half of array, before wrapping around
        DRHODX(i,jmt,k) = ( ( RHO(i+1    ,jmt,k) + RHO(imt-i  ,jmt,k) )       &
                          - ( RHO(i      ,jmt,k) + RHO(imt-i+1,jmt,k) ) )     &
                          / 2 / DXU(i,jmt)
        DRHODY(i,jmt,k) = ( ( RHO(imt-i+1,jmt,k) + RHO(imt-i  ,jmt,k) )       &
                          - ( RHO(i      ,jmt,k) + RHO(i+1    ,jmt,k) ) )     &
                          / 2 / DYU(i,jmt)
      elseif ( i.gt.imt/2 ) then  ! second half, after wrap
        DRHODX(i,jmt,k) = DRHODX(imt-i,jmt,k)
        DRHODY(i,jmt,k) = DRHODY(imt-i,jmt,k)
      endif
    endif
  enddo !i
  
  do j = 1,jmt-1
  
  ! eastern boundary
    if ( DZU(imt,j,k).ne.c0 ) then
      DRHODX(imt,j,k) = ( ( RHO(1  ,j  ,k) + RHO(1  ,j+1,k) )               &
                        - ( RHO(imt,j  ,k) + RHO(imt,j+1,k) ) )             &
                        / 2 / DXU(imt,j)
      DRHODY(imt,j,k) = ( ( RHO(imt,j+1,k) + RHO(1  ,j+1,k) )               &
                        - ( RHO(imt,j  ,k) + RHO(1  ,j  ,k) ) )             &
                        / 2 / DYU(imt,j)
    endif
  ! rest of arrays
    do i = 1,imt-1
      if ( DZU(i,j,k).ne.c0 ) then
        DRHODX(i,j,k) = ( ( RHO(i+1,j  ,k) + RHO(i+1,j+1,k) )                &
                        - ( RHO(i  ,j  ,k) + RHO(i  ,j+1,k) ) )              &
                        / 2 / DXU(i,j)
        DRHODY(i,j,k) = ( ( RHO(i  ,j+1,k) + RHO(i+1,j+1,k) )                &
                        - ( RHO(i  ,j  ,k) + RHO(i+1,j  ,k) ) )              &
                        / 2 / DYU(i,j)
      endif
    enddo !i
  enddo !j
        
  ! Testing
  do i = 1,imt
    do j = 1,jmt 
      if ( ISNAN(DRHODX(i,j,k)) ) then
        write (*,*) i,j,k, 'DRHODX is NaN'
      endif
      if ( ISNAN(DRHODY(i,j,k)) ) then
        write (*,*) i,j,k, 'DRHODY is NaN'
        if ( i.gt.3595 ) then
          write (*,*) DRHODX(i,j,k), DRHODY(i,j,k), kmT(i,j)
          write (*,*) DZT(i,j,k), DZU(i,j,k), DYU(i,j)
          write (*,*) RHO(i,j+1,k), RHO(i+1,j+1,k)
          write (*,*) RHO(i,j,k), RHO(i+1,j,k)
          write (*,*) ' '
          !DRHODY(i,j,k) = c0
        endif
      endif
      if ( DRHODY(i,j,k) == 1.0/0 .or. DRHODY(i,j,k) == -1/0.0 ) then
        write (*,*) i,j,k, 'DRHODY is +/- Infinity'
        DRHODY(i,j,k) = c0
      endif
    enddo !j
  enddo !i
  
  call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,&
                     DRHODX(:,:,k),TT_DRHODX(:,:,k))
  call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,&
                     DRHODY(:,:,k),TT_DRHODY(:,:,k))
  
  ! Testing
  do j = 1,jmt
    do i = 1,imt 
      if ( ISNAN(TT_DRHODX(i,j,k)) ) then
        write (*,*) i,j,k, 'TTT_DRHODX is NaN'
        write (*,*)  DZT(i,j,k), DZU(i,j,k), kmT(i,j), DXU(i,j)
        write (*,*)  DRHODX(i-1,j,k), DRHODX(i,j,k)
        write (*,*)  DRHODX(i-1,j-1,k), DRHODX(i,j-1,k)
        TT_DRHODX(i,j,k) = c0
      endif
      if ( ISNAN(TT_DRHODY(i,j,k)) ) then
        write (*,*) i,j,k, 'TTT_DRHODY is NaN'
        write (*,*)  DZT(i,j,k), DZU(i,j,k), kmT(i,j), DYU(i,j)
        write (*,*)  DRHODY(i-1,j,k), DRHODY(i,j,k)
        write (*,*)  DRHODY(i-1,j-1,k), DRHODY(i,j-1,k)
        TT_DRHODY(i,j,k) = c0
      endif
    enddo !i
  enddo !j

enddo !k

end subroutine grad_rho

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
subroutine nabla_hvel(imt,jmt,km,VEL,TTT_VEL,dz,DXT,DYT,DXU,DYU,DZT,DZU,&
                            TAREA,DVELDX,DVELDY,DVELDZ)
!
!     calculates 3D gradients of a horizontal velocity field
!
implicit none

! input/output variables
integer,                                 intent(in)  :: imt,jmt,km
double precision, dimension(        km), intent(in)  :: dz
double precision, dimension(imt,jmt   ), intent(in)  :: DXT, DYT, DXU, DYU, TAREA
real            , dimension(imt,jmt,km), intent(in)  :: DZT, DZU
real,             dimension(imt,jmt,km), intent(in)  :: VEL, TTT_VEL
real,             dimension(imt,jmt,km), intent(out) :: DVELDX, DVELDY, DVELDZ ! TTT gradients
! internal variables
real,             parameter                         :: p5=0.5, c0=0.0

write (*,*) 'hdiv vel started'

! horizontal gradients, central difference,
! eq. 3.6 page 16 in 2010 POP Reference Manual
do k = 1,km
! southernmost gridpoints (j=1) on TTT-grid must be 0,
! as there would not be UU gridpoints south of it otherwise
  do j = 2,jmt
    ! westernmost gridpoints (i=1)
    if ( DZT(1,j,k).ne.c0 ) then
      DVELDX(1,j,k) = p5*( ( DYU(1  ,j  )*DZU(1  ,j  ,k)*VEL(1  ,j  ,k)     &
                           + DYU(1  ,j-1)*DZU(1  ,j-1,k)*VEL(1  ,j-1,k) )   &
                         - ( DYU(imt,j  )*DZU(imt,j  ,k)*VEL(imt,j  ,k)     &
                           + DYU(imt,j-1)*DZU(imt,j-1,k)*VEL(imt,j-1,k) ) ) &
                         / TAREA(1,j) / DZT(1,j,k)
      DVELDY(1,j,k) = p5*( ( DXU(imt,j  )*DZU(imt,j  ,k)*VEL(imt,j  ,k)     &
                           + DXU(1  ,j  )*DZU(1  ,j  ,k)*VEL(1  ,j  ,k) )   &
                         - ( DXU(imt,j-1)*DZU(imt,j-1,k)*VEL(imt,j-1,k)     &
                           + DXU(1  ,j-1)*DZU(1  ,j-1,k)*VEL(1  ,j-1,k) ) ) &
                         / TAREA(1,j) / DZT(1,j,k)
    endif
    ! all gridpoints except western- and southernmost (i=1, j=1)
    do i = 2,imt
      if ( DZT(i,j,k).ne.c0 ) then
        DVELDX(i,j,k) = p5*( ( DYU(i  ,j  )*DZU(i  ,j  ,k)*VEL(i  ,j  ,k)    &
                             + DYU(i  ,j-1)*DZU(i  ,j-1,k)*VEL(i  ,j-1,k) )  &
                           - ( DYU(i-1,j  )*DZU(i-1,j  ,k)*VEL(i-1,j  ,k)    &
                             + DYU(i-1,j-1)*DZU(i-1,j-1,k)*VEL(i-1,j-1,k) )) &
                           / TAREA(i,j) / DZT(i,j,k)
        DVELDY(i,j,k) = p5*( ( DXU(i-1,j  )*DZU(i-1,j  ,k)*VEL(i-1,j  ,k)    &
                             + DXU(i  ,j  )*DZU(i  ,j  ,k)*VEL(i  ,j  ,k) )  &
                           - ( DXU(i-1,j-1)*DZU(i-1,j-1,k)*VEL(i-1,j-1,k)    &
                             + DXU(i  ,j-1)*DZU(i  ,j-1,k)*VEL(i  ,j-1,k) )) &
                           / TAREA(i,j) / DZT(i,j,k)
      endif
    enddo !i
  enddo !j
enddo !k

! vertical gradient
! top/bottom layer
do j = 1,jmt
  do i = 1,imt
    if ( DZT(i,j, 1).ne.c0 .and. DZT(i,j,2).ne.c0 ) then
      DVELDZ(i,j, 1) = 2*( TTT_VEL(i,j,1) - TTT_VEL(i,j,2) )                 &
                        / ( DZT(i,j,1) + DZT(i,j,2) )
    endif
    if ( DZT(i,j,km).ne.c0 ) then
      DVELDZ(i,j,km) = 2*( TTT_VEL(i,j,km-1) - TTT_VEL(i,j,km) )             &
                        / ( DZT(i,j,km-1) + DZT(i,j,km) )
    endif
  enddo !i 
enddo !j

do k = 2,km-1
  do j = 1,jmt
    do i = 1,imt
      if ( DZT(i,j,k).ne.c0 .and. DZT(i,j,k+1).ne.c0 ) then
        DVELDZ(i,j,k) = ( ( TTT_VEL(i,j,k-1) - TTT_VEL(i,j,k) )              &
                          / ( DZT(i,j,k-1) + DZT(i,j,k) )                    &
                        + ( TTT_VEL(i,j,k) - TTT_VEL(i,j,k+1) )              &
                          / (DZT(i,j,k) + DZT(i,j,k+1) ) )
      elseif ( DZT(i,j,k).ne.c0 .and. DZT(i,j,k+1).eq.c0 ) then
        DVELDZ(i,j,k) = 2*( TTT_VEL(i,j,k-1) - TTT_VEL(i,j,k) )              &
                          / ( DZT(i,j,k-1) + DZT(i,j,k) )
      endif
    enddo !i
  enddo !j
enddo !k

end subroutine nabla_hvel

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
real                                         :: A,B     ! to test 
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
      !T_FIELD(i,j,km) = p5 * (c3*W_FIELD(i,j,km) - W_FIELD(i,j,km-1) )
      T_FIELD(i,j,km) = p5 * W_FIELD(i,j,km)
    endif
  enddo !i
enddo !j

do k = 1,km
  call surf_int(1,1,imt,jmt,W_FIELD(:,:,:),TAREA,DZT(:,:,k),A)
  call surf_int(1,1,imt,jmt,T_FIELD(:,:,k),TAREA,DZT(:,:,k),B)
  write (*,*) k, ' wtt2ttt orig', A
  write (*,*) k, ' wtt2ttt new ', B
  write (*,*) k, ' wtt2ttt factor      ', B/A
enddo !k
end subroutine wtt2ttt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine interp_mom_fluxes(imt,jmt,k,DZT_k,DZU_k,DXT,DXU,DYT,DYU,&
                                   UE_FIELD,VN_FIELD,UE_NEW_FIELD,VN_NEW_FIELD)
!
!     interpolates momentum fluxes from east side of U-cell to momentum fluxes
!     at center of T-cells, for this a weighted average of the flux North/South
!     and East/West of the T-cell center point is used
!
implicit none

! input/output variables
integer,                              intent(in)  :: imt,jmt,k
real,             dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
double precision, dimension(imt,jmt), intent(in)  :: DXT,DXU,DYT,DYU
real,             dimension(imt,jmt), intent(in)  :: UE_FIELD ! (T)UU
real,             dimension(imt,jmt), intent(in)  :: VN_FIELD ! (T)UU
real,             dimension(imt,jmt), intent(out) :: UE_NEW_FIELD ! (T)TT
real,             dimension(imt,jmt), intent(out) :: VN_NEW_FIELD ! (T)TT
real :: A, B ! surface integrals, used for comparing old and new fields 

! southernmost gridpoints (j=1) on TTT-grid must be 0, 
! as there would not be UU gridpoints south of it otherwise

! east flux
do j = 2, jmt-1
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    if ( DZT_k(1,j-1).ne.c0 .or. DZT_k(1,j+1).ne.c0  ) then
      UE_NEW_FIELD(1,j) = 0.5 *                                            &
       ( UE_FIELD(imt,j  ) * DYU(imt,j  )                                  &
       + UE_FIELD(imt,j-1) * DYU(imt,j-1) )                                &
       / DYT(1,j)
    else
      UE_NEW_FIELD(1,j) = 0.0
    endif
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      if ( DZT_k(i,j-1).ne.c0 .or. DZT_k(i,j+1).ne.c0  ) then
        UE_NEW_FIELD(i,j) = 0.5 *                                           &
         ( UE_FIELD(i-1,j  ) * DYU(i-1,j  )                                 &
         + UE_FIELD(i-1,j-1) * DYU(i-1,j-1) )                               &
         / DYT(i,j)
      else
        UE_NEW_FIELD(i,j) = 0.0
      endif
    endif
  enddo !i
enddo !j

! north flux
do j = 2, jmt
  ! western edge
  if ( DZT_k(1,j).ne.c0 ) then
    if ( DZT_k(imt,j).ne.c0 .or. DZT_k(2,j).ne.c0  ) then
      VN_NEW_FIELD(1,j) = 0.5 *                                            &
       ( VN_FIELD(imt,j-1) * DXU(imt,j-1)                                  &
       + VN_FIELD(1  ,j-1) * DXU(1  ,j-1) )                                &
       / DXT(1,j)
    else
      VN_NEW_FIELD(1,j) = 0.0
    endif
  endif
  ! rest of array
  do i = 2, imt
    if ( DZT_k(i,j).ne.c0 ) then
      if ( DZT_k(i-1,j).ne.c0 .or. DZT_k(i+1,j).ne.c0  ) then
        VN_NEW_FIELD(i,j) = 0.5 *                                           &
         ( VN_FIELD(i-1,j-1) * DXU(i-1,j-1)                                 &
         + VN_FIELD(i  ,j-1) * DXU(i  ,j-1) )                               &
         / DXT(i,j)
      else
        VN_NEW_FIELD(i,j) = 0.0
      endif
    endif
  enddo !i
enddo !j

call surf_int(1,1,imt,jmt,UE_FIELD,    UAREA,DZT_k,A)
call surf_int(1,1,imt,jmt,UE_NEW_FIELD,TAREA,DZT_k,B)
!      write (*,*) 'A: ', A, 'B: ', B
write (*,*) k, 'rel. diff. of UE flux interp.: ', (A-B)/B
call surf_int(1,1,imt,jmt,VN_FIELD,    UAREA,DZT_k,A)
call surf_int(1,1,imt,jmt,VN_NEW_FIELD,TAREA,DZT_k,B)
!      write (*,*) 'A: ', A, 'B: ', B
write (*,*) k, 'rel. diff. of VN flux interp.: ', (A-B)/B

end subroutine interp_mom_fluxes

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

! southernmost gridpoints (j=1) on TTT-grid must be 0,
! as there would not be UU gridpoints south of it otherwise
do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
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
      NEW_FIELD_k(i,j) = 0.25                                               &
                   * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  ) * DZU_k(i  ,j  )   &
                     + FIELD_k(i-1,j  ) * UAREA(i-1,j  ) * DZU_k(i-1,j  )   &
                     + FIELD_k(i  ,j-1) * UAREA(i  ,j-1) * DZU_k(i  ,j-1)   &
                     + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) * DZU_k(i-1,j-1) ) &
                   / TAREA(i,j) / DZT_k(i,j)
    endif
  enddo !i
enddo !j

end subroutine uu2tt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,TUU_FIELD,TTT_FIELD)
implicit none
!
!  interpolates a 3D field from UU to TT
!

integer                                  :: imt,jmt,km
real, dimension(imt,jmt)                 :: WORK
double precision, dimension(imt,jmt),    intent(in)  :: TAREA, UAREA
real, dimension(imt,jmt,km), intent(in)  :: TUU_FIELD, DZU, DZT
real, dimension(imt,jmt,km), intent(out) :: TTT_FIELD

do k = 1,km
  call uu2tt(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,TUU_FIELD(:,:,k),WORK)
  TTT_FIELD(:,:,k) = WORK
enddo

end subroutine uu2tt_3D

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

! southernmost gridpoints (j=1) on TTT-grid must be 0,
! as there would not be UU gridpoints south of it otherwise
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

subroutine uu2tt_scalar_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,TUU_FIELD,TTT_FIELD)
implicit none
!
!  interpolates a 3D field from UU to TT
!

integer                                  :: imt,jmt,km
real, dimension(imt,jmt)                 :: WORK
double precision, dimension(imt,jmt),    intent(in)  :: TAREA, UAREA
real, dimension(imt,jmt,km), intent(in)  :: TUU_FIELD, DZU, DZT
real, dimension(imt,jmt,km), intent(out) :: TTT_FIELD

do k = 1,km
  call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,TUU_FIELD(:,:,k),WORK)
  TTT_FIELD(:,:,k) = WORK
enddo

end subroutine uu2tt_scalar_3D

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!
!  Below are some useful (model-native?) subroutines
!

!***********************************************************************

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

!***********************************************************************

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

!***********************************************************************

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

!***********************************************************************

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

!***********************************************************************


!***********************************************************************
!BOP
! !IROUTINE: state
! !INTERFACE:

subroutine state(S_orig, T, P_orig, RHO)

!  ADJUSTED/SLIMMED DOWN FROM state_mod.F90 June 2016


! !DESCRIPTION:
!  Returns the density of water at level k from equation of state
!  $\rho = \rho(d,\theta,S)$ where $d$ is depth, $\theta$ is
!  potential temperature, and $S$ is salinity. the density can be
!  returned as a perturbation (RHOOUT) or as the full density
!  (RHOFULL). Note that only the polynomial EOS choice will return
!  a perturbation density; in other cases the full density is returned
!  regardless of which argument is requested.
!
!  This routine also computes derivatives of density with respect
!  to temperature and salinity at level k from equation of state
!  if requested (ie the optional arguments are present).
!
!  If $k = kk$ are equal the density for level k is returned.
!  If $k \neq kk$ the density returned is that for a parcel
!  adiabatically displaced from level k to level kk.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:


double precision, intent(in) :: & 
  T,               &! pot. temperature at level k  [degC]
  S_orig,          &! salinity at level k         [kg/kg]
  P_orig            ! pressure at reference level   [bar]


! !OUTPUT PARAMETERS:

double precision, intent(out) :: &
  RHO               ! density of water [kg m^-3]

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

double precision :: &
  S,                &! [PSU]
  P,                &! [dbar]
  SQR,DENOMK,       &! work arrays
  WORK1, WORK2, T2

double precision, parameter :: &
  c10   = 1.0E01,              &
  c1000 = 1.0E03


!-----------------------------------------------------------------------
!
!  first check for valid range if requested
!
!-----------------------------------------------------------------------

! THIS IS SKIPPED HERE

!-----------------------------------------------------------------------
!
!  MWJF EOS coefficients
!
!-----------------------------------------------------------------------

!*** these constants will be used to construct the numerator
   
double precision, parameter ::   &
  mwjfnp0s0t0 =   9.99843699e+2, &
  mwjfnp0s0t1 =   7.35212840e+0, &
  mwjfnp0s0t2 =  -5.45928211e-2, &
  mwjfnp0s0t3 =   3.98476704e-4, &
  mwjfnp0s1t0 =   2.96938239e+0, &
  mwjfnp0s1t1 =  -7.23268813e-3, &
  mwjfnp0s2t0 =   2.12382341e-3, &
  mwjfnp1s0t0 =   1.04004591e-2, &
  mwjfnp1s0t2 =   1.03970529e-7, &
  mwjfnp1s1t0 =   5.18761880e-6, &
  mwjfnp2s0t0 =  -3.24041825e-8, &
  mwjfnp2s0t2 =  -1.23869360e-11

!*** these constants will be used to construct the denominator

double precision, parameter ::    &
  mwjfdp0s0t0 =   1.0e+0,         &
  mwjfdp0s0t1 =   7.28606739e-3,  &
  mwjfdp0s0t2 =  -4.60835542e-5,  &
  mwjfdp0s0t3 =   3.68390573e-7,  &
  mwjfdp0s0t4 =   1.80809186e-10, &
  mwjfdp0s1t0 =   2.14691708e-3,  &
  mwjfdp0s1t1 =  -9.27062484e-6,  &
  mwjfdp0s1t3 =  -1.78343643e-10, &
  mwjfdp0sqt0 =   4.76534122e-6,  &
  mwjfdp0sqt2 =   1.63410736e-9,  &
  mwjfdp1s0t0 =   5.30848875e-6,  &
  mwjfdp2s0t3 =  -3.03175128e-16, &
  mwjfdp3s0t1 =  -1.27934137e-17

!*** MWJF numerator coefficients including pressure

double precision ::                                                     &
  mwjfnums0t0, mwjfnums0t1, mwjfnums0t2, mwjfnums0t3,                   &
  mwjfnums1t0, mwjfnums1t1, mwjfnums2t0,                                &
  mwjfdens0t0, mwjfdens0t1, mwjfdens0t2, mwjfdens0t3, mwjfdens0t4,      &
  mwjfdens1t0, mwjfdens1t1, mwjfdens1t3,                                &
  mwjfdensqt0, mwjfdensqt2
!-----------------------------------------------------------------------

!  McDougall, Wright, Jackett, and Feistel EOS
!  test value : rho = 1.033213242 for
!  S = 35.0 PSU, theta = 20.0, pressz = 200.0
!
!  from paper:
!  [ (PSU, degC, dbar) -> kg m^-3 ]
!  (35, 25, 2000) -> 1031.654 229
!  ( 0, 20, 1000) -> 1017.726 743
!  (40, 12, 8000) -> 1062.928 258
!
!-----------------------------------------------------------------------

!   unit conversion
P   = c10*P_orig   ! [bar] -> [dbar]
S   = c1000*S_orig ! [kg/kg] -> [PSU]
SQR = sqrt(S)      ! square root 

!      write (*,*) 'S =', S, 'PSU'
!      write (*,*) 'T =', T, 'degC'
!      write (*,*) 'P =', P, 'dbar'

!***
!*** first calculate numerator of MWJF density [P_1(S,T,p)]
!***

mwjfnums0t0 = mwjfnp0s0t0 + P*(mwjfnp1s0t0 + P*mwjfnp2s0t0)
mwjfnums0t1 = mwjfnp0s0t1 
mwjfnums0t2 = mwjfnp0s0t2 + P*(mwjfnp1s0t2 + P*mwjfnp2s0t2)
mwjfnums0t3 = mwjfnp0s0t3
mwjfnums1t0 = mwjfnp0s1t0 + P*mwjfnp1s1t0
mwjfnums1t1 = mwjfnp0s1t1
mwjfnums2t0 = mwjfnp0s2t0

WORK1 = mwjfnums0t0 + T * (mwjfnums0t1 + T * (mwjfnums0t2 +              &
        mwjfnums0t3 * T)) + S * (mwjfnums1t0 +                           &
        mwjfnums1t1 * T + mwjfnums2t0 * S)

!***
!*** now calculate denominator of MWJF density [P_2(S,T,p)]
!***

mwjfdens0t0 = mwjfdp0s0t0 + P*mwjfdp1s0t0
mwjfdens0t1 = mwjfdp0s0t1 + P**3 * mwjfdp3s0t1
mwjfdens0t2 = mwjfdp0s0t2
mwjfdens0t3 = mwjfdp0s0t3 + P**2 * mwjfdp2s0t3
mwjfdens0t4 = mwjfdp0s0t4
mwjfdens1t0 = mwjfdp0s1t0
mwjfdens1t1 = mwjfdp0s1t1
mwjfdens1t3 = mwjfdp0s1t3
mwjfdensqt0 = mwjfdp0sqt0
mwjfdensqt2 = mwjfdp0sqt2

WORK2 = mwjfdens0t0 + T * (mwjfdens0t1 + T * (mwjfdens0t2 +              &
        T   * (mwjfdens0t3 + mwjfdens0t4 * T))) +                           &
        S   * (mwjfdens1t0 + T * (mwjfdens1t1 + T*T*mwjfdens1t3)+           &
        SQR * (mwjfdensqt0 + T*T*mwjfdensqt2))

RHO   = WORK1/WORK2

! SKIPPED RHOFULL, DRHODT, DRHODS CALCULATIONS HERE
! SKIPPED OTHER EOS FUNCTIONS

end subroutine state


end program
