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
nrec_PDW, nrec_UPD, nrec_VPD, nrec_WPD, nrec_T_STRONG_REST,nrec_S_STRONG_REST, &
nrec_S_WEAK_REST
! obsolete: RHO2
character*120 ::                                                            &
  LEC_folder,geometry1_file,geometry2_file,ref_state_file,projects_folder,  &
  inputfile,bin_file,global_file,levels_file
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
cPem_sint,    cKem_sint,    cPKm_sint,    cPKe_sint,                        &
tdz,p_z,vol

! internal model variables 
integer :: i, j, k

integer, dimension(:,:), allocatable :: kmT

real, dimension(:,:), allocatable ::                                        &
ALPHA01, SHF, BETA01, SFWF, STFRHO, SSFRHO, TAUX, TAUY, UTAUX, VTAUY,       &
RHOW_k, RHOU_k, RHOV_k,                                                     &
gPmh, gPms, gPm, gPeh, gPes, gPe, gKm, gKe, TT_gKe, TT_gKm,                 &
EVAP, PRECIP, SWREST, SSREST, RUNOFF, EVAPRHO, PRECIPRHO, SWRESTRHO,        &
SSRESTRHO, RUNOFFRHO, SWNET, LWNET, LATENT, SENSIBLE, TSREST,               &
gPt, gKt, TT_gKt, rKt_k,                                                    &
UEU_k, UEV_k, VNU_k, VNV_k, WTU_k, WTV_k, UUE_k, VUN_k, UTE_k, VTN_k,       &
UEU_eddy_k, UEV_eddy_k, VNU_eddy_k, VNV_eddy_k,                             &
TTT_UEU_k, TTT_UEV_k, TTT_VNU_k, TTT_VNV_k, TTT_WTU_k, TTT_WTV_k,           &
TTT_UEU_eddy_k, TTT_UEV_eddy_k, TTT_VNU_eddy_k, TTT_VNV_eddy_k,             &
UPD_k, VPD_k, TTT_UPD_k, TTT_VPD_k,                                         &
UPD_eddy_k, VPD_eddy_k, TTT_UPD_eddy_k, TTT_VPD_eddy_k,                     &
WUU_WTU_eddy_k,  WUU_WTV_eddy_k, WTT_WTU_eddy_k, WTT_WTV_eddy_k,            &
geometry2,ref_state

real, dimension(:,:,:), allocatable ::                                      &
UVEL, VVEL, KE, RHO, Q, RHOU, RHOV, WVEL, UVEL2, VVEL2, UV, UW, VW,         &
RHOW, PD, PD2, PDU, PDV, rKm, rKe, rPm, rPe, cPem, cKem, cPKm, cPKe, SALT,  &
TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV, TTT_UW, TTT_VW,           &
DRHODX, DRHODY, DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DPDDX, DPDDY,           &
TTT_WVEL, TTT_RHOW, TTT_PDW, TEMP, PDW,                                     &
WUU_WVEL, WTT_WTU_eddy, WTT_WTV_eddy, TTT_WTU_eddy, TTT_WTV_eddy

real, dimension(:), allocatable ::                                          &
z1, z2, rho_ref, salt_ref, temp_ref,                                        &
pressz, rho_new_ref, pd_new_ref, n0_new

real, dimension(:),     allocatable ::                                      &
dz,area,tdepth,n0,n0_inv,pd_ref

real, dimension(:,:),   allocatable ::                                      &
 DXT, DYT, TAREA, DXU, DYU, UAREA

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

read  (*,*)        ntavg
read  (*,'(a3)')   year

imt             = 3600
jmt             = 2400
km              =   42

LEC_folder      = '/home/dijkbio/andre/LEC/'
projects_folder = '/projects/0/samoc/jan/Andree/'
geometry1_file  = trim(LEC_folder)//'input/geometry1'
geometry2_file  = trim(LEC_folder)//'input/geometry2'
ref_state_file  = trim(LEC_folder)//'input/ref_state'

if ( ntavg==1 ) then
  inputfile   = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.year.'//year
  bin_file    = trim(projects_folder)//'LEC_bin_1_'//year
  levels_file = trim(LEC_folder)//'results/levels/level_LEC_1_'//year//'.out'
  global_file = trim(LEC_folder)//'results/global/global_LEC_1_'//year//'.out'
elseif ( ntavg==5 ) then
  inputfile   = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.5.year.'//year
  bin_file    = trim(projects_folder)//'LEC_bin_5_'//year
  levels_file = trim(LEC_folder)//'results/levels/level_LEC_5_'//year//'.out'
  global_file = trim(LEC_folder)//'results/global/global_LEC_5_'//year//'.out'
elseif ( ntavg==11 ) then
  inputfile   = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.11.year.'//year
  bin_file    = trim(projects_folder)//'LEC_bin_11_'//year
  levels_file = trim(LEC_folder)//'results/levels/level_LEC_11_'//year//'.out'
  global_file = trim(LEC_folder)//'results/global/global_LEC_11_'//year//'.out'
elseif ( ntavg==51 ) then
  inputfile   = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.51.year.'//year
  bin_file    = trim(projects_folder)//'LEC_bin_51_'//year
  levels_file = trim(LEC_folder)//'results/levels/level_LEC_51_'//year//'.out'
  global_file = trim(LEC_folder)//'results/global/global_LEC_51_'//year//'.out'
endif

write(*,*) inputfile

! (time averaged) input file
open(1,file=inputfile,  form='unformatted',access='direct',recl=imt*jmt,       &
        status='old')
! output files
open(3,file=bin_file,   form='unformatted',access='direct',recl=imt*jmt)
open(4,file=levels_file,form='formatted')
open(5,file=global_file,form='formatted')

write(*,*) 'input file:  ', inputfile
write(*,*) 'output file: ', bin_file

nrec_UVEL          =    1
nrec_VVEL          =   43
nrec_UVEL2         =   85
nrec_VVEL2         =  127
nrec_STFRHO        =  169
nrec_SSFRHO        =  170
nrec_EVAPRHO       =  171
nrec_PRECIPRHO     =  172
nrec_SWRESTRHO     =  173
nrec_SSRESTRHO     =  174
nrec_RUNOFFRHO     =  175
nrec_UTAUX         =  176
nrec_VTAUY         =  177
nrec_ALPHA01       =  178
nrec_BETA01        =  179
nrec_KE            =  180
nrec_TEMP          =  222
nrec_SALT          =  264
nrec_RHO           =  306
nrec_UV            =  348
nrec_SSH           =  390
nrec_SWNET         =  391
nrec_LWNET         =  392
nrec_LATENT        =  393
nrec_SENSIBLE      =  394
nrec_T_STRONG_REST =  395
nrec_EVAP          =  396
nrec_PRECIP        =  397
nrec_S_WEAK_REST   =  398
nrec_S_STRONG_REST =  399
nrec_RUNOFF        =  400
nrec_SHF           =  401
nrec_SFWF          =  402
nrec_TAUX          =  403
nrec_TAUY          =  404
nrec_WVEL          =  405
nrec_UEU           =  447
nrec_VNU           =  489
nrec_UEV           =  531
nrec_VNV           =  573
nrec_WTU           =  615
nrec_WTV           =  657
nrec_UET           =  699
nrec_VNT           =  741
nrec_WTT           =  783
nrec_UES           =  825
nrec_VNS           =  867
nrec_WTS           =  909
nrec_ADVT          =  951
nrec_ADVS          =  952
nrec_Q             =  953
nrec_PD            =  995
nrec_PD2           = 1037
nrec_RHO2          = 1079
nrec_UW            = 1121
nrec_VW            = 1163
nrec_RHOU          = 1205
nrec_RHOV          = 1247
nrec_RHOW          = 1289
nrec_URHO          = 1331
nrec_VRHO          = 1373
nrec_WRHO          = 1415
nrec_PDW           = 1457
nrec_PDU           = 1499
nrec_PDV           = 1541
nrec_UPD           = 1583
nrec_VPD           = 1625
nrec_WPD           = 1667

!===============================================================================
!  geometry
!===============================================================================

! read 2D geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km) )
open(2,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,    &
       status='old')
read(2,rec=1) DXT ! [m]
read(2,rec=2) DXU ! [m]
read(2,rec=3) DYT ! [m]
read(2,rec=4) DYU ! [m]
read(2,rec=5) TAREA ! [m^2]
read(2,rec=6) UAREA ! [m^2]
do k=1,km
  read(2,rec=6+   k) DZT(:,:,k) ! [m]
  read(2,rec=6+km+k) DZU(:,:,k) ! [m]
enddo
close(2)

! read 1D geometry fields
! all real: k,dz[m],tdepth[m],area[m^2],p[bar]
allocate( geometry2(6,km),tdz(km),dz(km),tdepth(km),area(km),p_z(km),vol(km) )
open(2,file=geometry2_file,access='direct',form='unformatted',recl=6,     &
       status='old')
do k=1,km
  read(2,rec=k) geometry2(:,k) 
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
p_z(:)    = geometry2(5,:)
vol(:)    = geometry2(6,:)
close(2)

! read referance state fields
! k, dz, tdepth, area, p,
! <T>, <S>, <RHO>, <PD>, <D(D,T,0)>, <D(D,T,p)>,
! <Q>, dz(<RHO>), dz(<PD>), dz(<D(S,T,0)>), dz(<D(S,T,p)>)
! <Q_w>, dz(<RHO>_w), dz(<PD>_w), dz(<D(S,T,0)>_w), dz(<D(S,T,p)>_w)

allocate( ref_state(21,km), rho_ref(km), pd_ref(km), n0(km), n0_inv(km) )
open(2,file=ref_state_file,access='direct',form='unformatted',recl=21,     &
       status='old')
do k=1,km
  read(2,rec=k) ref_state(:,k) ! all real
enddo
rho_ref(:) = ref_state( 8,:) ![kg/m^3]; <RHO>
pd_ref(:)  = ref_state( 9,:) ![kg/m^3]; <PD>
n0(:)      = ref_state(12,:) ![kg/m^4]; <Q>
n0(1)      = 2.0*n0(1) ! at the surface layer Q is wrongly calculated in POP
do k=1,km
  n0_inv(k) = 1.0/n0(k) ![m^4/kg]
enddo
close(2)
!===============================================================================
!  read fields
!===============================================================================

! 3D
allocate( UVEL(imt,jmt,km), VVEL(imt,jmt,km), PD(imt,jmt,km), TEMP(imt,jmt,km) )
call load_3D_field(1,nrec_PD,  PD  ) ! [g/cm^3]; (actual pot. density) - rho0 
call load_3D_field(1,nrec_UVEL,UVEL) ! [cm/s]
call load_3D_field(1,nrec_VVEL,VVEL) ! [cm/s]
call load_3D_field(1,nrec_TEMP,TEMP) ! [degC]

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
call load_3D_field(1,nrec_PD2, PD2) ! [g^2/cm^6]
call load_3D_field(1,nrec_KE,  KE ) ! [cm^2/s^2]

allocate(                                                                      &
  rKm(imt,jmt,km), rKe(imt,jmt,km), rPm(imt,jmt,km), rPe(imt,jmt,km),          &
  rKm_sint(km),    rKe_sint(km),    rPm_sint(km),    rPe_sint(km) )

do k = 1, km
  rPm(:,:,k) = -p5*g*n0_inv(k) * ( PD(:,:,k)*1.0E03  -  pd_ref(k) )**2 
  rPe(:,:,k) = -p5*g*n0_inv(k) * (   PD2(:,:,k)      -  PD(:,:,k)**2 )* 1.0E06
  rKm(:,:,k) =  p5*rho0 * ( UVEL(:,:,k)**2 + VVEL(:,:,k)**2 ) * 1.0E-04 
  rKe(:,:,k) =  rho0 * KE(:,:,k) * 1.0E-04 - rKm(:,:,k)
enddo !k

! integrate
call vol_int_full(1,1,1,imt,jmt,km,rPm,TAREA,dz,DZT,rPm_int)
call vol_int_full(1,1,1,imt,jmt,km,rPe,TAREA,dz,DZT,rPe_int)
call vol_int_part(1,1,1,imt,jmt,km,rKm,TAREA,dz,DZU,rKm_int)
call vol_int_part(1,1,1,imt,jmt,km,rKe,TAREA,dz,DZU,rKe_int)

! surface integrals
do k = 1,km
  call surf_int_3D(1,1,imt,jmt,rPm(:,:,k),TAREA,DZT(:,:,k),rPm_sint(k))
  call surf_int_3D(1,1,imt,jmt,rPe(:,:,k),TAREA,DZT(:,:,k),rPe_sint(k))
  call surf_int_3D(1,1,imt,jmt,rKm(:,:,k),TAREA,DZU(:,:,k),rKm_sint(k))
  call surf_int_3D(1,1,imt,jmt,rKe(:,:,k),TAREA,DZU(:,:,k),rKe_sint(k))
enddo !k

! output
300 FORMAT(A10, ES14.6E2, A2, ES10.2E2)
301 FORMAT(A10, ES14.6E2, A2, ES10.2E2, A15, ES14.6E2, A2, ES10.2E2)
write (*,300) 'rPm     =',    rPm_int, 'J',    rPm_int/1.66E23
write (*,300) 'rPe     =',    rPe_int, 'J',    rPe_int/6.38E18
write (*,300) 'rKm     =',    rKm_int, 'J',    rKm_int/1.27E18
write (*,300) 'rKe     =',    rKe_int, 'J',    rKe_int/3.55E18

do i = 1,4
  do k = 1,km
    if ( i==1 ) write (3,rec=(i-1)*km+k) rPm(:,:,k) 
    if ( i==2 ) write (3,rec=(i-1)*km+k) rPe(:,:,k) 
    if ( i==3 ) write (3,rec=(i-1)*km+k) rKm(:,:,k) 
    if ( i==4 ) write (3,rec=(i-1)*km+k) rKe(:,:,k) 
  enddo !k
enddo !i

! rKm not deallocated
deallocate( rKm, rKe, rPm, rPe, PD2, KE)

!===============================================================================
!  2. GENERATION
!===============================================================================

write (*,*) ''
write (*,*) '2. GENERATION'
!  all PE terms quantities are on TT-grid, all KE terms on UU-grid

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

allocate(                                                                      &
  gPmh(imt,jmt), gPms(imt,jmt), gPm(imt,jmt), gPeh(imt,jmt),                   &
  gPes(imt,jmt), gPe(imt,jmt),  gKm(imt,jmt), gKe(imt,jmt),                    &
  gPt(imt,jmt), gKt(imt,jmt) )

!  2.1/2.2 Mean/Eddy Potential Energy Generation
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

!  2.3/2.4 Mean/Eddy Kinetic Energy Generation
gKm  = ( UVEL(:,:,1)*TAUX + VVEL(:,:,1)*TAUY )*1.0E-03
gKe  = ( UTAUX + VTAUY )*1.0E-03 - gKm

call surf_int_2D(1,1,imt,jmt,gPmh,TAREA,DZT(:,:,1),gPmh_int)
call surf_int_2D(1,1,imt,jmt,gPms,TAREA,DZT(:,:,1),gPms_int)
call surf_int_2D(1,1,imt,jmt, gPm,TAREA,DZT(:,:,1), gPm_int)
call surf_int_2D(1,1,imt,jmt,gPeh,TAREA,DZT(:,:,1),gPeh_int)
call surf_int_2D(1,1,imt,jmt,gPes,TAREA,DZT(:,:,1),gPes_int)
call surf_int_2D(1,1,imt,jmt, gPe,TAREA,DZT(:,:,1), gPe_int)
call surf_int_2D(1,1,imt,jmt, gKm,TAREA,DZU(:,:,1), gKm_int)
call surf_int_2D(1,1,imt,jmt, gKe,TAREA,DZU(:,:,1), gKe_int)

write (3,rec=169) gPm(:,:) 
write (3,rec=170) gPe(:,:) 
write (3,rec=171) gKm(:,:) 
write (3,rec=172) gKe(:,:) 
write (3,rec=173) gPmh(:,:) 
write (3,rec=174) gPms(:,:) 
write (3,rec=175) gPeh(:,:) 
write (3,rec=176) gPes(:,:) 


deallocate( gPmh, gPms, gPm, gPeh, gPes, gPe, gKm, gKe,                        &
ALPHA01, SHF, BETA01, SFWF, STFRHO, SSFRHO, TAUX, TAUY, UTAUX, VTAUY, gPt, gKt )

!  output
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

call load_3D_field(1,nrec_WVEL, WVEL ) ! [cm/s]
call load_3D_field(1,nrec_UVEL2,UVEL2) ! [cm^2/s^2]
call load_3D_field(1,nrec_VVEL2,VVEL2) ! [cm^2/s^2]
call load_3D_field(1,nrec_UV,   UV   ) ! [cm^2/s^2]
call load_3D_field(1,nrec_UW,   UW   ) ! [cm^2/s^2]
call load_3D_field(1,nrec_VW,   VW   ) ! [cm^2/s^2]
call load_3D_field(1,nrec_PDU,  PDU  ) ! [g/cm^2/s]
call load_3D_field(1,nrec_PDV,  PDV  ) ! [g/cm^2/s]

! calculate non-zonality parameter
! VVEL/UVEL integrated over [0E,40E]x[90S,45S]=[1100,1500]x[0,676]
!call vol_int(1,1100,1,676,1500,17,VVEL2/UVEL2,UAREA,DZU,V2U2    )
!call vol_int(1,1100,1,676,1500,17,DZU/DZU    ,UAREA,DZU,V2U2_vol)
!V2U2 = V2U2/V2U2_vol
!write (*,*) 'V2U2:', V2U2

!  first create new fields:
!  3.1: (tuu->ttt) uvel,vvel;  (quadratic? central difference) drho/dx, drho/dy,
!       (nabla, horizontally linear central difference) 
!          du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz
!  3.2: (tuu->ttt) uvel2,vvel2,uvelvel
!  3.3: (wtt->ttt) wvel
!  4.4: (wtt->ttt) wvelrho

! Interpolation
! > TUU->TTT: TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV
allocate(                                                                      &
  TTT_UVEL(imt,jmt,km), TTT_VVEL(imt,jmt,km), TTT_UVEL2(imt,jmt,km),           &
  TTT_VVEL2(imt,jmt,km), TTT_UV(imt,jmt,km) )
call uu2tt_3D(DZT,DZU,TAREA,UAREA,UVEL ,TTT_UVEL )
call uu2tt_3D(DZT,DZU,TAREA,UAREA,VVEL ,TTT_VVEL )
call uu2tt_3D(DZT,DZU,TAREA,UAREA,UVEL2,TTT_UVEL2)
call uu2tt_3D(DZT,DZU,TAREA,UAREA,VVEL2,TTT_VVEL2)
call uu2tt_3D(DZT,DZU,TAREA,UAREA,UV   ,TTT_UV   )
deallocate( UVEL2, VVEL2, UV )

! > WTT->TTT: TTT_UW, TTT_VW, TTT_WVEL, TTT_PDW
allocate( TTT_UW(imt,jmt,km), TTT_VW(imt,jmt,km), TTT_WVEL(imt,jmt,km) )
call wtt2ttt(  UW,DZT,  TTT_UW)
call wtt2ttt(  VW,DZT,  TTT_VW)
call wtt2ttt(WVEL,DZT,TTT_WVEL)

deallocate( UW, VW )

! Derivatives
! > DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ
allocate( DUDX(imt,jmt,km), DUDY(imt,jmt,km), DUDZ(imt,jmt,km),                &
          DVDX(imt,jmt,km), DVDY(imt,jmt,km), DVDZ(imt,jmt,km) )
call nabla_hvel(UVEL,TTT_UVEL,dz,DXT,DYT,DXU,DYU,DZT,DZU,TAREA,     &
                DUDX,DUDY,DUDZ)
call nabla_hvel(VVEL,TTT_VVEL,dz,DXT,DYT,DXU,DYU,DZT,DZU,TAREA,     &
                DVDX,DVDY,DVDZ)
! > DRHODX, DRHODY
allocate( DPDDX(imt,jmt,km), DPDDY(imt,jmt,km) )
!call gradient(DXU,DYU,PD,DPDDX,DPDDY)
call grad_rho(DZT,DZU,DXU,DYU,TAREA,UAREA,PD,DPDDX,DPDDY)
!===============================================================================

allocate(                                                                      &
  cPem(imt,jmt,km), cKem(imt,jmt,km), cPKm(imt,jmt,km), cPKe(imt,jmt,km),      &
  cPem_sint(km),    cKem_sint(km),    cPKm_sint(km),    cPKe_sint(km) )

!  3.1/3.2 Eddy to Mean Energy
do k = 1, km
  cPem(:,:,k) = - g * n0_inv(k) * 1.0E04                                       &
              * ( ( PDU(:,:,k) - PD(:,:,k) * TTT_UVEL(:,:,k)) * DPDDX(:,:,k)   &
              +   ( PDV(:,:,k) - PD(:,:,k) * TTT_VVEL(:,:,k)) * DPDDY(:,:,k) )
  cKem(:,:,k) = rho0 * 1.0E-06                                                 &
              * ( ( TTT_UVEL2(:,:,k) - TTT_UVEL(:,:,k)**2 ) * DUDX(:,:,k)      &
              +   ( TTT_VVEL2(:,:,k) - TTT_VVEL(:,:,k)**2 ) * DVDY(:,:,k)      &
              +   ( TTT_UV(:,:,k)    - TTT_UVEL(:,:,k) * TTT_VVEL(:,:,k) )     &
              *                              ( DUDY(:,:,k) + DVDX(:,:,k) )     &
              +   ( TTT_UW(:,:,k)    - TTT_UVEL(:,:,k) * TTT_WVEL(:,:,k) )     &
              *                                              DUDZ(:,:,k)       &
              +   ( TTT_VW(:,:,k)    - TTT_VVEL(:,:,k) * TTT_WVEL(:,:,k) )     &
              *                                              DVDZ(:,:,k)   )
enddo !k

deallocate( PDU, PDV, TTT_UVEL2, TTT_VVEL2, TTT_UV, TTT_UW, TTT_VW, TTT_WVEL )

!  3.3/3.4 Potential to Kinetic Energy
allocate( PDW(imt,jmt,km), TTT_PDW(imt,jmt,km) )
call load_3D_field(1,nrec_PDW,PDW)
call wtt2ttt( PDW,DZT, TTT_PDW)

do k = 1, km
  if (k==1) then
    cPKm(:,:,k) = ( WVEL(:,:,k+1)*(PD(:,:,k)+PD(:,:,k+1)) )
  elseif (k==km) then
    cPKm(:,:,k) = ( WVEL(:,:,k  )*(PD(:,:,k)+PD(:,:,k-1)) )
  else
    cPKm(:,:,k) = ( WVEL(:,:,k+1)*(PD(:,:,k)+PD(:,:,k+1))           &
                  + WVEL(:,:,k  )*(PD(:,:,k)+PD(:,:,k-1)) )
  endif
  cPKm(:,:,k) = -g * p25 * cPKm(:,:,k) * 1.0E01
  cPKe(:,:,k) = -g * TTT_PDW(:,:,k)*1.0E01 - cPKm(:,:,k)
enddo !k

deallocate( PDW, TTT_PDW )

!  volume integrals
call vol_int_full(1,1,1,imt,jmt,km,cPem,TAREA,dz,DZT,cPem_int)
call vol_int_part(1,1,1,imt,jmt,km,cKem,TAREA,dz,DZT,cKem_int)
call vol_int_full(1,1,1,imt,jmt,km,cPKm,TAREA,dz,DZT,cPKm_int)
call vol_int_full(1,1,1,imt,jmt,km,cPKe,TAREA,dz,DZT,cPKe_int)

!  surface integrals
do k = 1,km
  call surf_int_3D(1,1,imt,jmt,cPem(:,:,k),TAREA,DZT(:,:,k),cPem_sint(k))
  call surf_int_3D(1,1,imt,jmt,cKem(:,:,k),TAREA,DZT(:,:,k),cKem_sint(k))
  call surf_int_3D(1,1,imt,jmt,cPKm(:,:,k),TAREA,DZT(:,:,k),cPKm_sint(k))
  call surf_int_3D(1,1,imt,jmt,cPKe(:,:,k),TAREA,DZT(:,:,k),cPKe_sint(k))
enddo !k

!  output
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

!  =============================================================================
!  4. DISSIPATION
!  =============================================================================

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


!===============================================================================
!   FILE OUTPUT
!===============================================================================

! closing the output binary file that includes all fields
close(3)

write (*,*) ''
write (*,*) 'OUTPUT'

! global surface integrals
! (for reservoirs/exchange/dissipation terms)
200 FORMAT (13(A,","),A)
201 FORMAT (13(E15.7,","),E15.7)
write (4,200) 'k','depth_k','area_k','rho_ref_k','pd_ref_k','n0_k',            &
              'rPm_k','rPe_k','rKm_k','rKe_k',                                 &
              'cPem_k','cKem_k','cPKm_k','cPKe_k'
do k = 1,km
  !surface integrals of energy budget terms
  write (4,201) real(k), tdepth(k), area(k), rho_ref(k), pd_ref(k), n0(k),     &
                rPm_sint(k),  rPe_sint(k),  rKm_sint(k),  rKe_sint(k),         &
                cPem_sint(k), cKem_sint(k), cPKm_sint(k), cPKe_sint(k)
enddo
close (4)

! global integrals: 
! volume (for reservoirs/exchange/dissipation terms)
! surface (for generation terms)
100 FORMAT (20(A,","),A)
101 FORMAT (20(E15.7,","),E15.7)
write (5,100) 'pd_avg',                                                        &
              'rPm','rPe','rKm','rKe',                                         &
              'gPmh','gPms','gPm',                                             &
              'gPeh','gPes','gPe',                                             &
              'gKm','gKe',                                                     &
              'cPem','cKem','cPKm','cPKe',                                     &
              'dPm','dPe','dKm','dKe'
write (5,101) pd_avg,                                                          &
              rPm_int,     rPe_int,     rKm_int,     rKe_int,                  &
              gPmh_int,    gPms_int,    gPm_int,                               &
              gPeh_int,    gPes_int,    gPe_int,                               &
              gKm_int,     gKe_int,                                            &
              cPem_int,    cKem_int,    cPKm_int,    cPKe_int,                 &
              dPm,          dPe,        dKm,         dKe
close (5)


write (*,*) '>>> done writing output'

end program
