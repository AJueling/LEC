program ref_state
implicit none

!===============================================================================
!
!  this script calculates the reference state variable for the last 51 years of
!  the extravars_viebahn high-resolution POP run
!  
!===============================================================================

!===============================================================================
!  variables
!===============================================================================


character*120 ::                                                               &
  input_folder,in_depths,geometry1_file,geometry2_file,ref_state_file,tavg_file

integer                                         ::                             &
  imt,jmt,km,rec_length,i,j,k,                                                 &
  nrec_TEMP,nrec_SALT,nrec_RHO,nrec_Q,nrec_PD

double precision :: D_test
real,             dimension(:),     allocatable ::                             &
  dz,tdepth,area,p_z,vol,                                                      &
  T_avg,S_avg,RHO_avg,PD_avg,Q_avg,D0_avg,Dp_avg,PD_ddz,RHO_ddz,D0_ddz,Dp_ddz, &
  T_wavg,S_wavg,RHO_wavg,PD_wavg,Q_wavg,D0_wavg,Dp_wavg,PD_wddz,RHO_wddz,      &
  D0_wddz,Dp_wddz
real,             dimension(:,:),   allocatable ::                             &
  DXT, DYT, TAREA, DXU, DYU, UAREA, geometry2
real,             dimension(:,:,:), allocatable ::                             &
  RHO,Q,PD,TEMP,SALT,DZT,DZU,D0,Dp
double precision, dimension(:),     allocatable ::                             &
  tdz
double precision, parameter                     ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,rho0 = 4.1/3.996*1000

imt            = 3600
jmt            = 2400
km             =   42

input_folder   = '/home/dijkbio/andre/LEC/input/'
in_depths      = trim(input_folder)//'in_depths.42.dat'
geometry1_file = trim(input_folder)//'geometry1'
geometry2_file = trim(input_folder)//'geometry2'
ref_state_file = trim(input_folder)//'ref_state'

tavg_file      = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.51.year.301'

nrec_TEMP      = 222
nrec_SALT      = 264
nrec_RHO       = 306
nrec_Q         = 953
nrec_PD        = 995

write (*,*) ''
write (*,*) '--- REFERENCE STATE ---'
write (*,*) ''

! read 2D geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km) )
open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,     &
       status='old')
read(1,rec=1) DXT ! [m]
read(1,rec=2) DXU ! [m]
read(1,rec=3) DYT ! [m]
read(1,rec=4) DYU ! [m]
read(1,rec=5) TAREA ! [m^2]
read(1,rec=6) UAREA ! [m^2]
do k=1,km
  read(1,rec=6+   k) DZT(:,:,k) ! [m]
  read(1,rec=6+km+k) DZU(:,:,k) ! [m]
enddo
close(1)

! read 1D geometry fields
allocate( geometry2(6,km),tdz(km),dz(km),tdepth(km),area(km),p_z(km),vol(km) )
open(1,file=geometry2_file,access='direct',form='unformatted',recl=6,     &
       status='old')
do k=1,km
  read(1,rec=k) geometry2(:,k) ! all real: k, dz[m], tdepth[m], area[m^2], p[bar]
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
p_z(:)    = geometry2(5,:)
vol(:)    = geometry2(6,:)
close(1)

!===============================================================================
!  read fields
!===============================================================================

! 3D
allocate( RHO(imt,jmt,km),  Q(imt,jmt,km), PD(imt,jmt,km),                     &
          SALT(imt,jmt,km), TEMP(imt,jmt,km), D0(imt,jmt,km), Dp(imt,jmt,km) ) 

! open file
open(1,file=tavg_file,access='direct',form='unformatted',recl=imt*jmt,         &
       status='old')

! read 3-D fields
call load_3D_field(imt,jmt,km,1,nrec_RHO, RHO ) ! [g/cm^3]
call load_3D_field(imt,jmt,km,1,nrec_PD,  PD  ) ! [g/cm^3]
call load_3D_field(imt,jmt,km,1,nrec_Q,   Q   ) ! [g/cm^4]
call load_3D_field(imt,jmt,km,1,nrec_SALT,SALT) ! [g/kg]
call load_3D_field(imt,jmt,km,1,nrec_TEMP,TEMP) ! [degC]

close(1)

! testing state
!  (35, 25, 2000) -> 1031.654 229
!  ( 0, 20, 1000) -> 1017.726 743 something wrong here
!  (40, 12, 8000) -> 1062.928 258
call state(dble(35)/1.0e3,dble(25),dble(200),D_test)
write(*,*) D_test+rho0
call state(dble(0)/1.0e3,dble(20),dble(100),D_test)
write(*,*) D_test+rho0
call state(dble(40)/1.0e3,dble(12),dble(800),D_test)
write(*,*) D_test+rho0
!  test value : rho = 1.033213242 for
!  S = 35.0 PSU, theta = 20.0, pressz = 200.0
call state(dble(35)/1.0e3,dble(20),dble(200),D_test)
write(*,*) D_test+rho0


!stop

! create new potential density fields
do k=1,km
  do j=1,jmt
    do i=1,imt
      if ( DZT(i,j,k).ne.0.0 ) then
        call state(dble(SALT(i,j,k)),dble(TEMP(i,j,k)),dble(0)     ,D_test)
        D0(i,j,k) = real(D_test)
        call state(dble(SALT(i,j,k)),dble(TEMP(i,j,k)),dble(p_z(k)),D_test)
        Dp(i,j,k) = real(D_test)
      endif
    enddo
  enddo
enddo

!===============================================================================
!  calculate averaged quantities
!===============================================================================

allocate( T_avg(km),  S_avg(km), RHO_avg(km), PD_avg(km), Q_avg(km), D0_avg(km),&
          Dp_avg(km),  PD_ddz(km), RHO_ddz(km), D0_ddz(km), Dp_ddz(km),         &
          T_wavg(km), S_wavg(km), RHO_wavg(km), PD_wavg(km), Q_wavg(km), D0_wavg(km),&
          Dp_wavg(km), PD_wddz(km), RHO_wddz(km), D0_wddz(km), Dp_wddz(km) )

! k, dz, tdepth, area, p,
! <T>, <S>, <RHO>, <PD>, <D(D,T,0)>, <D(D,T,p)>, <Q>,
! dz(<PD>), dz(<D(S,T,0)>), dz(<D(S,T,p)>)
call area_avg(TEMP       ,DZT,TAREA,area,  T_avg)
call area_avg(SALT       ,DZT,TAREA,area,  S_avg)
call area_avg( RHO*1.0E03,DZT,TAREA,area,RHO_avg)
call area_avg(  PD*1.0E03,DZT,TAREA,area, PD_avg)
call area_avg(   Q*1.0E05,DZT,TAREA,area,  Q_avg)
call area_avg(  D0       ,DZT,TAREA,area, D0_avg)
call area_avg(  Dp       ,DZT,TAREA,area, Dp_avg)

call area_avg_weighted(TEMP       ,DZT,TAREA,vol,  T_wavg)
call area_avg_weighted(SALT       ,DZT,TAREA,vol,  S_wavg)
call area_avg_weighted( RHO*1.0E03,DZT,TAREA,vol,RHO_wavg)
call area_avg_weighted(  PD*1.0E03,DZT,TAREA,vol, PD_wavg)
call area_avg_weighted(   Q*1.0E05,DZT,TAREA,vol,  Q_wavg)
call area_avg_weighted(  D0       ,DZT,TAREA,vol, D0_wavg)
call area_avg_weighted(  Dp       ,DZT,TAREA,vol, Dp_wavg)

do k=1,km
  write(*,*) k,T_avg(k)/T_wavg(k),S_avg(k)/S_wavg(k),RHO_avg(k)/RHO_wavg(k),PD_avg(k)/PD_wavg(k),Q_avg(k)/Q_wavg(k),D0_avg(k)/D0_wavg(k),Dp_avg(k)/Dp_wavg(k)
enddo

!stop

call vert_der( PD_avg,dz, PD_ddz)
call vert_der(RHO_avg,dz,RHO_ddz)
call vert_der( D0_avg,dz, D0_ddz)
call vert_der( Dp_avg,dz, Dp_ddz)

call vert_der( PD_wavg,dz, PD_wddz)
call vert_der(RHO_wavg,dz,RHO_wddz)
call vert_der( D0_wavg,dz, D0_wddz)
call vert_der( Dp_wavg,dz, Dp_wddz)
!write (*,*) '-- n0 at surface and bottom corrected by factor 2 --'

!===============================================================================
!  OUTPUT
!===============================================================================

! create readme
! >>>

! k, dz, tdepth, area, p,
! <T>, <S>, <RHO>, <PD>, <D(D,T,0)>, <D(D,T,p)>,
! <Q>, dz(<RHO>), dz(<PD>), dz(<D(S,T,0)>), dz(<D(S,T,p)>)
! <Q_w>, dz(<RHO>_w), dz(<PD>_w), dz(<D(S,T,0)>_w), dz(<D(S,T,p)>_w)

open(1,file=ref_state_file,access='direct',form='unformatted',recl=21,         &
       status='unknown')
do k=1,km
  write(1,rec=k) real(k),dz(k),tdepth(k),area(k),p_z(k),                       &
                 T_avg(k),S_avg(k),RHO_avg(k),PD_avg(k),D0_avg(k),Dp_avg(k),   &
                 Q_avg(k), RHO_ddz(k), PD_ddz(k), D0_ddz(k), Dp_ddz(k),        &
                 Q_wavg(k),RHO_wddz(k),PD_wddz(k),D0_wddz(k),Dp_wddz(k)
enddo
close(1)
!===============================================================================
!===============================================================================
contains

!===============================================================================
subroutine vert_der(avg,dz,ddz)
! vertical derivative

integer                           :: k
real,   dimension(:), intent(in)  :: avg, dz
real,   dimension(:), intent(out) :: ddz

do k = 1, km
  if ( k==1 ) then
    ddz(k) = ( avg(k  ) - avg(k+1) ) / dz(k)
  else if ( k==km ) then
    ddz(k) = ( avg(k-1) - avg(k  ) ) / dz(k)
  else
    ddz(k) = ( avg(k-1) - avg(k+1) ) / 2.0 / dz(k)
  endif
enddo !k

end subroutine vert_der
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
!real, intent(in) :: & 
  T,               &! pot. temperature at level k  [degC]
  S_orig,          &! salinity at level k         [kg/kg]
  P_orig            ! pressure at reference level   [bar]


! !OUTPUT PARAMETERS:

double precision, intent(out) :: &
!real, intent(out) :: &
  RHO               ! density of water [kg m^-3]

!-----------------------------------------------------------------------
!  local variables:
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
!  first check for valid range if requested
!-----------------------------------------------------------------------

! THIS IS SKIPPED HERE

!-----------------------------------------------------------------------
!  MWJF EOS coefficients
!-----------------------------------------------------------------------

!*** these constants will be used to construct the numerator
   
double precision, parameter ::   &
!real, parameter ::                                                     &
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
!real, parameter ::                                                     &
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
!real ::                                                     &
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
!-----------------------------------------------------------------------

!   unit conversion
P   = c10*P_orig   ! [bar] -> [dbar]
S   = c1000*S_orig ! [kg/kg] -> [PSU]
SQR = sqrt(S)      ! square root 

!      write (*,*) 'S =', S, 'PSU'
!      write (*,*) 'T =', T, 'degC'
!      write (*,*) 'P =', P, 'dbar'

!*** first calculate numerator of MWJF density [P_1(S,T,p)]
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

!*** now calculate denominator of MWJF density [P_2(S,T,p)]
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

RHO   = WORK1/WORK2-rho0

! SKIPPED RHOFULL, DRHODT, DRHODS CALCULATIONS HERE
! SKIPPED OTHER EOS FUNCTIONS

end subroutine state
!===============================================================================

!===============================================================================
!===============================================================================
end program ref_state
