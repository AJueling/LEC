!===============================================================================
function pressure(depth) result(press)

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

real :: press   ! pressure [bars]

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

press = pc1*(exp(-pc2*depth) - c1) + pc3*depth + pc4*depth**2

!return

end function pressure
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

double precision, parameter :: rho0=4.1/3.996*1000

! !INPUT PARAMETERS:

double precision, intent(in) :: &
  T,               &! pot. temperature at level k  [degC]
  S_orig,          &! salinity at level k         [kg/kg]
  P_orig            ! pressure at reference level   [bar]


! !OUTPUT PARAMETERS:

double precision, intent(out) :: &
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

double precision ::&
  mwjfnums0t0, mwjfnums0t1, mwjfnums0t2, mwjfnums0t3,                          &
  mwjfnums1t0, mwjfnums1t1, mwjfnums2t0,                                       &
  mwjfdens0t0, mwjfdens0t1, mwjfdens0t2, mwjfdens0t3, mwjfdens0t4,             &
  mwjfdens1t0, mwjfdens1t1, mwjfdens1t3,                                       &
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

WORK1 = mwjfnums0t0 + T * (mwjfnums0t1 + T * (mwjfnums0t2 +                    &
        mwjfnums0t3 * T)) + S * (mwjfnums1t0 +                                 &
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

WORK2 = mwjfdens0t0 + T * (mwjfdens0t1 + T * (mwjfdens0t2 +                    &
        T   * (mwjfdens0t3 + mwjfdens0t4 * T))) +                              &
        S   * (mwjfdens1t0 + T * (mwjfdens1t1 + T*T*mwjfdens1t3)+              &
        SQR * (mwjfdensqt0 + T*T*mwjfdensqt2))

RHO   = WORK1/WORK2-rho0

! SKIPPED RHOFULL, DRHODT, DRHODS CALCULATIONS HERE
! SKIPPED OTHER EOS FUNCTIONS

end subroutine state
!===============================================================================

subroutine hydrostatic_pressure(g,dz,RHO,SSH,HSP)
implicit none
!
!  calculate hydrostatic pressure with in-situ density
!
integer                                     :: k
integer, parameter                          :: imt=3600, jmt=2400, km=42
double precision,               intent(in)  :: g
real,    dimension(km),         intent(in)  :: dz
real,    dimension(imt,jmt,km), intent(in)  :: RHO
real,    dimension(imt,jmt),    intent(in)  :: SSH
real,    dimension(imt,jmt,km), intent(out) :: HSP
real,    parameter                          :: p5 = 0.5

do k=1,km
  if (k==1) then
    HSP(:,:,k) = g*(RHO(:,:,k)*1.0E3)*(dz(k)/2.0 + SSH(:,:)*1.0E-02)
  else
    HSP(:,:,k) = HSP(:,:,k-1) + p5*g*((RHO(:,:,k-1)+RHO(:,:,k))*1.0E03)&
                                 *(dz(k-1)+dz(k))/2.0
  endif
enddo !k
end subroutine hydrostatic_pressure

