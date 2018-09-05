!  Compile on Huygens with
!    xlf90 -O3 -o SALT_BUDGET SALT_BUDGET.f90 -lnetcdf -lnetcdff
!
   program ENERGY_BUDGET
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Program based on SALT_BUDGET.f90 version from January 2014
!  by ( Michael  Kliphuis & Matthijs den Toom (IMAU, August 2011), with comments by Dewi (IMAU, January 2014):
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   use netcdf 
   implicit none

!=====================================================================================================================================================
!  variables
!=====================================================================================================================================================

   ! user input 
   integer :: ncase, imt,jmt,km

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

   character*120 :: grid_file,kmt_file,in_depths,pbc_file,tavg_file
   character*3   :: year
   ! program output
   real ::                                                                     &
   rho_avg, pd_avg, mKE_int, eKE_int, mPE_int, ePE_int,                        &
   GPmh_int, GPms_int, GPeh_int, GPes_int, GPm_int, GPe_int, GKm_int, GKe_int, &
   cPem_int, cKem_int, cPKm_int, cPKe_int, DPm, DPe, DKm, DKe, DPm_pd, DPe_pd, &
   DKm_pd, DKe_pd, mPE_pd_int, ePE_pd_int, cPem_pd_int, cPKm_pd_int,           &
   cPKe_pd_int, test_var, t1, t2, t3,                                          &
   rPt_int, rKt_int, gPt_int, gKt_int, cPKt_int, dPt, dKt

   double precision :: volume

   real, dimension(:), allocatable ::                                          &
   mKE_sint, eKE_sint, mPE_sint, ePE_sint, cPem_sint, cKem_sint, cPKm_sint,    &
   cPKe_sint

   ! internal model variables 
   integer :: rec_length, i, ip, j, k, ncid, maskid

   integer, dimension(:,:), allocatable :: kmT

   integer, dimension(:,:,:), allocatable :: TUU_mask

   real, dimension(:,:), allocatable ::                                        &
   ALPHA01, SHF, BETA01, SFWF, STFRHO, SSFRHO, TAUX, TAUY, UTAUX, VTAUY,       &
   UVEL_k, VVEL_k, KE_k, RHO_k, RHO2_k, Q_k, PDU_k, PDV_k, WVEL_k, UVEL2_k,    &
   VVEL2_k, UV_k, UW_k,  VW_k, RHOW_k, PD_k, PD2_k, RHOU_k, RHOV_k,            &
   GPmh, GPms, GPm, GPeh, GPes, GPe, GKm, GKe, TT_GKe, TT_GKm, mKE_k, eKE_k,   &
   EVAP, PRECIP, SWREST, SSREST, RUNOFF, EVAPRHO, PRECIPRHO, SWRESTRHO,        &
   SSRESTRHO, RUNOFFRHO, SWNET, LWNET, LATENT, SENSIBLE, TSREST, TEMP_k, PDW_k,&
   gPt, gKt, TT_gKt, rKt_k

   real, dimension(:,:,:), allocatable ::                                      &
   UVEL, VVEL, KE, RHO, RHO2, Q, RHOU, RHOV, WVEL, UVEL2, VVEL2, UV, UW, VW,   &
   RHOW, PD, PD2, PDU, PDV, mKE, eKE, mPE, ePE, cPem, cKem, cPKm, cPKe,        &
   TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV, TTT_UW, TTT_VW,           &
   DRHODX, DRHODY, DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DPDDX, DPDDY,           &
   TTT_WVEL, TTT_RHOW, TTT_PDW,  TTT_mKE, TTT_eKE, mPE_pd, ePE_pd,             &
   cPem_pd, cPKm_pd, cPKe_pd, TEMP, PDW,                                       &
   rPt, rKt, cPKt, TTT_rKt

   double precision, dimension(:), allocatable ::                              &
   dz,z1, z2, area, rhoref, n0, n0_inv, pdref, tdepth

   double precision, dimension(:,:), allocatable ::                            &
   HTN, HTE, WORK, DXT, DYT, TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK2,    &
   WORK3  

   double precision, dimension(:,:,:), allocatable :: DZT

   ! Testing
   real :: WVEL_int, RHO_int

   ! WTT-calculations
   real, dimension(:,:,:), allocatable ::                                      &
   WTT_RHO, WTT_UVEL, WTT_VVEL, WTT_UVEL2, WTT_VVEL2, WTT_UV, WTT_DUDX,        &
   WTT_DUDY, WTT_DUDZ, WTT_DVDX, WTT_DVDY, WTT_DVDZ, cKem_WTT, cPKm_WTT,       &
   cPKe_WTT

   real :: cKem_WTT_int, CPKm_WTT_int, cPKe_WTT_int

   ! Parameters
   double precision, parameter ::                                              &
   c0 = 0., p5 = 0.5, c1 = 1., S0ref = 0.035, rho0 = 4.1/3.996*1000, c = 3996, &
   g = 9.806, salinity_factor = -34.7*1.0E-04, RE = 6.37E06
   ! salinity factor from http://www.cesm.ucar.edu/models/cesm1.0/cesm/
   ! cesmBbrowser/html_code/pop/constants.F90.html "flux (kg/m^2/s) to salt flux
   ! (msu*cm/s)"


   !netCDF
   include 'netcdf.inc'

   write (*,*) ''
   write (*,*) '--- ENERGY BUDGET ---'
   write (*,*) ''

!=====================================================================================================================================================
!  user input
!=====================================================================================================================================================

   write (*,*) ' enter imt, jmt, km'
   read  (*,*) imt, jmt, km
   write (*,*) ' enter name of GRID file'
   read  (*,'(a120)') grid_file
   write (*,*) ' enter name of kmT file'
   read  (*,'(a120)') kmt_file
   write (*,*) ' enter name of in_depths file'
   read  (*,'(a120)') in_depths
   write (*,*) ' enter name of DZBC file'
   read  (*,'(a120)') pbc_file
   write (*,*) ' enter name of tavg file'
   read  (*,'(a120)') tavg_file
   write (*,*) ' enter year'
   read  (*,'(a3)') year
   write (*,*) ' enter first records of fields'

   read  (*,*) nrec_UVEL, nrec_VVEL, nrec_KE, nrec_TEMP, nrec_TAUX,    &
   nrec_TAUY, nrec_PD

   write(*,*) 'Year:', year

!=====================================================================================================================================================
!  read and create horizontal grid spacing, define TAREA
!=====================================================================================================================================================

   allocate(                                                                   &
   HTN(imt,jmt), HTE(imt,jmt), HUW(imt,jmt), HUS(imt,jmt), WORK(imt,jmt),      &
   WORK2(imt,jmt), WORK3(imt,jmt), DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), &
   DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt) )

   inquire (iolength = rec_length) HTN
   open (1,file=grid_file,access='direct',form='unformatted',recl=rec_length,status='unknown')
   read (1,rec=3) HTN ! [cm]
   read (1,rec=4) HTE ! [cm]
   close(1)
   write(*,*)' read file: ',grid_file

   ! Why do we do this?
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

!=====================================================================================================================================================
!  read bathymetry and depth level variables, define DZT
!=====================================================================================================================================================

   allocate( dz(km), tdepth(km), kmT(imt,jmt), DZBC(imt,jmt) )
   allocate( DZT(imt,jmt,km), TUU_mask(imt,jmt,km) )

   ! kmT
   inquire (iolength=rec_length) kmT
   open    (1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,status='unknown')
   read    (1,rec=1) kmT
   close   (1)
   write   (*,*) 'kmT file: ', kmt_file

   ! dz
   open    (1,file=in_depths,status='old')
   do k = 1, km
   read  (1,*) dz(k) ! [cm]
   enddo
   close   (1)
   write   (*,*) 'dz(k) file: ', in_depths

   ! partial bottom cell depths
   inquire (iolength=rec_length) DZBC
   open    (1,file=pbc_file,access='direct',form='unformatted',recl=rec_length,status='unknown')
   read    (1,rec=1) DZBC ! [cm]
   close   (1)
   write   (*,*) 'PBC file: ', pbc_file

   ! create DZT-field, also used as 3D TTT-mask
   do k=1,km
     where(kmT>k)
       DZT(:,:,k) = dz(k)
     elsewhere(kmT==k)
       DZT(:,:,k) = DZBC
     elsewhere
       DZT(:,:,k) = c0
     endwhere
   enddo
   deallocate( kmT, DZBC )
   
   ! create 3D TUU-mask (1=ocean, 0=land)
   do k = 1,km
   do j = 1,jmt
   do i = 1,imt
   if ( DZT(i,j,k).ne.c0 ) then
     TUU_mask(i,j,k) = 1
   elseif ( i.lt.imt .and. DZT(i+1,j,k).ne.c0 ) then
     TUU_mask(i,j,k) = 1
   elseif ( j.lt.jmt .and. DZT(i,j+1,k).ne.c0 ) then
     TUU_mask(i,j,k) = 1
   elseif ( i.lt.imt .and. j.lt.jmt .and. DZT(i+1,j+1,k).ne.c0 ) then
     TUU_mask(i,j,k) = 1
   else
     TUU_mask(i,j,k) = 0
   endif
   enddo
   enddo
   enddo

!  create tdepth 1D array
   tdepth(1) = dz(1)/2.0
   do k = 2,km
     tdepth(k) = tdepth(k-1)+p5*(dz(k-1)+dz(k))
   enddo

!=====================================================================================================================================================
!  read fields
!=====================================================================================================================================================


   ! 3D
   allocate(                                                                   &
   UVEL_k(imt,jmt),  VVEL_k(imt,jmt),                          &
   UVEL(imt,jmt,km), VVEL(imt,jmt,km),                        &
   PD_k(imt,jmt),    PD(imt,jmt,km) )

   !open file
   inquire (iolength=rec_length) UVEL_k
   open (1,file=tavg_file,access='direct',form='unformatted',recl=rec_length,status='unknown')

   !then read 3-D fields, storing full salinity field, and velocity field for sections
   do k = 1, km
!     read (1,rec=nrec_RHO+k-1)   RHO_k      ! [g/cm^3]
     read (1,rec=nrec_PD+k-1)    PD_k       ! [g/cm^3]
!     read (1,rec=nrec_Q+k-1)     Q_k        ! [g/cm^4]
     read (1,rec=nrec_UVEL+k-1)  UVEL_k     ! [cm/s]
     read (1,rec=nrec_VVEL+k-1)  VVEL_k     ! [cm/s]
!     RHO(:,:,k)                = RHO_k
     PD(:,:,k)                 = PD_k
!     Q(:,:,k)                  = Q_k
     UVEL(:,:,k)               = UVEL_k
     VVEL(:,:,k)               = VVEL_k
   enddo

   deallocate( UVEL_k, VVEL_k, PD_k )

   write (*,*) 'file: ', tavg_file

!=====================================================================================================================================================
!  calculate averaged quantities
!=====================================================================================================================================================

   allocate( area(km), pdref(km), n0(km), n0_inv(km) )

   volume  = sum(TAREA*sum(DZT,3))*1.0E-06                                ! volume between bottom and z=0 [m^3]
!   call vol_int(imt,jmt,km,RHO,TAREA,DZT,rho_avg)
!   rho_avg = rho_avg/volume
   call vol_int(imt,jmt,km, PD,TAREA,DZT, pd_avg)
   pd_avg  =  pd_avg/volume
!   write (*,"(A10, ES13.4E2, A7)") 'rho_avg =', rho_avg, 'g/cm^3'
   write (*,"(A10, ES13.4E2, A7)") 'pd_avg  =',  pd_avg, 'g/cm^3'
   write (*,"(A10, ES13.4E2, A7)") 'rho0    =',    rho0, 'kg/m^3'
   write (*,"(A10, ES13.4E2, A4)") 'volume  =',  volume, 'm^3'
   write (*,*) ''
   write (*,*) 'k           area          pdref             n0'

   do k = 1, km
     area(k)     = sum(TAREA,            DZT(:,:,k).ne.c0)*1.0E-04         ! (T)area per level [m^2]
!     rhoref(k)   = sum(TAREA*RHO(:,:,k), DZT(:,:,k).ne.0)/area(k)*1.0E-01 ! rho_ref: average in situ density at level k [kg/m^3]
     pdref(k)    = sum(TAREA*PD(:,:,k),  DZT(:,:,k).ne.c0)/area(k)*1.0E-01 ! rho_ref: average in situ density at level k [kg/m^3]

!     n0(k)       = sum(TAREA*Q(:,:,k),   DZT(:,:,k).ne.0)/area(k)*1.0E01  ! n0: average Q at level k [kg/m^4]
!     n0_inv(k)   = c1/n0(k)                                               ! inverse n0 [m^4/kg]
   enddo

   open(3,file="n0",form='formatted',status='old',action='read')

   read (3,fmt="(41(E12.5,1x),E12.5)") n0(1), n0(2), n0(3), n0(4), n0(5), n0(6), n0(7), n0(8), n0(9), n0(10), n0(11), n0(12), n0(13), n0(14), n0(15), n0(16), n0(17), n0(18), n0(19), n0(20), n0(21), n0(22), n0(23), n0(24), n0(25), n0(26), n0(27), n0(28), n0(29), n0(30), n0(31), n0(32), n0(33), n0(34), n0(35), n0(36), n0(37), n0(38), n0(39), n0(40), n0(41), n0(42)

   do k = 1, km
     n0_inv(k)   = c1/n0(k)                                               ! inverse n0 [m^4/kg]
     write (*,"( I2.1, 3ES15.6E2)") k, area(k), pdref(k), n0(k)
   enddo

   write (*,*) ''

!=====================================================================================================================================================
!  ENERGY CALCULATIONS
!  ( from now everything will be converted into SI units )
!=====================================================================================================================================================

!=====================================================================================================================================================
!  1. RESERVOIRS
!=====================================================================================================================================================

   write (*,*) ''
   write (*,*) '1. RESERVOIRS'

!  all PE terms are on TTT-grid, all KE terms on TUU grid
   allocate( KE_k(imt,jmt), KE(imt,jmt,km) )


   do k = 1, km
!     read (1,rec=nrec_RHO2+k-1)  RHO2_k     ! [g^2/cm^6]
!     read (1,rec=nrec_PD2+k-1)   PD2_k      ! [g^2/cm^6]
     read (1,rec=nrec_KE+k-1)    KE_k       ! [cm^2/s^2]
     KE(:,:,k)                 = KE_k
!     RHO2(:,:,k)               = RHO2_k
!     PD2(:,:,k)                = PD2_k
   enddo

   allocate( mKE(imt,jmt,km), eKE(imt,jmt,km), mPE(imt,jmt,km),                &
             mPE_pd (imt,jmt,km), mKE_sint(km), eKE_sint(km), mPE_sint(km) )

   do k = 1, km
!     rPt(:,:,k)      = -p5*g*n0_inv(k) * PD2(:,:,k) * 1.0E06
!     rKt(:,:,k)      =  KE(:,:,k) 
!     mPE(:,:,k)      = -p5*g*n0_inv(k) * ( RHO(:,:,k)*1.0E03 - rhoref(k) )**2
     mPE_pd(:,:,k)   = -p5*g*n0_inv(k) * ( PD(:,:,k)*1.0E03 - pdref(k) )**2 
!     ePE(:,:,k)      = -p5*g*n0_inv(k) * (  RHO2(:,:,k)      - RHO(:,:,k)**2 ) * 1.0E06
!     ePE_pd(:,:,k)   = -p5*g*n0_inv(k) * (   PD2(:,:,k)      -  PD(:,:,k)**2 ) * 1.0E06
     mKE(:,:,k)      =  p5*rho0 * ( UVEL(:,:,k)**2 + VVEL(:,:,k)**2 ) * 1.0E-04 
     eKE(:,:,k)      =  rho0 * KE(:,:,k) * 1.0E-04 - mKE(:,:,k)
   enddo

   ! interpolate KE terms onto TTT-grid
   allocate( TTT_mKE(imt,jmt,km), TTT_eKE(imt,jmt,km),    &
             mKE_k(imt,jmt),      eKE_k(imt,jmt) )

   do k = 1, km
     call uu2tt(imt,jmt,k,mKE(:,:,k),DZT(:,:,k),mKE_k)
     call uu2tt(imt,jmt,k,eKE(:,:,k),DZT(:,:,k),eKE_k)
!     call uu2tt(imt,jmt,k,rKt(:,:,k),DZT(:,:,k),rKt_k)
     TTT_mKE(:,:,k) = mKE_k
     TTT_eKE(:,:,k) = eKE_k
!     TTT_rKt(:,:,k) = rKt_k
   enddo

   ! integrate
!   call vol_int(imt,jmt,km,    rPt,TAREA,DZT,   rPt_int)
!   call vol_int(imt,jmt,km,TTT_rKt,TAREA,DZT,   rKt_int)
   call vol_int(imt,jmt,km, mPE_pd,TAREA,DZT,mPE_pd_int)
!   call vol_int(imt,jmt,km, ePE_pd,TAREA,DZT,ePE_pd_int)
   call vol_int(imt,jmt,km,TTT_mKE,TAREA,DZT,   mKE_int)
   call vol_int(imt,jmt,km,TTT_eKE,TAREA,DZT,   eKE_int)

   ! surface integrals
   do k = 1,km
   call surf_int(imt,jmt, mPE_pd(:,:,k),TAREA,DZT(:,:,k),mPE_sint(k))
!   call surf_int(imt,jmt, ePE_pd(:,:,k),TAREA,DZT(:,:,k),ePE_sint(k))
   call surf_int(imt,jmt,TTT_mKE(:,:,k),TAREA,DZT(:,:,k),mKE_sint(k))
   call surf_int(imt,jmt,TTT_eKE(:,:,k),TAREA,DZT(:,:,k),eKE_sint(k))
   enddo

   deallocate( mKE, mKE_k, eKE, eKE_k, mPE_pd, KE_k, KE, TTT_mKE, TTT_eKE)
   ! output
300 FORMAT(A10, ES14.6E2, A2, ES10.2E2)
301 FORMAT(A10, ES14.6E2, A2, ES10.2E2, A15, ES14.6E2, A2, ES10.2E2)
   write (*,300) 'mPE_pd  =', mPE_pd_int, 'J', mPE_pd_int/1.66E23
!   write (*,301) 'ePE     =',    ePE_int, 'J',    ePE_int/6.38E18,&
!                 'ePE_pd  =', ePE_pd_int, 'J', ePE_pd_int/6.38E18
   write (*,300) 'mKE     =',    mKE_int, 'J',    mKE_int/1.27E18
   write (*,300) 'eKE     =',    eKE_int, 'J',    eKE_int/3.55E18

!=====================================================================================================================================================
!  2. GENERATION
!=====================================================================================================================================================

   write (*,*) ''
   write (*,*) '2. GENERATION'

   ! 2D
   allocate( TAUX(imt,jmt), TAUY(imt,jmt) )
!,                                             &
!   EVAP(imt,jmt), PRECIP(imt,jmt), SWREST(imt,jmt), SSREST(imt,jmt),           &
!   RUNOFF(imt,jmt), EVAPRHO(imt,jmt), PRECIPRHO(imt,jmt), SWRESTRHO(imt,jmt),  &
!   SSRESTRHO(imt,jmt), RUNOFFRHO(imt,jmt), SWNET(imt,jmt), LWNET(imt,jmt),     &
!   LATENT(imt,jmt), SENSIBLE(imt,jmt), TSREST(imt,jmt) )
   !read 2-D fields first
!   read (1,rec=nrec_ALPHA01)     ALPHA01   ! [g/cm^3/K]
!   read (1,rec=nrec_BETA01)      BETA01    ! [g/cm^3/msu]
!   read (1,rec=nrec_SHF)         SHF       ! [W/m^2]
!   read (1,rec=nrec_SFWF)        SFWF      ! [kg/m^2/s]
!   read (1,rec=nrec_STFRHO)      STFRHO    ! [C*g/cm^2/s]
!   read (1,rec=nrec_SSFRHO)      SSFRHO    ! [msu*g/cm^2/s]
   read (1,rec=nrec_TAUX)        TAUX      ! [dyne/cm^2]
   read (1,rec=nrec_TAUY)        TAUY      ! [dyne/cm^2]
!   read (1,rec=nrec_UTAUX)       UTAUX     ! [g/s^3]
!   read (1,rec=nrec_VTAUY)       VTAUY     ! [g/^3]
!   read (1,rec=nrec_EVAP)        EVAP      ! [kg/m^2/s]
!   read (1,rec=nrec_PRECIP)      PRECIP    ! [kg/m^2/s]
!   read (1,rec=nrec_SWREST)      SWREST    ! [kg/m^2/s]
!   read (1,rec=nrec_SSREST)      SSREST    ! [kg/m^2/s]
!   read (1,rec=nrec_RUNOFF)      RUNOFF    ! [kg/m^2/s]
!   read (1,rec=nrec_EVAPRHO)     EVAPRHO   ! [kg^2/m^5/s]
!   read (1,rec=nrec_PRECIPRHO)   PRECIPRHO ! [kg^2/m^5/s]
!   read (1,rec=nrec_SWRESTRHO)   SWRESTRHO ! [kg^2/m^5/s]
!   read (1,rec=nrec_SSRESTRHO)   SSRESTRHO ! [kg^2/m^5/s]
!   read (1,rec=nrec_RUNOFFRHO)   RUNOFFRHO ! [kg^2/m^5/s]
!   read (1,rec=nrec_SWNET)       SWNET     ! [W/m^2]
!   read (1,rec=nrec_LWNET)       LWNET     ! [W/m^2]
!   read (1,rec=nrec_LATENT)      LATENT    ! [W/m^2]
!   read (1,rec=nrec_SENSIBLE)    SENSIBLE  ! [W/m^2]
!   read (1,rec=nrec_TSREST)      TSREST    ! [W/m^2]

!  all PE terms quantities are on TT-grid, all KE terms on UU-grid

   allocate( GKm(imt,jmt), TT_GKm(imt,jmt) )

!   ! 2.1/2.2 Mean/Eddy Potential Energy Generation
!   gPt  = ( -g*ALPHA01*1.0E03*n0_inv(1)*(STFRHO*1.0E01 - SHF/c/rho0*rhoref(1)) ) + &
!          ( -g* BETA01*1.0E03*n0_inv(1)*(SSFRHO*1.0E01 - SFWF*salinity_factor*1.0E-02*rhoref(1)) )
!
!   GPmh = -g*ALPHA01*1.0E03*n0_inv(1) *  SHF /c/rho0                  * ( RHO(:,:,1)*1.0E03 - rhoref(1) )  ! []
!   GPms = -g* BETA01*1.0E03*n0_inv(1) * SFWF *salinity_factor*1.0E-02 * ( RHO(:,:,1)*1.0E03 - rhoref(1) )
!   GPm  = GPmh + GPms
!
!   GPeh = -g*ALPHA01*1.0E03*n0_inv(1) * ( STFRHO*1.0E01 -  SHF /c/rho0                  * rhoref(1) ) - GPmh
!   GPes = -g* BETA01*1.0E03*n0_inv(1) * ( SSFRHO*1.0E01 - SFWF *salinity_factor*1.0E-02 * rhoref(1) ) - GPms
!   GPe  = GPeh + GPes

   ! 2.3/2.4 Mean/Eddy Kinetic Energy Generation
   GKm  = ( UVEL(:,:,1)*TAUX + VVEL(:,:,1)*TAUY )*1.0E-03
!   GKe  = ( UTAUX + VTAUY )*1.0E-03 - GKm

!   call surf_int(imt,jmt,   gPt,TAREA,DZT(:,:,1), gPt_int)
!
!   call surf_int(imt,jmt,  GPmh,TAREA,DZT(:,:,1),GPmh_int)
!   call surf_int(imt,jmt,  GPms,TAREA,DZT(:,:,1),GPms_int)
!   call surf_int(imt,jmt,   GPm,TAREA,DZT(:,:,1), GPm_int)
!
!   call surf_int(imt,jmt,  GPeh,TAREA,DZT(:,:,1),GPeh_int)
!   call surf_int(imt,jmt,  GPes,TAREA,DZT(:,:,1),GPes_int)
!   call surf_int(imt,jmt,   GPe,TAREA,DZT(:,:,1), GPe_int)

!   call uu2tt(imt,jmt,k,gKt,DZT(:,:,1),TT_gKt)
   call uu2tt(imt,jmt,k,GKm,DZT(:,:,1),TT_GKm)
!   call uu2tt(imt,jmt,k,GKe,DZT(:,:,1),TT_GKe)

!   call surf_int(imt,jmt,TT_gKt,TAREA,DZT(:,:,1), gKt_int)
   call surf_int(imt,jmt,TT_GKm,TAREA,DZT(:,:,1), GKm_int)
!   call surf_int(imt,jmt,TT_GKe,TAREA,DZT(:,:,1), GKe_int)

   deallocate( GKm, TT_GKm, TAUX, TAUY )
!   EVAP, PRECIP, SWREST, SSREST, RUNOFF, EVAPRHO, PRECIPRHO, SWRESTRHO,        &
!   SSRESTRHO, RUNOFFRHO, SWNET, LWNET, LATENT, SENSIBLE, TSREST )

!   write (*,301) '|GPmh   =', GPmh_int, 'W', GPmh_int/2.00E12,&
!                 '{|GPmh  =', GPmh_int, 'W', GPmh_int/2.00E12
!   write (*,301) '|GPms   =', GPms_int, 'W', GPms_int/2.00E12,&
!                 '{|GPms  =', GPms_int, 'W', GPms_int/2.00E12
!   write (*,301) 'GPm     =',  GPm_int, 'W',  GPm_int/2.00E12,&
!                 '{GPm    =',  GPm_int, 'W',  GPm_int/2.00E12

!   write (*,301) '|GPeh   =', GPeh_int, 'W', GPeh_int/ 5.8E11,&
!                 '{|GPeh  =', GPeh_int, 'W', GPeh_int/ 5.8E11
!   write (*,301) '|GPes   =', GPes_int, 'W', GPes_int/ 5.8E11,&
!                 '{|GPes  =', GPes_int, 'W', GPes_int/ 5.8E11
!   write (*,301) 'GPe     =',  GPe_int, 'W',  GPe_int/ 5.8E11,&
!                 '{GPe    =',  GPe_int, 'W',  GPe_int/ 5.8E11


   write (*,300) 'GKm     =',  GKm_int, 'W',  GKm_int/1.85E12
!   write (*,301) 'GKe     =',  GKe_int, 'W',  GKe_int/2.19E12,&
!                 '{GKe    =',  GKe_int, 'W',  GKe_int/2.19E12


!=====================================================================================================================================================
!  3. CONVERSION
!=====================================================================================================================================================

!   write (*,*) ''
!   write (*,*) '3. CONVERSION'
!
!   allocate(                                                                   &
!   UVEL2_k(imt,jmt),  VVEL2_k(imt,jmt),  UV_k(imt,jmt),                        &
!   VW_k(imt,jmt),     UW_k(imt,jmt),     WVEL_k(imt,jmt),                      &
!   RHOU_k(imt,jmt),   RHOV_k(imt,jmt),   RHOW_k(imt,jmt),                      &
!   PDU_k(imt,jmt),    PDV_k(imt,jmt),    PDW_k(imt,jmt),                       &
!   UVEL2(imt,jmt,km), VVEL2(imt,jmt,km), UV(imt,jmt,km),                       &
!   VW(imt,jmt,km),    UW(imt,jmt,km) ,   WVEL(imt,jmt,km),                     &
!   RHOU(imt,jmt,km),  RHOV(imt,jmt,km),  RHOW(imt,jmt,km),                     &
!   PDU(imt,jmt,km),   PDV(imt,jmt,km),   PDW(imt,jmt,km) )
!
!   do k = 1, km
!     read (1,rec=nrec_WVEL+k-1)  WVEL_k     ! [cm/s]
!     read (1,rec=nrec_UVEL2+k-1) UVEL2_k    ! [cm^2/s^2]
!     read (1,rec=nrec_VVEL2+k-1) VVEL2_k    ! [cm^2/s^2]
!     read (1,rec=nrec_UV+k-1)    UV_k       ! [cm^2/s^2]
!     read (1,rec=nrec_UW+k-1)    UW_k       ! [cm^2/s^2]
!     read (1,rec=nrec_VW+k-1)    VW_k       ! [cm^2/s^2]
!     read (1,rec=nrec_RHOU+k-1)  RHOU_k     ! [g/cm^2/s]
!     read (1,rec=nrec_RHOV+k-1)  RHOV_k     ! [g/cm^2/2]
!     read (1,rec=nrec_RHOW+k-1)  RHOW_k     ! [g/cm^2/s]
!     read (1,rec=nrec_PDU+k-1)   PDU_k      ! [g/cm^2/s]
!     read (1,rec=nrec_PDV+k-1)   PDV_k      ! [g/cm^2/s]
!     read (1,rec=nrec_PDW+k-1)   PDW_k      ! [g/cm^2/s]
!     WVEL(:,:,k)               = WVEL_k
!     UVEL2(:,:,k)              = UVEL2_k
!     VVEL2(:,:,k)              = VVEL2_k
!     UV(:,:,k)                 = UV_k
!     UW(:,:,k)                 = UW_k
!     VW(:,:,k)                 = VW_k
!     RHOU(:,:,k)               = RHOU_k
!     RHOV(:,:,k)               = RHOV_k
!     RHOW(:,:,k)               = RHOW_k
!     PDU(:,:,k)                = PDU_k
!     PDV(:,:,k)                = PDV_k
!     PDW(:,:,k)                = PDW_k
!   enddo
!   deallocate(                                                                 &
!   UVEL2_k, VVEL2_k, UV_k,   VW_k,  UW_k,  WVEL_k,                             &
!   RHOU_k,  RHOV_k,  RHOW_k, PDU_k, PDV_k, PDW_k )
!
!! first create new fields:
!! 3.1: (tuu->ttt) uvel,vvel;  (quadratic? central difference) drho/dx, drho/dy,
!!      (nabla, horizontally linear central difference) du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz
!! 3.2: (tuu->ttt) uvel2,vvel2,uvelvel
!! 3.3: (wtt->ttt) wvel
!! 4.4: (wtt->ttt) wvelrho
!
!! Interpolation
!
!! > TUU->TTT: TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV
!   allocate(                                                                   &
!   TTT_UVEL(imt,jmt,km), TTT_VVEL(imt,jmt,km), TTT_UVEL2(imt,jmt,km),          &
!   TTT_VVEL2(imt,jmt,km), TTT_UV(imt,jmt,km), UVEL_k(imt,jmt),                 &
!   VVEL_k(imt,jmt), UVEL2_k(imt,jmt), VVEL2_k(imt,jmt), UV_k(imt,jmt) )
!   do k = 1, km
!   call uu2tt(imt,jmt,k,UVEL(:,:,k),DZT(:,:,k),UVEL_k)
!   call uu2tt(imt,jmt,k,VVEL(:,:,k),DZT(:,:,k),VVEL_k)
!   call uu2tt(imt,jmt,k,UVEL2(:,:,k),DZT(:,:,k),UVEL2_k)
!   call uu2tt(imt,jmt,k,VVEL2(:,:,k),DZT(:,:,k),VVEL2_k)
!   call uu2tt(imt,jmt,k,UV(:,:,k),DZT(:,:,k),UV_k)
!   TTT_UVEL(:,:,k)  = UVEL_k
!   TTT_VVEL(:,:,k)  = VVEL_k
!   TTT_UVEL2(:,:,k) = UVEL2_k
!   TTT_VVEL2(:,:,k) = VVEL2_k
!   TTT_UV(:,:,k)    = UV_k
!   enddo
!   deallocate( UVEL_k, VVEL_k, UVEL2_k, VVEL2_k, UV_k )
!
!! > WTT->TTT: TTT_UW, TTT_VW, TTT_WVEL, TTT_RHOW
!   allocate(                                                                   &
!   TTT_UW(imt,jmt,km), TTT_VW(imt,jmt,km), TTT_WVEL(imt,jmt,km),               &
!   TTT_RHOW(imt,jmt,km), TTT_PDW(imt,jmt,km) )
!   call wtt2ttt(imt,jmt,km,  UW,DZT,  TTT_UW)
!   call wtt2ttt(imt,jmt,km,  VW,DZT,  TTT_VW)
!   call wtt2ttt(imt,jmt,km,WVEL,DZT,TTT_WVEL)
!   call wtt2ttt(imt,jmt,km,RHOW,DZT,TTT_RHOW)
!   call wtt2ttt(imt,jmt,km, PDW,DZT, TTT_PDW)
!
!! Derivatives
!
!! > DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ
!   allocate(                                                                   &
!   DUDX(imt,jmt,km), DUDY(imt,jmt,km), DUDZ(imt,jmt,km), DVDX(imt,jmt,km),     &
!   DVDY(imt,jmt,km), DVDZ(imt,jmt,km) )
!   call nabla_hvel(imt,jmt,km,UVEL,TTT_UVEL,dz,DXT,DYT,DZT,DUDX,DUDY,DUDZ) ![1/s]
!   call nabla_hvel(imt,jmt,km,VVEL,TTT_VVEL,dz,DXT,DYT,DZT,DVDX,DVDY,DVDZ)
!
!! > DRHODX, DRHODY
!   allocate(                                                                   &
!   DRHODX(imt,jmt,km), DRHODY(imt,jmt,km), DPDDX(imt,jmt,km), DPDDY(imt,jmt,km))
!   call hdiv_rho(imt,jmt,km,RHO,DXT,DYT,DZT,HUW,HUS,DRHODX,DRHODY)
!   call hdiv_rho(imt,jmt,km, PD,DXT,DYT,DZT,HUW,HUS, DPDDX, DPDDY)
!
!!=====================================================================================================================================================
!   allocate(                                                                   &
!   cPem(imt,jmt,km), cKem(imt,jmt,km), cPKm(imt,jmt,km), cPKe(imt,jmt,km),     &
!   cPem_pd(imt,jmt,km), cPKm_pd(imt,jmt,km) , cPKe_pd(imt,jmt,km),             &
!   cPem_sint(km), cKem_sint(km), cPKm_sint(km), cPKe_sint(km) ,                &
!   cPKt(imt,jmt,km) )
!
!   do k = 1, km
!   ! 3.1/3.2 Eddy to Mean Energy
!     cPem(:,:,k)    = -g*n0_inv(k)*1.0E06 *                                    &
!                   ((RHOU(:,:,k)-RHO(:,:,k)*TTT_UVEL(:,:,k))*DRHODX(:,:,k)+    &
!                    (RHOV(:,:,k)-RHO(:,:,k)*TTT_VVEL(:,:,k))*DRHODY(:,:,k))
!     cPem_pd(:,:,k) = -g*n0_inv(k)*1.0E06 *                                    &
!                   ((PDU(:,:,k)-PD(:,:,k)*TTT_UVEL(:,:,k))*DPDDX(:,:,k)+       &
!                    (PDV(:,:,k)-PD(:,:,k)*TTT_VVEL(:,:,k))*DPDDY(:,:,k))
!     cKem(:,:,k)    = rho0*1.0E-04 *                                           &
!                   ((TTT_UVEL2(:,:,k) - TTT_UVEL(:,:,k)**2)*DUDX(:,:,k) +      &
!                    (TTT_VVEL2(:,:,k) - TTT_VVEL(:,:,k)**2)*DVDY(:,:,k) +      &
!                    (TTT_UV(:,:,k)    - TTT_UVEL(:,:,k)*TTT_VVEL(:,:,k)) *     &
!                                                   (DUDY(:,:,k)+DVDX(:,:,k)) + &
!                    (TTT_UW(:,:,k)    - TTT_UVEL(:,:,k)*TTT_WVEL(:,:,k))*      &
!                                                                   DUDZ(:,:,k)+&
!                    (TTT_VW(:,:,k)    - TTT_VVEL(:,:,k)*TTT_WVEL(:,:,k)) *     &
!                                                                   DVDZ(:,:,k))
!
!   ! 3.3/3.4 Potential to Kinetic Energy
!     cPKt(:,:,k)    = -g * TTT_PDW(:,:,k)*1.0E01 
!
!     cPKm(:,:,k)    = -g * RHO(:,:,k) * TTT_WVEL(:,:,k) * 1.0E01 
!     cPKm_pd(:,:,k) = -g *  PD(:,:,k) * TTT_WVEL(:,:,k) * 1.0E01 
!     cPKe(:,:,k)    = -g * TTT_RHOW(:,:,k)*1.0E01  - cPKm(:,:,k)
!     cPKe_pd(:,:,k) = -g *  TTT_PDW(:,:,k)*1.0E01  - cPKm_pd(:,:,k)
!   enddo
!
!   deallocate(                                                                 &
!   UVEL2,    VVEL2,    UV, VW, UW, WVEL, RHOU, RHOV, RHOW, PDU, PDV, PDW,      &
!   TTT_UVEL, TTT_VVEL, TTT_UVEL2, TTT_VVEL2, TTT_UV,                           &
!   TTT_UW,   TTT_VW, TTT_WVEL, TTT_RHOW, TTT_PDW,                              &
!   DUDX,     DUDY, DUDZ, DVDX, DVDY, DVDZ, DRHODX, DRHODY, DPDDX, DPDDY )
!
!   CALL vol_int(imt,jmt,km,   CPEM,TAREA,DZT,   CPEM_INT)
!   CALL vol_int(imt,jmt,km,CPEM_PD,TAREA,DZT,CPEM_PD_INT)
!   CALL vol_int(imt,jmt,km,   CKEM,TAREA,DZT,   CKEM_INT)
!   call vol_int(imt,jmt,km,   cPKt,TAREA,DZT,   cPKt_int)
!   call vol_int(imt,jmt,km,   cPKm,TAREA,DZT,   cPKm_int)
!   CALL vol_int(imt,jmt,km,CPkm_PD,TAREA,DZT,CPkm_PD_INT)
!   CALL vol_int(imt,jmt,km,   CPKE,TAREA,DZT,   CPKE_INT)
!   CALL vol_int(imt,jmt,km,CPKE_PD,TAREA,DZT,CPKE_PD_INT)
!
!   ! surface integrals
!   do k = 1,km
!   call surf_int(imt,jmt,cPem_pd(:,:,k),TAREA,DZT(:,:,k),cPem_sint(k))
!   call surf_int(imt,jmt,   cKem(:,:,k),TAREA,DZT(:,:,k),cKem_sint(k))
!   call surf_int(imt,jmt,cPKm_pd(:,:,k),TAREA,DZT(:,:,k),cPKm_sint(k))
!   call surf_int(imt,jmt,cPKe_pd(:,:,k),TAREA,DZT(:,:,k),cPKe_sint(k))
!   enddo
!
!   deallocate( cPem, cPem_pd, cKem, cPKm, cPKm_pd, cPKe, cPKe_pd )
!
!   write (*,301) 'cPem    =',    cPem_int, 'W',    CPem_int/-8.3E11,&
!                 'cPem_pd =', cPem_pd_int, 'W', CPem_pd_int/-8.3E11
!   write (*,301) 'cKem    =',    cKem_int, 'W',    CKem_int/-1.1E11,&
!                 '{cKem   =',    cKem_int, 'W',    CKem_int/-1.1E11
!   write (*,301) 'cPKm    =',    cPKm_int, 'W',    CPKm_int/-4.9E11,&
!                 'cPKm_pd =', cPKm_pd_int, 'W', CPKm_pd_int/-4.9E11
!   write (*,301) 'cPKe    =',    cPKe_int, 'W',    CPKe_int/ 7.3E11,&
!                 'cPKe_pd =', cPKe_pd_int, 'W', CPKe_pd_int/ 7.3E11
!
!!  clsoing the output file
!   close(1)
!
!!=====================================================================================================================================================
!!  4. DISSIPATION
!!=====================================================================================================================================================
!
!   write (*,*) ''
!   write (*,*) '4. DISSIPATION'
!
!   dPt    = gPt_int - cPKt_int
!   DPm    = GPm_int + cPem_int    - cPKm_int
!   DPm_pd = GPm_int + cPem_pd_int - cPKm_pd_int
!   DPe    = GPe_int - cPem_int    - cPKe_int
!   DPe_pd = GPe_int - cPem_pd_int - cPKe_pd_int
!   dKt    = gKt_int + cPKt_int
!   DKm    = GKm_int + cKem_int    + cPKm_int
!   DKm_pd = GKm_int + cKem_int    + cPKm_pd_int
!   DKe    = GKe_int - CKem_int    + cPKe_int
!   DKe_pd = GKe_int - CKem_int    + cPKe_pd_int
!
!   write (*,301) 'DPm     =', DPm   , 'W',    DPm/1.66E12,&
!                 'DPm_pd  =', DPm_pd, 'W', DPm_pd/1.66E12
!   write (*,301) 'DPe     =', DPe   , 'W',    DPe/6.80E11,&
!                 'DPe_pd  =', DPe_pd, 'W', DPe_pd/6.80E11
!   write (*,301) 'DKm     =', DKm   , 'W',    DKm/1.36E12,&
!                 'DKm_pd  =', DKm_pd, 'W', DKm_pd/1.36E12
!   write (*,301) 'DKe     =', DKe   , 'W',    DKe/3.03E12,&
!                 'DKe_pd  =', DKe_pd, 'W', DKe_pd/3.03E12
!
!
!   write (*,*) ''
!   write (*,*) ' TOTAL ENERGY'
!   write (*,300) 'gPt     =',  gPt_int, 'W',  gPt_int/2.58E12
!   write (*,301) 'rPt     =',  rPt_int, 'J',  rPt_int/1.66E23,&
!                 'dPt     =',  dPt    , 'W',      dPm/2.34E12
!   write (*,300) 'cPKt    =', cPKt_int, 'W', cPKt_int/ 2.4E11
!   write (*,301) 'rKt     =',  rKt_int, 'J',  rKt_int/4.83E18,&
!                 'dKt     =',  dKt    , 'W',      dKt/3.39E12
!   write (*,300) 'gKt     =',  gKt_int, 'W',  gKt_int/4.04E12
!!=====================================================================================================================================================
!!  CALCULATION ON WTT-GRID
!!=====================================================================================================================================================
!
!! > TTT->WTT: WTT_RHO, WTT_UVEL, WTT_VVEL,
!!             WTT_UVEL2, WTT_VVEL2, WTT_UV,
!!             WTT_DUDX, WTT_DUDY, WTT_DUDZ,
!             WTT_DVDX, WTT_DVDY, WTT_DVDZ,
!   to create cKem_WTT, cPKm_WTT, cPKe_WTT,
!             cKem_WTT_int, CPKm_WTT_int, cPKe_WTT_int
!   allocate( WTT_RHO(imt,jmt,km), WTT_UVEL(imt,jmt,km), WTT_VVEL(imt,jmt,km),  &
!             WTT_UVEL2(imt,jmt,km), WTT_VVEL2(imt,jmt,km), WTT_UV(imt,jmt,km), &
!             WTT_DUDX(imt,jmt,km), WTT_DUDY(imt,jmt,km), WTT_DUDZ(imt,jmt,km), &
!             WTT_DVDX(imt,jmt,km), WTT_DVDY(imt,jmt,km), WTT_DVDZ(imt,jmt,km), &
!             cKem_WTT(imt,jmt,km), cPKm_WTT(imt,jmt,km), cPKe_WTT(imt,jmt,km) )
!
!   call ttt2wtt(imt,jmt,km,      RHO,dz,DZT,  WTT_RHO)
!   call ttt2wtt(imt,jmt,km, TTT_UVEL,dz,DZT, WTT_UVEL)
!   call ttt2wtt(imt,jmt,km, TTT_VVEL,dz,DZT, WTT_VVEL)
!   call ttt2wtt(imt,jmt,km,TTT_UVEL2,dz,DZT,WTT_UVEL2)
!   call ttt2wtt(imt,jmt,km,TTT_VVEL2,dz,DZT,WTT_VVEL2)
!   call ttt2wtt(imt,jmt,km,   TTT_UV,dz,DZT,   WTT_UV)
!   call ttt2wtt(imt,jmt,km,     DUDX,dz,DZT, WTT_DUDX)
!   call ttt2wtt(imt,jmt,km,     DUDY,dz,DZT, WTT_DUDY)
!   call ttt2wtt(imt,jmt,km,     DUDZ,dz,DZT, WTT_DUDZ)
!   call ttt2wtt(imt,jmt,km,     DVDX,dz,DZT, WTT_DVDX)
!   call ttt2wtt(imt,jmt,km,     DVDY,dz,DZT, WTT_DVDY)
!   call ttt2wtt(imt,jmt,km,     DVDZ,dz,DZT, WTT_DVDZ)
!
!   do k = 1, km
!   do j = 1, jmt
!   do i = 1, imt
!   if ( DZT(i,j,k).ne.0 ) then
!   ! 3.1/3.2 Eddy to Mean Potential/Kinetic Energy
!     cKem_WTT(i,j,k) = rho0*1.0E-04 *                                          &
!                   ((WTT_UVEL2(i,j,k)-WTT_UVEL(i,j,k)**2)*WTT_DUDX(i,j,k)+     &
!                    (WTT_VVEL2(i,j,k)-WTT_VVEL(i,j,k)**2)*WTT_DVDY(i,j,k)+     &
!                    (   WTT_UV(i,j,k)-WTT_UVEL(i,j,k)*WTT_VVEL(i,j,k))*        &
!                                            (WTT_DUDY(i,j,k)+ WTT_DVDX(i,j,k))+&
!                    ( UW(i,j,k)-WTT_UVEL(i,j,k)*WVEL(i,j,k))* WTT_DUDZ(i,j,k) +&
!                    ( VW(i,j,k)-WTT_VVEL(i,j,k)*WVEL(i,j,k))* WTT_DVDZ(i,j,k))
!   ! 3.3/3.4 Mean/Eddy Potential to Kinetic Energy
!     cPKm_WTT(i,j,k) = -g * WTT_RHO(i,j,k) * WVEL(i,j,k) * 1.0E01 
!     cPKe_WTT(i,j,k) = -g * RHOW(i,j,k)*1.0E01  - cPKm(i,j,k)
!   endif
!   enddo
!   enddo
!   enddo
!
!   call wtt2ttt(imt,jmt,km,cKem_WTT,DZT,cKem)
!   call wtt2ttt(imt,jmt,km,cPKm_WTT,DZT,cPKm)
!   call wtt2ttt(imt,jmt,km,cPKe_WTT,DZT,cPKe)
!
!   call vol_int(imt,jmt,km,cKem,TAREA,DZT,cKem_WTT_int)
!   call vol_int(imt,jmt,km,cPKm,TAREA,DZT,cPKm_WTT_int)
!   call vol_int(imt,jmt,km,cPKe,TAREA,DZT,cPKe_WTT_int)
!
!   write (*,*) '> Calculations on WTT-grid:'
!   write (*,*) 'cKem_WTT=',    cKem_WTT_int, 'W', cKem_WTT_int/-1.1E11
!   write (*,*) 'cPKm_WTT=',    cPKm_WTT_int, 'W', cPKm_WTT_int/-4.9E11
!   write (*,*) 'cPKe_WTT=',    cPKe_WTT_int, 'W', cPKe_WTT_int/ 7.3E11

!=====================================================================================================================================================
!  TESTING
!=====================================================================================================================================================

!   write (*,*) ''
!   write (*,*) 'TESTING'

!   write (*,*) 'RHO(1800,1200,1) =', RHO(1800,1200,1)
!   write (*,*) 'UVEL(1800,1200,1)  =', UVEL(1800,1200,1)
!   write (*,*) 'VVEL(1800,1200,1)  =', VVEL(1800,1200,1)
!   write (*,*) 'PD(1800,1200,1)    =', PD(1800,1200,1)
!   write (*,*) 'KE(1800,1200,1)    =', KE(1800,1200,1)
!   write (*,*) 'TAUY(1800,1200)    =', TAUY(1800,1200)
!   write (*,*) 'TAUX(1800,1200)    =', TAUX(1800,1200)

!   write (*,*) 'testing if the averaging procedure works'
!   write (*,*) 'RHO(1800,1200,1) =', RHO(1800,1200,1)
!   write (*,*) 'RHO2(1800,1200,1) =', RHO2(1800,1200,1)
!   >>> works

!   write (*,*) 'k     rhoref     n0'
!   do k =1,km
!   write (*,*) k, rhoref(k), n0(k)
!   enddo
!   >>> rhoref steadily increases downward from 1024.483 to 1054.354
!   >>> n0 is negativ with a maximumeat level 5 (n0=-1.69E-02) with larger
!   values at the top (n0(1)=-6.3E-04) and bottom (n0(km)=-1.35E-06)

!   write (*,*) 'southernmost t-cells must be zero, so that u-points lie at land boundary'
!   do i =  1,imt
!   if ( DZT(i,1,1).ne.0 ) then
!   write (*,*) 'nonzero at i =', i, '; RHO(i,1,1) = ', RHO(i,1,1)
!   endif
!   enddo
!   >>> they are indeed all of zero depth

!   write (*,*) 'k=km t-cells must be zero, so that u-points lie at land boundary'
!   WVEL_int = 0.0
!   do j =  1,jmt
!   do i =  1,imt
!   if ( DZT(i,j,km).ne.0.0 ) then
!   WVEL_int = WVEL_int + 1.0
!   endif
!   enddo
!   enddo
!   write (*,*) WVEL_int, 't-cells are nonzero at bottom k=km'
!   >>> there are 96238 t-cells at k=km that are nonzero in depth

!   write (*,*) 'surface integrals of vertical velocities'
!   do k = 1,km
!   call surf_int(imt,jmt,WVEL(:,:,k),TAREA,DZT(:,:,k),WVEL_int)
!   write (*,*) k, 'WVEL_int =', WVEL_int/area(k)/1.0E02, 'm/s'
!   enddo
!   >>> ca. 10^-13 m/s, seems correct

!   write (*,*) 'surface integrals of TTT-vertical velocities'
!   do k = 1,km
!   call surf_int(imt,jmt,TTT_WVEL(:,:,k),TAREA,DZT(:,:,k),WVEL_int)
!   write (*,*) k, 'WVEL_int =', WVEL_int/area(k)/1.0E02, 'm/s'
!   enddo
!   >>> same as above

!   write (*,*) 'volume integral of vertical velocity'
!   call vol_int(imt,jmt,km,TTT_WVEL,TAREA,DZT,WVEL_int)
!   write (*,*) 'WVEL_int =', WVEL_int/volume/1.0E02, 'm/s'
!   >>> 8.6*10^-8 m/s, seems correct
!   >>> about the same for volume integral of WVEL, instead of TTT_WVEL

!   write (*,*) 'average density'
!   call vol_int(imt,jmt,km,RHO,TAREA,DZT,WVEL_int)
!   write (*,*) 'rho_avg =', WVEL_int/volume, 'g/cm^3'
!   >>> 1.03761968245201 g/cm^3

!   write (*,*) 'volume integral of density'
!   call vol_int(imt,jmt,km,RHO,TAREA,DZT,RHO_int)
!   write (*,*) 'RHO_int =', RHO_int*1.0E03, 'kg'
!   >>> 1.38*10^21 kg, which is correct

!   write (*,*) 'testing whether DUVELDX works'
!   write (*,*) 'UVEL(1800,1200,1) =', UVEL(1800,1200,1)
!   write (*,*) 'UVEL(1801,1200,1) =', UVEL(1801,1200,1)
!   write (*,*) 'DXT(1800,1200) ='   , DXT(1800,1200)
!   write (*,*) 'DUDX(1801,1200,1) =', DUDX(1801,1200,1)
!   write (*,*) 'by hand DUDX(1801,1200,1) =', (p5*(UVEL(1801,1200,1)+UVEL(1801,1199,1))-p5*(UVEL(1800,1200,1)+UVEL(1800,1199,1)))/DXT(1801,1200)
!   >>> works

!   write (*,*) 'TESTING creation of HUS/HUW'
!   write (*,*) 'HTE(1800,1200) =', HTE(1800,1200)
!   write (*,*) 'HTN(1800,1200) =', HTN(1800,1200)
!   write (*,*) 'HUS(1800,1200) =', HUS(1800,1200)
!   write (*,*) 'HUW(1800,1200) =', HUW(1800,1200)
!   >>> works

!=====================================================================================================================================================
!   FILE OUTPUT
!=====================================================================================================================================================

   write (*,*) ''
   write (*,*) 'OUTPUT'
100 FORMAT (6(A,","),A)
101 FORMAT (6(E15.7,","),E15.7)
   open  (1,file=('output/total_energy_budget_'//year)//'.out',form='formatted')
   write (1,100) 'volume','rho0','pd_avg','1.1','1.3','1.4','2.3'
   write (1,101) volume, rho0, pd_avg, mPE_pd_int,  mKE_int, eKE_int, GKm_int
   close (1)
200 FORMAT (7(A,","),A)
201 FORMAT (7(E15.7,","),E15.7)
   open  (2,file='output/level_energy_budget_'//(year//'.out'),form='formatted')
   write (2,200) 'k','depth','area','pdref','n0','1.1','1.3','1.4'
   do k = 1,km
   !surface integrals of energy budget terms
   write (2,201)  real(k), tdepth(k), area(k),pdref(k), n0(k),mPE_sint(k), mKE_sint(k), eKE_sint(k)
   enddo
   close (2)

   write (*,*) '>>> done writing output'

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
      subroutine surf_int(imt,jmt,FIELD,TAREA,MASK,INTEGRAL)
!
!     calculates surface integral (*[m^2]) for TT-surface, using Kahan Summation
!
      implicit none

      ! input/output variables
      integer,                                 intent(in)  :: imt,jmt
      real,                dimension(imt,jmt), intent(in)  :: FIELD
      double precision,    dimension(imt,jmt), intent(in)  :: TAREA, MASK !=DZT(:,:,1)
      real,                                    intent(out) :: INTEGRAL
      real                                                 :: c, y, t

      INTEGRAL = 0.0
      c = 0.0
      do j = 1,jmt
      do i = 1,imt
      if (MASK(i,j).ne.0.0) then
        y = TAREA(i,j) * FIELD(i,j) - c
        t = INTEGRAL + y
        c = ( t -INTEGRAL ) - y
        INTEGRAL = t
      endif
      enddo
      enddo

      INTEGRAL = INTEGRAL/1.0E04

!      INTEGRAL = sum(TAREA*FIELD, MASK.ne.0.0)/1.0E04

      end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
      subroutine vol_int(imt,jmt,km,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field, using Kahan Summation
!
      implicit none

      ! input/output variables
      integer,                                    intent(in)  :: imt,jmt,km
      real,                dimension(imt,jmt,km), intent(in)  :: FIELD
      double precision,    dimension(imt,jmt   ), intent(in)  :: TAREA
      double precision,    dimension(imt,jmt,km), intent(in)  :: DZT
      real,                                       intent(out) :: INTEGRAL
      real                                                    :: c, y, t

      INTEGRAL = 0.0
      c = 0.0
      do k = 1,km
      do j = 1,jmt
      do i = 1,imt
      if (DZT(i,j,k).ne.0.0) then
        y = TAREA(i,j) *DZT(i,j,k) *  FIELD(i,j,k) - c
        t = INTEGRAL + y
        c = ( t -INTEGRAL ) - y
        INTEGRAL = t
      endif
      enddo
      enddo
      enddo

      INTEGRAL = INTEGRAL/1.0E06

!      INTEGRAL = sum(TAREA*sum(DZT*FIELD,3))/1.0E06

      end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
      subroutine hdiv_rho(imt,jmt,km,RHO,DXT,DYT,DZT,HUW,HUS,DRHODX,DRHODY)
!
!     calculates horizontal divergence of RHO
!     via average of 1st order forward and backward finite differences
!     if only one of the finite differences is available, only that is taken
!     units [g/cm^4]
!
      implicit none

      ! input/output variables
      integer,                                 intent(in)  :: imt,jmt,km
      double precision, dimension(imt,jmt   ), intent(in)  :: DXT, DYT, HUW, HUS
      double precision, dimension(imt,jmt,km), intent(in)  :: DZT
      real,             dimension(imt,jmt,km), intent(in)  :: RHO
      real,             dimension(imt,jmt,km), intent(out) :: DRHODX, DRHODY
      ! internal variables
      double precision,  parameter                         :: p5=0.5, c0=0.0


      do k = 1,km
      do j = 2,jmt
      do i = 1,imt
      if ( DZT(i,j,k).ne.c0 ) then

       ! d/dx
       if ( i==1 ) then
       ! westernmost T-points
        if     ( DZT(imt,j,k).ne.c0 .and. DZT(i+1,j,k).ne.c0 ) then
        DRHODX(i,j,k) = p5 * ( (RHO(i+1,j,k)-RHO(i,j,k))/HUS(i,j) + (RHO(i,j,k)-RHO(imt,j,k))/HUS(imt,j) )
        elseif ( DZT(imt,j,k).ne.c0                          ) then
        DRHODX(i,j,k) =                                             (RHO(i,j,k)-RHO(imt,j,k))/HUS(imt,j)
        elseif (                          DZT(i+1,j,k).ne.c0 ) then
        DRHODX(i,j,k) =        (RHO(i+1,j,k)-RHO(i,j,k))/HUS(i,j)
        else
        DRHODX(i,j,k) = 0.0
        endif

       elseif ( i==imt) then
       ! easternmost T-points
        if     ( DZT(i-1,j,k).ne.c0 .and. DZT(  1,j,k).ne.c0 ) then
        DRHODX(i,j,k) = p5 * ( (RHO(  1,j,k)-RHO(i,j,k))/HUS(i,j) + (RHO(i,j,k)-RHO(i-1,j,k))/HUS(i-1,j) )
        elseif ( DZT(i-1,j,k).ne.c0                          ) then
        DRHODX(i,j,k) =                                             (RHO(i,j,k)-RHO(i-1,j,k))/HUS(i-1,j)
        elseif (                          DZT(  1,j,k).ne.c0 ) then
        DRHODX(i,j,k) =        (RHO(  1,j,k)-RHO(i,j,k))/HUS(i,j)
        else
        DRHODX(i,j,k) = 0.0
        endif

       else
       ! other points
        if     ( DZT(i-1,j,k).ne.c0 .and. DZT(i+1,j,k).ne.c0 ) then
        DRHODX(i,j,k) = p5 * ( (RHO(i+1,j,k)-RHO(i,j,k))/HUS(i,j) + (RHO(i,j,k)-RHO(i-1,j,k))/HUS(i-1,j) )
        elseif ( DZT(i-1,j,k).ne.c0                          ) then
        DRHODX(i,j,k) =                                             (RHO(i,j,k)-RHO(i-1,j,k))/HUS(i-1,j)
        elseif (                          DZT(i+1,j,k).ne.c0 ) then
        DRHODX(i,j,k) =        (RHO(i+1,j,k)-RHO(i,j,k))/HUS(i,j)
        else
        DRHODX(i,j,k) = 0.0
        endif

       endif
 
       ! d/dy
       if ( j==jmt ) then
       ! northernmost T-points
        if ( DZT(i,j-1,k).ne.c0 ) then
        DRHODY(i,j,k) = (RHO(i,j,k)-RHO(i,j-1,k))/HUW(i,j-1)
        else
        DRHODY(i,j,k) = 0.0
        endif

       ! elseif ( j==1 ) not necessary as southernmost T-points must be on land anyways, tested and found true

       else
       ! non-northermost T-points
        if     ( DZT(i,j-1,k).ne.c0 .and. DZT(i,j+1,k).ne.c0 ) then
        DRHODY(i,j,k) = p5 * ( (RHO(i,j+1,k)-RHO(i,j,k))/HUW(i,j) + (RHO(i,j,k)-RHO(i,j-1,k))/HUW(i,j-1) )
        elseif ( DZT(i,j-1,k).ne.c0                          ) then
        DRHODY(i,j,k) =                                             (RHO(i,j,k)-RHO(i,j-1,k))/HUW(i,j-1)
        elseif (                          DZT(i,j+1,k).ne.c0 ) then
        DRHODY(i,j,k) =        (RHO(i,j+1,k)-RHO(i,j,k))/HUW(i,j)
        else
        DRHODY(i,j,k) = 0.0
        endif

       endif

      endif
      enddo
      enddo

      enddo

      end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
      subroutine nabla_hvel(imt,jmt,km,VEL,TTT_VEL,dz,DXT,DYT,DZT,DVELDX,DVELDY,DVELDZ)
!
!     calculates 3D gradients of a horizontal velocity field
!
      implicit none

      ! input/output variables
      integer,                                 intent(in)  :: imt,jmt,km
      double precision, dimension(        km), intent(in)  :: dz
      double precision, dimension(imt,jmt   ), intent(in)  :: DXT, DYT
      double precision, dimension(imt,jmt,km), intent(in)  :: DZT
      real,             dimension(imt,jmt,km), intent(in)  :: VEL, TTT_VEL
      real,             dimension(imt,jmt,km), intent(out) :: DVELDX, DVELDY, DVELDZ ! TTT gradients
      ! internal variables
      double precision                                     :: dwa_inv, dwb_inv ! inverse distance between T-points above and below
      double precision,  parameter                         :: p5=0.5, c0=0.0

      ! horizontal gradients, central difference, eq. 3.6 page 16 in 2010 POP Reference Manual
      do k = 1,km
      ! southernmost gridpoints (j=1) on TTT-grid must be 0, as there would not be UU gridpoints south of it otherwise
      do j = 2,jmt
      ! westernmost gridpoints (i=1)
        if ( DZT(1,j,k).ne.c0 ) then
        DVELDX(1,j,k) = ( p5*(VEL(1,j,k)+VEL(  1,j-1,k)) - p5*(VEL(imt,  j,k)+VEL(imt,j-1,k)) ) / DXT(i,j)
        DVELDY(1,j,k) = ( p5*(VEL(1,j,k)+VEL(imt,  j,k)) - p5*(VEL(  1,j-1,k)+VEL(imt,j-1,k)) ) / DYT(i,j)
        endif
      ! all gridpoints except western- and southernmost (i=1, j=1)
      do i = 2,imt
        if ( DZT(i,j,k).ne.c0 ) then
        DVELDX(i,j,k) = ( p5*(VEL(i,j,k)+VEL(  i,j-1,k)) - p5*(VEL(i-1,  j,k)+VEL(i-1,j-1,k)) ) / DXT(i,j)
        DVELDY(i,j,k) = ( p5*(VEL(i,j,k)+VEL(i-1,  j,k)) - p5*(VEL(  i,j-1,k)+VEL(i-1,j-1,k)) ) / DYT(i,j)
        endif
      enddo
      enddo
      enddo

      ! vertical gradient
      ! top/bottom layer
      dwb_inv = 1.0/(p5*(dz(   1)+dz( 2))) ! distance b/w first 2 T-points
      dwa_inv = 1.0/(p5*(dz(km-1)+dz(km))) ! distance b/w last 2 T-points
      do j = 1,jmt
      do i = 1,imt
        if ( DZT(i,j, 1).ne.c0 ) then
        DVELDZ(i,j, 1) = p5 * (TTT_VEL(i,j,   1)-TTT_VEL(i,j, 2))*dwb_inv
        endif
        if ( DZT(i,j,km).ne.c0 ) then
        DVELDZ(i,j,km) = p5 * (TTT_VEL(i,j,km-1)-TTT_VEL(i,j,km))*dwa_inv
        endif
      enddo 
      enddo 

      do k = 2,km-1
      dwa_inv = 1.0/(p5*(dz(k-1)+dz(  k)))
      dwb_inv = 1.0/(p5*(dz(k  )+dz(k+1)))
      do j = 1,jmt
      do i = 1,imt
        if ( DZT(i,j,k).ne.c0 ) then
        DVELDZ(i,j,k) = p5 * ( (TTT_VEL(i,j,k-1)-TTT_VEL(i,j,k))*dwa_inv + (TTT_VEL(i,j,k)-TTT_VEL(i,j,k+1))*dwb_inv )
        endif
      enddo 
      enddo 
      enddo 

      end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
      subroutine wtt2ttt(imt,jmt,km,W_FIELD,DZT,T_FIELD)
!
!     interpolated real field from WTT-grid to TTT-grid
!
      implicit none

      ! input/output variables
      integer,                                 intent(in)  :: imt,jmt,km
      double precision, dimension(imt,jmt,km), intent(in)  :: DZT
      real,             dimension(imt,jmt,km), intent(in)  :: W_FIELD ! WTT
      real,             dimension(imt,jmt,km), intent(out) :: T_FIELD ! TTT
      ! internal variables
      double precision,  parameter                         :: p5=0.5, c0=0.0, c3=3.0

      do k = 1,km-1
      do j = 1,jmt
      do i = 1,imt
      if ( DZT(i,j,k).ne.c0 ) then
      ! interpolating between two W-levels
      ! T-point is located exactly in the middle of two W-points
      ! T_FIELD(i,j,k) = W_FIELD(i,j,k) + (W_FIELD(i,j,k+1)-W_FIELD(i,j,k)) / dz(k) * dz(k)/2
        T_FIELD(i,j,k) = ( W_FIELD(i,j,k) + W_FIELD(i,j,k+1) )
      endif
      enddo
      enddo
      enddo

      ! bottom level
      do j = 1,jmt
      do i = 1,imt
      if ( DZT(i,j,km).ne.c0 ) then
      ! extrapolating to last level with gradient between last two W-levels
      ! T_FIELD(i,j,km) = W_FIELD(i,j,km) + (W_FIELD(i,j,km)-W_FIELD(i,j,km-1)) / dz(k) * dz(k)/2
        T_FIELD(i,j,km) = p5 * (c3*W_FIELD(i,j,km) - W_FIELD(i,j,km-1) )
      endif
      enddo
      enddo

      end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
      subroutine ttt2wtt(imt,jmt,km,T_FIELD,dz,DZT,W_FIELD)
!
!     interpolated real field from TTT-grid to WTT-grid, opposite of above
!
      implicit none

      ! input/output variables
      integer,                                 intent(in)  :: imt,jmt,km
      double precision, dimension(        km), intent(in)  :: dz
      double precision, dimension(imt,jmt,km), intent(in)  :: DZT
      real,             dimension(imt,jmt,km), intent(in)  :: T_FIELD ! TTT
      real,             dimension(imt,jmt,km), intent(out) :: W_FIELD ! WTT
      ! internal variables
      double precision                                     :: distance
      double precision,  parameter                         :: p5=0.5, c0=0.0, c3=3.0

      ! top level
      distance = p5*(dz(1)+dz(2))
      do j = 1,jmt
      do i = 1,imt
      if ( DZT(i,j,1).ne.c0 .and. DZT(i,j,2).ne.c0 ) then
      ! extrapolating to first level with gradient between first two W-levels
        W_FIELD(i,j,k) = T_FIELD(i,j,1) - (T_FIELD(i,j,2)-T_FIELD(i,j,1)) / distance * p5*dz(1)
      endif
      enddo
      enddo


      do k = 2,km
      distance = p5*(dz(k-1)+dz(k))
      do j = 1,jmt
      do i = 1,imt
      if ( DZT(i,j,k).ne.c0 ) then
      ! interpolating between two T-levels
      ! W-point is NOT exactly located between two T-points
      ! W_FIELD(i,j,k) = T_FIELD(i,j,k-1) + (T_FIELD(i,j,k)-T_FIELD(i,j,k-1)) / distance * dz(k-1)/2
        W_FIELD(i,j,k) = T_FIELD(i,j,k-1) + (T_FIELD(i,j,k)-T_FIELD(i,j,k-1)) / distance * p5*dz(k-1)
      endif
      enddo
      enddo
      enddo

      end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine uu2tt(imt,jmt,k,FIELD_k,DZT_k,NEW_FIELD_k)
!
!     interpolates real 2D field from UU-grid to TT-grid on T-depth levels
!
      implicit none

      ! input/output variables
      integer,                              intent(in)  :: imt,jmt,k
      double precision, dimension(imt,jmt), intent(in)  :: DZT_k
      real,             dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
      real,             dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT

      ! southernmost gridpoints (j=1) on TTT-grid must be 0, as there would not be UU gridpoints south of it otherwise
      do j = 2, jmt

      ! westernmost gridpoints
      if ( DZT_k(1,j).ne.c0 ) then
        NEW_FIELD_k(1,j) = 0.25 * ( FIELD_k(1,j) + FIELD_k(imt,j) + FIELD_k(1,j-1) + FIELD_k(imt,j-1) )
      endif

      ! most of the ocean
      do i = 2,imt
      if ( DZT_k(i,j).ne.c0 ) then
        NEW_FIELD_k(i,j) = 0.25 * ( FIELD_k(i,j) + FIELD_k(i-1,j) + FIELD_k(i,j-1) + FIELD_k(i-1,j-1) )
      endif
      enddo

      enddo

      end subroutine uu2tt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine check(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
        write(*,*) trim(nf90_strerror(status))
        stop
      endif
 
      end subroutine check
 
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   end program

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
        end do
      end do

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
        end do
      end do

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
