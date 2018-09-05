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
   integer          :: imt,jmt,km
   integer          :: nrec_TEMP, nrec_PD, nrec_SHF
   character*120    :: grid_file,kmt_file,in_depths,pbc_file,tavg_file
   character*3      :: year

   ! program output
   real             :: OHC_int, OHC_int_SS, OHC_int_CS, OHC_int_NS, pd_avg,    &
   OHC_int_700, OHC_int_SS_700, OHC_int_CS_700, OHC_int_NS_700,                &
   OHC_A,OHC_B,OHC_C,OHC_D, OHC_A_700,OHC_B_700,OHC_C_700,OHC_D_700,           &
   SST_A,SST_B,SST_C,SST_D, SHF_A,SHF_B,SHF_C,SHF_D
   double precision :: volume

   ! internal model variables 
   integer                                :: rec_length, i, ip, j, k, ncid, maskid
   integer,   dimension(:,:), allocatable :: kmT
   integer, dimension(:,:,:), allocatable :: TUU_mask
   real,        dimension(:), allocatable :: pdref
   real,      dimension(:,:), allocatable :: PD_k, TEMP_k, SHF_k
   real,    dimension(:,:,:), allocatable :: PD, TEMP, SHF, OHC

   double precision,     dimension(:), allocatable :: dz, z1, z2, area, tdepth
   double precision,   dimension(:,:), allocatable ::                          &
   HTN, HTE, WORK, DXT, DYT, TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK2,    &
   WORK3  
   double precision, dimension(:,:,:), allocatable :: DZT,DZT_SS,DZT_CS,DZT_NS,&
   DZT_700,DZT_SS_700,DZT_CS_700,DZT_NS_700,DZT_SOMA,DZT_SOMB,DZT_SOMC,DZT_SOMD


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
   write (*,*) '--- OHC ---'
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
   read  (*,*) nrec_PD, nrec_TEMP, nrec_SHF
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

!  create masks for Southern/Central/Northern Southern Hemisphere (78.5S to 55/35/0S)
!  create masks for temperature index
   allocate( DZT_SS(imt,jmt,km), DZT_CS(imt,jmt,km), DZT_NS(imt,jmt,km),       &
   DZT_700(imt,jmt,km), DZT_SS_700(imt,jmt,km), DZT_CS_700(imt,jmt,km),        &
   DZT_NS_700(imt,jmt,km), DZT_SOMA(imt,jmt,km), DZT_SOMB(imt,jmt,km),         &
   DZT_SOMC(imt,jmt,km), DZT_SOMD(imt,jmt,km) )

   do j=1,jmt
   if (j.le.520) then
     DZT_SS(:,j,:) = DZT(:,j,:)
   else
     DZT_SS(:,j,:) = c0
   endif 

   if (j.le.807) then
     DZT_CS(:,j,:) = DZT(:,j,:)
   else
     DZT_CS(:,j,:) = c0
   endif 

   if (j.le.1181) then
     DZT_NS(:,j,:) = DZT(:,j,:)
   else
     DZT_NS(:,j,:) = c0
   endif 
   enddo

   do k=1,km
   if (k.le.21) then
     DZT_700(:,:,k)    = DZT(:,:,k)
     DZT_SS_700(:,:,k) = DZT_SS(:,:,k)
     DZT_CS_700(:,:,k) = DZT_CS(:,:,k)
     DZT_NS_700(:,:,k) = DZT_NS(:,:,k)
   else
     DZT_700(:,:,k)    = c0
     DZT_SS_700(:,:,k) = c0
     DZT_CS_700(:,:,k) = c0
     DZT_NS_700(:,:,k) = c0
   endif 
   enddo

   DZT_SOMA = DZT
   DZT_SOMB = DZT
   DZT_SOMC = DZT
   DZT_SOMD = DZT

   DZT_SOMA(1:800,:,:) = c0
   DZT_SOMA(1100:imt,:,:) = c0
   DZT_SOMA(:,1:465,:) = c0
   DZT_SOMA(:,745:jmt,:) = c0

   DZT_SOMB(1:1100,:,:) = c0
   DZT_SOMB(1400:imt,:,:) = c0
   DZT_SOMB(:,1:465,:) = c0
   DZT_SOMB(:,745:jmt,:) = c0

   DZT_SOMC(1:1100,:,:) = c0
   DZT_SOMC(1400:imt,:,:) = c0
   DZT_SOMC(:,1:200,:) = c0
   DZT_SOMC(:,465:jmt,:) = c0
   
   DZT_SOMD(1:800,:,:) = c0
   DZT_SOMD(1100:imt,:,:) = c0
   DZT_SOMD(:,1:200,:) = c0
   DZT_SOMD(:,465:jmt,:) = c0

!  create tdepth 1D array
   tdepth(1) = dz(1)/2.0
   do k = 2,km
     tdepth(k) = tdepth(k-1)+p5*(dz(k-1)+dz(k))
   enddo

!=====================================================================================================================================================
!  read fields
!=====================================================================================================================================================


   ! 3D
   allocate( TEMP_k(imt,jmt), PD_k(imt,jmt), SHF_k(imt,jmt),&
             TEMP(imt,jmt,km), PD(imt,jmt,km), SHF(imt,jmt,km) )

   !open file
   inquire (iolength=rec_length) TEMP_k
   open (1,file=tavg_file,access='direct',form='unformatted',recl=rec_length,status='unknown')

   !then read 3-D fields, storing full salinity field, and velocity field for sections
   do k = 1, km
     read (1,rec=nrec_PD+k-1)     PD_k       ! [g/cm^3]
     read (1,rec=nrec_TEMP+k-1)   TEMP_k     ! [deg C]
     read (1,rec=nrec_SHF+k-1)    SHF_k      ! [deg C]
     PD(:,:,k)                  = PD_k
     TEMP(:,:,k)                = TEMP_k
     SHF(:,:,k)                 = SHF_k
   enddo

   deallocate( TEMP_k, PD_k )

   write (*,*) 'file: ', tavg_file

!=====================================================================================================================================================
!  calculate averaged quantities
!=====================================================================================================================================================

   allocate( area(km), pdref(km) )

   volume  = sum(TAREA*sum(DZT,3))*1.0E-06                                ! volume between bottom and z=0 [m^3]
   call vol_int(imt,jmt,km, PD,TAREA,DZT, pd_avg)
   pd_avg  =  pd_avg/volume
   write (*,"(A10, ES13.4E2, A7)") 'pd_avg  =',  pd_avg, 'g/cm^3'
   write (*,"(A10, ES13.4E2, A4)") 'volume  =',  volume, 'm^3'
   write (*,*) ''
   write (*,*) 'k       area     pdref'

   do k = 1, km
     area(k)     = sum(TAREA,            DZT(:,:,k).ne.0)*1.0E-04         ! (T)area per level [m^2]
     pdref(k)    = sum(TAREA*PD(:,:,k),  DZT(:,:,k).ne.0)/area(k)*1.0E-01 ! rho_ref: average in situ density at level k [kg/m^3]
     write (*,"( I2.1, 3ES11.2E2)") k, tdepth(k)/1.0E02, area(k), pdref(k)
   enddo

   write (*,*) 'setup done'

!=====================================================================================================================================================
!  CALCULATION
!=====================================================================================================================================================

   allocate( OHC(imt,jmt,km) )

   write (*,*) PD(1800,1200,1)

   if (nrec_PD.eq.301) then
   do k = 1, km
     OHC(:,:,k) = c * PD(:,:,k) * TEMP(:,:,k) * 1.0E03 ! now in [J/m^3]
   enddo
   elseif (nrec_PD.eq.995) then
   do k = 1, km
     OHC(:,:,k) = c * (rho0+PD(:,:,k)*1.0E03) * TEMP(:,:,k) ! now in [J/m^3]
   enddo
   endif

   write (*,*) 'setting up OHC done'

   ! integrate
   call vol_int(imt,jmt,km,OHC,TAREA,   DZT    ,   OHC_int    )
   call vol_int(imt,jmt,km,OHC,TAREA,   DZT_700,   OHC_int_700)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SS    ,OHC_int_SS    )
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SS_700,OHC_int_SS_700)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_CS    ,OHC_int_CS    )
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_CS_700,OHC_int_CS_700)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_NS    ,OHC_int_NS    )
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_NS_700,OHC_int_NS_700)

   ! OHC index
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMA,OHC_A) ! [J]
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMB,OHC_B)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMC,OHC_C)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMD,OHC_D)

   ! SST index
   call surf_int(imt,jmt,TEMP(:,:,1),TAREA,DZT_SOMA(:,:,1),SST_A)
   call surf_int(imt,jmt,TEMP(:,:,1),TAREA,DZT_SOMB(:,:,1),SST_B)
   call surf_int(imt,jmt,TEMP(:,:,1),TAREA,DZT_SOMC(:,:,1),SST_C)
   call surf_int(imt,jmt,TEMP(:,:,1),TAREA,DZT_SOMD(:,:,1),SST_D)

   ! SHF index
   call surf_int(imt,jmt,SHF,TAREA,DZT_SOMA(:,:,1),SHF_A)
   call surf_int(imt,jmt,SHF,TAREA,DZT_SOMB(:,:,1),SHF_B)
   call surf_int(imt,jmt,SHF,TAREA,DZT_SOMC(:,:,1),SHF_C)
   call surf_int(imt,jmt,SHF,TAREA,DZT_SOMD(:,:,1),SHF_D)

   !OHC_A = OHC_A/( sum(TAREA*sum(DZT_SOMA,3))*1.0E-06 ) ! [J]
   !OHC_B = OHC_B/( sum(TAREA*sum(DZT_SOMB,3))*1.0E-06 )
   !OHC_C = OHC_C/( sum(TAREA*sum(DZT_SOMC,3))*1.0E-06 )
   !OHC_D = OHC_D/( sum(TAREA*sum(DZT_SOMD,3))*1.0E-06 )

   SST_A = SST_A/( sum(TAREA(800:1100,465:745))*1.0E-04 ) ! [deg C]
   SST_B = SST_B/( sum(TAREA(1100:1400,465:745))*1.0E-04 )
   SST_C = SST_C/( sum(TAREA(1100:1400,200:465))*1.0E-04 )
   SST_D = SST_D/( sum(TAREA(800:1100,200:465))*1.0E-04 )

   SHF_A = SHF_A/( sum(TAREA(800:1100,465:745))*1.0E-04 ) ! [W/m^2]
   SHF_B = SHF_B/( sum(TAREA(1100:1400,465:745))*1.0E-04 )
   SHF_C = SHF_C/( sum(TAREA(1100:1400,200:465))*1.0E-04 )
   SHF_D = SHF_D/( sum(TAREA(800:1100,200:465))*1.0E-04 )

   write (*,*) 'Volume of DZT_A:', sum(TAREA*sum(DZT_SOMA,3))*1.0E-06 , 'm^3'
   write (*,*) 'Volume of DZT_B:', sum(TAREA*sum(DZT_SOMB,3))*1.0E-06 , 'm^3'
   write (*,*) 'Volume of DZT_C:', sum(TAREA*sum(DZT_SOMC,3))*1.0E-06 , 'm^3'
   write (*,*) 'Volume of DZT_D:', sum(TAREA*sum(DZT_SOMD,3))*1.0E-06 , 'm^3'
   write (*,*) 'Area   of DZT_A:', sum(TAREA(800:1100,465:745))*1.0E-04 , 'm^2'
   write (*,*) 'Area   of DZT_B:', sum(TAREA(1100:1400,465:745))*1.0E-04 , 'm^2'
   write (*,*) 'Area   of DZT_C:', sum(TAREA(1100:1400,200:465))*1.0E-04 , 'm^2'
   write (*,*) 'Area   of DZT_D:', sum(TAREA(800:1100,200:465))*1.0E-04 , 'm^2'

   do k = 1,km
   if ( k.lt.21 ) then
    DZT_SOMA(:,:,k) = c0
    DZT_SOMB(:,:,k) = c0
    DZT_SOMC(:,:,k) = c0
    DZT_SOMD(:,:,k) = c0
   endif
   enddo

   ! OHC index
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMA,OHC_A_700) ! [J]
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMB,OHC_B_700)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMC,OHC_C_700)
   call vol_int(imt,jmt,km,OHC,TAREA,DZT_SOMD,OHC_D_700)

   ! output
300 FORMAT(A10, ES14.6E2, A6)

   write (*,300) 'OHC        =', OHC_int       , 'J'
   write (*,300) 'OHC_700    =', OHC_int_700   , 'J'
   write (*,300) 'OHC_SS     =', OHC_int_SS    , 'J'
   write (*,300) 'OHC_700_SS =', OHC_int_SS_700, 'J'
   write (*,300) 'OHC_CS     =', OHC_int_CS    , 'J'
   write (*,300) 'OHC_700_CS =', OHC_int_CS_700, 'J'
   write (*,300) 'OHC_NS     =', OHC_int_NS    , 'J'
   write (*,300) 'OHC_700_NS =', OHC_int_NS_700, 'J'
   write (*,300) 'OHC_A      =', OHC_A, 'J'
   write (*,300) 'OHC_B      =', OHC_B, 'J'
   write (*,300) 'OHC_C      =', OHC_C, 'J'
   write (*,300) 'OHC_D      =', OHC_D, 'J'
   write (*,300) 'OHC_A_700  =', OHC_A_700, 'J'
   write (*,300) 'OHC_B_700  =', OHC_B_700, 'J'
   write (*,300) 'OHC_C_700  =', OHC_C_700, 'J'
   write (*,300) 'OHC_D_700  =', OHC_D_700, 'J'
   write (*,300) 'SST_A      =', SST_A, 'deg C'
   write (*,300) 'SST_B      =', SST_B, 'deg C'
   write (*,300) 'SST_C      =', SST_C, 'deg C'
   write (*,300) 'SST_D      =', SST_D, 'deg C'
   write (*,300) 'SHF_A      =', SHF_A, 'W/m^2'
   write (*,300) 'SHF_B      =', SHF_B, 'W/m^2'
   write (*,300) 'SHF_C      =', SHF_C, 'W/m^2'
   write (*,300) 'SHF_D      =', SHF_D, 'W/m^2'

!=====================================================================================================================================================
!  TESTING
!=====================================================================================================================================================

!   write (*,*) ''
!   write (*,*) 'TESTING'

!   write (*,*) 'RHO(1800,1200,1) =', RHO(1800,1200,1)
!   write (*,*) 'PD(1800,1200,1)  =', PD(1800,1200,1)

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
100 FORMAT (23(A,", "),A)
101 FORMAT (23(E15.7,","),E15.7)
   open  (1,file=('OHC_'//year)//'.out',form='formatted')
   write (1,100) 'OHC', 'OHC_SS', 'OHC_CS', 'OHC_NS',                          &
   'OHC_700', 'OHC_700_SS', 'OHC_700_CS', 'OHC_700_NS',                        &
   'OHC_A', 'OHC_B', 'OHC_C', 'OHC_D', 'OHC_A_700', 'OHC_B_700', 'OHC_C_700', 'OHC_D_700', &
   'SST_A', 'SST_B', 'SST_C', 'SST_D', 'SHF_A', 'SHF_B', 'SHF_C', 'SHF_D'
   write (1,101) OHC_int, OHC_int_SS, OHC_int_CS, OHC_int_NS ,                 &
   OHC_int_700, OHC_int_SS_700, OHC_int_CS_700, OHC_int_NS_700,                &
   OHC_A, OHC_B, OHC_C, OHC_D, OHC_A_700, OHC_B_700, OHC_C_700, OHC_D_700,     &
   SST_A, SST_B, SST_C, SST_D, SHF_A, SHF_B, SHF_C, SHF_D
 
   close (1)

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
