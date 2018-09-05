program rKm_budget
implicit none

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  this script calculates various quantities from the original POP output files
!  and the LEC files created with LEC.f90
!
!  calculated quantities:
!  1. vertical integrals of cPKm/cPKe
!  2. Eulerian meridional overturning stream function
!  3. meridional transport of energy terms
!  4. mixed layer depth, vertically integrated w, ohc
!
!  generated output
!  1. two binary files with verticaly integrated cPKm/cPKe for each year 
!     + an average of said years
!  2. two binary file with Eul. stream functions as an output, one calculated
!     from v and one from w
!  3. meridional advection of energy terms
!  4. mixed layer depth, w, ohc
!
!  for curvilinear grid:
!  longitude (i) 0=250E, 1100=0E, 1500=40E
!  latitude  (j) 0=78.5S, 518=55S, 677=45S, n_30S=30S, 1181=0S
!  depth     (k) 1=5m, 9=96.9m, 10=112m, 14=198m, 17=318m, 35=4125, 42=5875m
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!===============================================================================
!  variables
!===============================================================================

! user input 
integer ::                                                                     &
  y, ny, start_year, end_year, nt, ni, nj,                                     &
  nrec_rPm,  nrec_rPe,  nrec_rKm,   nrec_rKe,                                  &
  nrec_gPm,  nrec_gPe,  nrec_gKm,   nrec_gKe,                                  &
  nrec_gPmh, nrec_gPms, nrec_gPeh,  nrec_gPes,                                 &
  nrec_cPem, nrec_cKem, nrec_cPKm,  nrec_cPKe,                                 &
  nrec_TEMP, nrec_VVEL, nrec_WVEL,                                             &
  nrec_HMXL, nrec_XMXL, nrec_TMXL,  nrec_RHO,                                  &
  nrec_UVEL, nrec_SSH,  nrec_VVEL2, nrec_UV,   nrec_PD

character*120 ::                                                               & 
  grid_file,kmt_file, in_depths,      pbc_file,                                &
  geometry1_file,     geometry2_file, input_folder, filename3
character*38  :: LEC_file
character*58  :: bin_file
character*3   :: yr
character*62  :: filename0,filename1
character*42  :: filename2

! internal model variables 
integer                                         :: rec_length,i,j,k
integer,          dimension(:,:),   allocatable :: KMU, NA_mask
integer,          dimension(:,:,:), allocatable :: NA_mask3D

real, dimension(:),     allocatable :: dz,area
real, dimension(:,:),   allocatable :: DXT, DYT, TAREA, DXU, DYU, UAREA

real                                            ::                             &
  rPm_int, rPe_int, rKm_int, rKe_int, cPKm_yint, cPKe_yint, cKem_yint,         &
  vzint_VVEL, vzint_TTTV, pressure,                                            &
  rKm_SO30, gKm_SO30, cPKm_SO30, cKem_SO30, bKm4, bKm5,vrKm_SO30,              &
  udpdx_SO30, vdpdy_SO30, vp_SO30, uAu_SO30, vAv_SO30,                         &
  cPKm_old_SO30, PUx_SO30, PVy_SO30, PWz_SO30, WPz_SO30, Pz_SO30, gPD_SO30,    &
  UPx_SO30, VPy_SO30,                                                          &
  ddxPU_SO30,ddyPV_SO30,ddzPW_SO30,cPKm_j,VPy_j,PVy_j,cPKm_jA,VPy_jA,PVy_jA,   &
  PD_SO30,PD_NO30,WVEL_SO30,WVEL_NO30,SSH_SO30,SSH_NO30,TAREA_SO30,TAREA_NO30, &
  DZT_SO30,DZT_NO30,                                                           &
  VPy_SA30, PVy_SA30, UPx_SA30, PUx_SA30,                                      &
  VPy_NA30, PVy_NA30, UPx_NA30, PUx_NA30

real,             dimension(:),     allocatable ::                             &
  tdepth, cPKm_zint, cPKe_zint, cKem_zint,                                     &
  p_prof_glob, p_prof_so, p_prof_fun, vp_glob_zint, vp_so_zint, vp_fun_zint

real,             dimension(:,:),   allocatable ::                             &
  cPKm_vint, cPKe_vint, cKem_vint, gKm,                                        &
  vrKm,                                                                        &
  ref_state, vp_glob_vint, vp_so_vint, vp_fun_vint, SSH,                       &
  WORK1,DUS,DUN,DUE,DUW,KXU,KYU,WORK2,DXKX,DYKX,DYKY,DXKY,                     &
  DUM,DMC,DME,DMN,DUC,DMW,DMS,CC,CN,CS,CE,CW,D2UK,D2VK,AMF,                    &
  HUS,HTE,HUW,HTN, DYUU, DXUV

real,             dimension(:,:,:), allocatable ::                             &
  DZT, DZU, cPKm, cPKe, cKem, rPm, rPe, rKm, rKe, VVEL, WVEL, TTT_VVEL,        &
  FBC, HSP, HSP1, HSP2, HSP3, UVEL, VVEL2, UV, DPDX, DPDY, HDUK, HDVK,         &
  cPKm_new, PD, cPKm_old, PUx, PVy, PWz, WPz,                                  &
  UPx, VPy, UPx4, VPy4

integer,parameter ::                                                           &
  spy = 24*3600*360, n_30S=866, imt = 3600, jmt = 2400, km = 42

double precision, parameter ::                                                 &
  c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1., c2=2.0,                &
  S0ref = 0.035, rho0 = 4.1/3.996*1000, c = 3996,                              &
  g = 9.806, salinity_factor = -34.7*1.0E-04, RE = 6.37E06,                    &
  am = -27.0e17 !from pop_in file


write (*,*) ''
write (*,*) '--- rKm BUDGET TERMS ANALYSIS ---'

!===============================================================================
!  INPUT
!===============================================================================



input_folder    = '/home/dijkbio/andre/LEC/input/'
in_depths       = trim(input_folder)//'in_depths.42.dat'
geometry1_file  = trim(input_folder)//'geometry1'
geometry2_file  = trim(input_folder)//'geometry2'

! read geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km),         &
          HUS(imt,jmt), HTE(imt,jmt), HUW(imt,jmt),   HTN(imt,jmt),            &
          dz(km),       area(km),     ref_state(6,km),FBC(imt,jmt,km) )
open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,    &
       status='old')
open(2,file=geometry2_file,access='direct',form='unformatted',recl=6,          &
       status='old')
! k, dz, tdepth, area, p,
! <T>, <S>, <RHO>, <PD>, <D(D,T,0)>, <D(D,T,p)>, <Q>,
! dz(<PD>), dz(<D(S,T,0)>), dz(<D(S,T,p)>)
read(1,rec= 1) DXT ! [m]
read(1,rec= 2) DXU ! [m]
read(1,rec= 3) DYT ! [m]
read(1,rec= 4) DYU ! [m]
read(1,rec= 5) TAREA ! [m^2]
read(1,rec= 6) UAREA ! [m^2]
read(1,rec=91) HTN ! [m]
read(1,rec=92) HTE ! [m]
read(1,rec=93) HUS ! [m]
read(1,rec=94) HUW ! [m]
do k=1,km
  read(1,rec=6+   k) DZT(:,:,k) ! [m]
  read(1,rec=6+km+k) DZU(:,:,k) ! [m]
  read(2,rec=k)      ref_state(:,k) ! [varying]
  dz(k)            = ref_state(2,k) ! [m]
  area(k)          = ref_state(4,k) ! [m^2]
  FBC(:,:,k)       = dz(k)
enddo

close(1)
close(2)

nrec_rPm   =    1
nrec_rPe   =   43
nrec_rKm   =   85
nrec_rKe   =  127
nrec_gPm   =  169
nrec_gPe   =  170
nrec_gKm   =  171
nrec_gKe   =  172
nrec_gPmh  =  173
nrec_gPms  =  174
nrec_gPeh  =  175
nrec_gPes  =  176
nrec_cPem  =  177
nrec_cKem  =  219
nrec_cPKm  =  261
nrec_cPKe  =  303

nrec_VVEL  =   43
nrec_WVEL  =  405
nrec_RHO   =  306
nrec_PD    =  995
nrec_UVEL  =    1
nrec_SSH   =  390
nrec_VVEL2 =  127
nrec_UV    =  348

!===============================================================================
!  CALCULATIONS
!===============================================================================

! years
start_year = 278
end_year   = 324
nt = end_year-start_year+1

! avg binary file
filename0 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.5.year.avg'  
open (1,file=trim(filename0),access='direct',form='unformatted',recl=imt*jmt,  &
        status='unknown')

! output file
filename3 = '../../results/rKm_budget/rKm_budget.out'
open (3,file=trim(filename3),form='formatted')

!===============================================================================
!  pressure profiles
!===============================================================================
allocate( HSP(imt,jmt,km), SSH(imt,jmt), PD(imt,jmt,km), p_prof_glob(km) )
call load_3D_field(1,nrec_PD,PD)
call load_2D_field(1,nrec_SSH,SSH)
call hydrostatic_pressure(g,dz,PD,SSH,HSP)
do k=1,km
  p_prof_glob(k) = sum(HSP(:,:,k)    *TAREA(:,:)    ,DZT(:,:,k).ne.0.0)        &
                 /sum(TAREA(:,:),DZT(:,:,k).ne.0.0 )
enddo !k
close(1)

!===============================================================================
!  horizontal viscosity fields
!===============================================================================
!allocate(CC(imt,jmt)  ,CN(imt,jmt)  ,CS(imt,jmt)  ,CE(imt,jmt)  ,CW(imt,jmt),  &
!         DMC(imt,jmt) ,DMN(imt,jmt) ,DMS(imt,jmt) ,DME(imt,jmt) ,DMW(imt,jmt), &
!         D2UK(imt,jmt),D2VK(imt,jmt),HDUK(imt,jmt,km),HDVK(imt,jmt,km),        &
!         AMF(imt,jmt))
!
!AMF(:,:) = (UAREA(:,:)/UAREA(1,1181))**1.5 ! (1,1181) at equator
!
!, HSP2(imt,jmt,km), HSP3(imt,jmt,km),              &
!          vp_glob_vint(imt,jmt),
!allocate( vp_so_vint(imt,jmt))!, vp_fun_vint(imt,jmt),  &
!          vp_glob_zint(jmt), 
!allocate( vp_so_zint(jmt) )!, vp_fun_zint(jmt),              &
!allocate(WORK1(imt,jmt),WORK2(imt,jmt),KMU(imt,jmt),                           &
!         DUS(imt,jmt),DUN(imt,jmt),DUE(imt,jmt),DUW(imt,jmt),                  &
!         KXU(imt,jmt),KYU(imt,jmt),                                            &
!         DXKX(imt,jmt),DYKX(imt,jmt),DYKY(imt,jmt),DXKY(imt,jmt),              &
!         DUM(imt,jmt),DUC(imt,jmt) )
!
!KMU(:,:) = 0
!do k=1,km
!  where (DZU(:,:,k)>0.0)
!    KMU(:,:) = KMU(:,:) + 1
!  endwhere
!enddo !k
!
!
!ni=500
!nj=500
!
!DUS = WORK1/UAREA
!DUN = eoshift(WORK1,dim=2,shift=1)/UAREA
!
!WORK1 = HUW/HTN
!DUW = WORK1/UAREA
!DUE = cshift(WORK1,dim=1,shift=1)/UAREA
!
!KXU = ( cshift(HUW,dim=1,shift=1) - HUW)/UAREA
!KYU = (eoshift(HUS,dim=2,shift=1) - HUS)/UAREA
!
!WORK1 = (HTE - cshift(HTE,dim=1,shift=-1))/TAREA  ! KXT
!WORK2 = cshift(WORK1,dim=1,shift=1) - WORK1
!DXKX  = p5*(WORK2 + eoshift(WORK2,dim=2,shift=1))/DXU
!WORK2 = eoshift(WORK1,dim=2,shift=1) - WORK1
!DYKX  = p5*(WORK2 + cshift(WORK2,dim=1,shift=1))/DYU
!WORK1 = (HTN - eoshift(HTN,dim=2,shift=-1))/TAREA  ! KYT
!WORK2 = eoshift(WORK1,dim=2,shift=1) - WORK1
!DYKY  = p5*(WORK2 + cshift(WORK2,dim=1,shift=1))/DYU
!WORK2 = cshift(WORK1,dim=1,shift=1) - WORK1
!DXKY = p5*(WORK2 + eoshift(WORK2,dim=2,shift=1))/DXU
!
!DUM = -(DXKX + DYKY + c2*(KXU**2 + KYU**2))
!DMC = DXKY - DYKX
!
!DME =  c2*KYU/(HTN + cshift(HTN,dim=1,shift=1))
!DMN = -c2*KXU/(HTE + eoshift(HTE,dim=2,shift=1))
!
!!write(*,*) DUN(ni,nj),DUS(ni,nj),DUE(ni,nj),DUW(ni,nj)
!DUC = -(DUN + DUS + DUE + DUW)               ! scalar laplacian
!DMW = -DME
!DMS = -DMN
!CC = DUC + DUM
!
!
!deallocate(HUS,HTE,HUW,HTN)
!deallocate(WORK1,KXU,KYU,WORK2,DXKX,DYKX,DYKY,DXKY)

allocate( VVEL(imt,jmt,km),     UVEL(imt,jmt,km),                              &
          DPDX(imt,jmt,km),     DPDY(imt,jmt,km),                              &
          cPKm_new(imt,jmt,km), WVEL(imt,jmt,km),                              &
          cPKm(imt,jmt,km),     vrKm(imt,jmt),        rKm(imt,jmt,km),         &
          PUx(imt,jmt,km),      PVy(imt,jmt,km),  PWz(imt,jmt,km),             &
          WPz(imt,jmt,km),      UPx(imt,jmt,km), VPy(imt,jmt,km),              &
          DYUU(imt,jmt),        DXUV(imt,jmt),                                 &
          UPx4(imt,jmt,km),     VPy4(imt,jmt,km) )

! creating a mask for   
allocate( NA_mask(imt,jmt), NA_mask3D(imt,jmt,km) )
NA_mask(   :    ,    :   ) = 0
write(*,*) sum(NA_mask)
NA_mask(401:1400,   1:jmt) = 1
NA_mask(267: 400,1268:jmt) = 1
NA_mask(167: 266,1338:jmt) = 1
NA_mask( 85: 166,1364:jmt) = 1
do k=1,km
  NA_mask3D(:,:,k) = NA_mask(:,:)
enddo !k
!===============================================================================
!  Looping through years
!===============================================================================

write (*,*) ''
100 FORMAT(A20, ES14.6E2)
ny=1

!do y = 290,315,25
do y = start_year, end_year
   write(*,*) 'year', y
  ! open file
  write(yr,"(I3)") y
  ! original binary POP output file
  filename1 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.5.year.'//yr
  write(*,*) filename1
  open (1,file=trim(filename1),access='direct',form='unformatted',recl=imt*jmt,&
        status='unknown')
  ! LEC.f90 output file
  filename2 = '/projects/0/samoc/jan/Andree/LEC_bin_5_'//yr  
  open (2,file=trim(filename2),access='direct',form='unformatted',recl=imt*jmt,&
        status='unknown')

  call load_3D_field(1,nrec_PD  ,  PD)
  call load_2D_field(1,nrec_SSH , SSH)
  call load_3D_field(1,nrec_VVEL,VVEL)
  call load_3D_field(1,nrec_UVEL,UVEL)
  call load_3D_field(1,nrec_WVEL,WVEL)
  call load_3D_field(2,nrec_cPKm,cPKm)
  call load_3D_field(2,nrec_rKm , rKm)

  !=============================================================================
  !  1. Testing volume transport across latitudes
  !=============================================================================

!  call vol_int_part(1,1,1,imt,n_30S,km,rKm,UAREA,dz,DZU,rKm_SO30)
!  call surf_int_2D(1,1,imt,n_30S,gKm,UAREA,DZU(:,:,1),gKm_SO30) 
!
!  write(*,100) 'rKm                :', rKm_SO30
!  write(*,100) 'gKm                :', gKm_SO30

  !if (y==start_year) then
  !write(*,*) 'j=',j 
  !vzint_VVEL = 0.0
  !vzint_TTTV = 0.0
  !do i=1,imt
  !  do k=1,km
  !     vzint_VVEL = vzint_VVEL + DXU(i,j) * DZU(i,j,k) *     VVEL(i,j,k)
       !vzint_TTTV = vzint_TTTV + DXT(i,j) * DZT(i,j,k) * TTT_VVEL(i,j,k)
  !  enddo !k
  !enddo !i

  !write(*,*) 'vzint_VVEL = ', vzint_VVEL/1e2, 'm^3/s'
  !write(*,*) 'vzint_TTTV = ', vzint_TTTV/1e2, 'm^3/s'
  !write(*,*) ' '
  !endif

  !=============================================================================
  !  2. cPKm/cPKe volume integral to latitude
  !=============================================================================

  do k = 1, km
    if (k==1) then
      cPKm_new(:,:,k) = ( WVEL(:,:,k+1)*(PD(:,:,k)+PD(:,:,k+1)) )
    elseif (k==km) then
      cPKm_new(:,:,k) = ( WVEL(:,:,k  )*(PD(:,:,k)+PD(:,:,k-1)) )
    else
      cPKm_new(:,:,k) = ( WVEL(:,:,k+1)*(PD(:,:,k)+PD(:,:,k+1))           &
                        + WVEL(:,:,k  )*(PD(:,:,k)+PD(:,:,k-1)) )
    endif
    cPKm_new(:,:,k) = -g * p25 * cPKm_new(:,:,k) * 1.0E01
  enddo !k
  call vol_int_part(1,1,1,imt,n_30S,km,cPKm,TAREA,dz,FBC,cPKm_SO30)
  write(*,100) 'cPKm_old           :', cPKm_SO30
  call vol_int_part(1,1,1,imt,n_30S,km,cPKm_new,TAREA,dz,FBC,cPKm_SO30)
  write(*,100) 'cPKm_new           :', cPKm_SO30

  !=============================================================================
  !  3. Meridional Advection Of Energy VVEL*rE
  !=============================================================================
  
  ! calculate the acdvection (vertically integrated and weighted by UAREA)
  vrKm(:,:) = sum(VVEL(:,:,:)*1.0E-02*rKm(:,:,:)*DZU(:,:,:),dim=3)*UAREA(:,:)
  vrKm(:,:) = p25*(              (vrKm(:,:)               )                +   &
                           cshift(vrKm(:,:),dim=1,shift=-1)                +   &
                   eoshift(       vrKm(:,:),dim=2,shift=-1)                +   &
                   eoshift(cshift(vrKm(:,:),dim=1,shift=-1),dim=2,shift=-1))/ &
                TAREA(:,:)
  write (*,100) 'v*rKm(n_30S)       : ', sum(vrKm(:,n_30S)*DXT(:,n_30S))



  !=============================================================================
  ! 4. Pressure work term (VVEL*HSP_anomaly)
  !=============================================================================

  ! calculate hydrostatic pressure everywhere
  call hydrostatic_pressure(g,dz,PD,SSH,HSP)
  call gradient(DXU,DYU,HSP,DPDX,DPDY)
  
  ! UPx, VPy, WPy
  do k=1,km
    UPx(:,:,k) = UVEL(:,:,k)*DPDX(:,:,k)*1.0E-02*                              &
                 (UAREA(:,:)*DZU(:,:,k))
    VPy(:,:,k) = VVEL(:,:,k)*DPDY(:,:,k)*1.0E-02*                              &
                 (UAREA(:,:)*DZU(:,:,k))
    UPx4(:,:,k) = p25*(         (UPx(:,:,k)               )                +   &
                          cshift(UPx(:,:,k),dim=1,shift=-1)                +   &
                  eoshift(       UPx(:,:,k),dim=2,shift=-1)                +   &
                  eoshift(cshift(UPx(:,:,k),dim=1,shift=-1),dim=2,shift=-1))/  &
                  (TAREA(:,:)*FBC(:,:,k))
  
    VPy4(:,:,k) = p25*(         (VPy(:,:,k)               )                +   &
                          cshift(VPy(:,:,k),dim=1,shift=-1)                +   &
                  eoshift(       VPy(:,:,k),dim=2,shift=-1)                +   &
                  eoshift(cshift(VPy(:,:,k),dim=1,shift=-1),dim=2,shift=-1))/  &
                  (TAREA(:,:)*FBC(:,:,k))

    if (k==km) then
      WPz(:,:,k) = WVEL(:,:,k  )*(HSP(:,:,k)  -HSP(:,:,k-1))/(dz(k)+dz(k-1))
    elseif (k==1) then
      WPz(:,:,k) = WVEL(:,:,k+1)*(HSP(:,:,k+1)-HSP(:,:,k  ))/(dz(k)+dz(k+1))
    else
      WPz(:,:,k) = WVEL(:,:,k+1)*(HSP(:,:,k+1)-HSP(:,:,k  ))/(dz(k)+dz(k+1)) &
                 + WVEL(:,:,k  )*(HSP(:,:,k)  -HSP(:,:,k-1))/(dz(k)+dz(k-1))
    endif
  PWz(:,:,k) = PWz(:,:,k)*1.0E-02
  WPz(:,:,k) = WPz(:,:,k)*1.0E-02
  enddo !k


  ! PUx, PVy, PWz
  do k=1,km
    DYUU(:,:)  = DYU(:,:)*DZU(:,:,k)*UVEL(:,:,k)*1.0E-02/dz(k) 
    PUx(:,:,k) = HSP(:,:,k)*p5/TAREA(:,:)*(                                  &
                                DYUU(:,:)                                    &
                +eoshift(       DYUU(:,:),dim=2,shift=-1)                    &
                -        cshift(DYUU(:,:),dim=1,shift=-1)                    &
                -eoshift(cshift(DYUU(:,:),dim=1,shift=-1),dim=2,shift=-1))

    DXUV(:,:)  = DXU(:,:)*DZU(:,:,k)*VVEL(:,:,k)*1.0E-02/dz(k) 
    PVy(:,:,k) = HSP(:,:,k)*p5/TAREA(:,:)*(                                  &
                                DXUV(:,:)                                    &
                -eoshift(       DXUV(:,:),dim=2,shift=-1)                    &
                +        cshift(DXUV(:,:),dim=1,shift=-1)                    &
                -eoshift(cshift(DXUV(:,:),dim=1,shift=-1),dim=2,shift=-1))

    if (k==km) then
      PWz(:,:,k) = HSP(:,:,k)*(WVEL(:,:,k))*1.0E-02/dz(k)
    else
      PWz(:,:,k) = HSP(:,:,k)*(WVEL(:,:,k)-WVEL(:,:,k+1))*1.0E-02/dz(k)
    endif
  enddo !k


! OUTPUT
  call vol_int_part(1,1,1,imt,n_30S,km,PUx,TAREA,dz,FBC,PUx_SO30)
  call vol_int_part(1,1,1,imt,n_30S,km,PVy,TAREA,dz,FBC,PVy_SO30)
  call vol_int_part(1,1,1,imt,n_30S,km,PWz,TAREA,dz,FBC,PWz_SO30)

  call vol_int_part(1,1,1,imt,n_30S,km,UPx4,TAREA,dz,FBC,UPx_SO30)
  call vol_int_part(1,1,1,imt,n_30S,km,VPy4,TAREA,dz,FBC,VPy_SO30)
  call vol_int_part(1,1,1,imt,n_30S,km, WPz,TAREA,dz,FBC,WPz_SO30)

  write(*,  *) 'Integrals to 30S'
  write(*,100) 'UPx                :', UPx_SO30 
  write(*,100) '    PUx            :',            PUx_SO30
  write(*,100) 'UPx+PUx            :', UPx_SO30 + PUx_SO30
  write(*,100) 'VPy                :', VPy_SO30 
  write(*,100) '    PVy            :',            PVy_SO30
  write(*,100) 'VPy+PVy            :', VPy_SO30 + PVy_SO30
  write(*,100) 'WPz                :', WPz_SO30
  write(*,100) '    PWz            :',            PWz_SO30
  write(*,100) 'WPz-PWz            :', WPz_SO30 - PWz_SO30
  write(*,100) 'PUx+PVy+PWz        :', PUx_SO30 + PVy_SO30 + PWz_SO30

  ! Southern Atlantic [70W, 20E] south of 30S
  call vol_int_part(400,1,1,1300,n_30S,km,VPy4,TAREA,dz,FBC,VPy_SA30)
  call vol_int_part(400,1,1,1300,n_30S,km,PVy ,TAREA,dz,FBC,PVy_SA30)
  call vol_int_part(400,1,1,1300,n_30S,km,UPx4,TAREA,dz,FBC,UPx_SA30)
  call vol_int_part(400,1,1,1300,n_30S,km,PUx ,TAREA,dz,FBC,PUx_SA30)

  ! successfully tested whether intragel over remaining Southern Ocean is the
  ! same as the difference between the SO30 and SA30  

  ! North Atlantic (masked) north of 30S
  call vol_int_part(1,1,1,imt,jmt,km,NA_mask3D*VPy4,TAREA,dz,FBC,VPy_NA30)
  call vol_int_part(1,1,1,imt,jmt,km,NA_mask3D*PVy ,TAREA,dz,FBC,PVy_NA30)
  call vol_int_part(1,1,1,imt,jmt,km,NA_mask3D*UPx4,TAREA,dz,FBC,UPx_NA30)
  call vol_int_part(1,1,1,imt,jmt,km,NA_mask3D*PUx ,TAREA,dz,FBC,PUx_NA30)

  !=============================================================================
  ! 5. Biharmonic friction terms (VVEL*A4*D3VDY3)
  !=============================================================================

  ! mistake in Wu et al.
  ! this term does not exist, rather, the horizontal friction is calculated 

  !=============================================================================
  ! 6. Eddy advection terms ()
  !=============================================================================


!  bKm4 = 0.0
!  bKm5 = 0.0
!  do k=1,km
!    do i=1,imt
!      bKm4 = bKm4 + (UV(i,n_30S,k)    - UVEL(i,n_30S,k)*VVEL(i,n_30S,k) )      &
!                    * UVEL(i,n_30S,k) * DZU(i,n_30S,k) * DXU(i,n_30S)
!      bKm5 = bKm5 + (VVEL2(i,n_30S,k) - VVEL(i,n_30S,k)*VVEL(i,n_30S,k) )      &
!                    * VVEL(i,n_30S,k) * DZU(i,n_30S,k) * DXU(i,n_30S)
!    enddo !i
!  enddo !k
!  bKm4 = bKm4 * 1.0E-06 * rho0
!  bKm5 = bKm5 * 1.0E-06 * rho0
!
!  write(*,100) 'bKm4               :', bKm4 
!  write(*,100) 'bKm5               :', bKm5



  !=============================================================================
  ! 7. Horizontal dissipation
  !=============================================================================
  
  ! from hmix_del4.F90 of POP source code
  ! partial_bottom_cells = .true.
  ! ltopostress          = .false.
  ! lvariable_hmixu      = .true.
  ! changed A(i+1) to cshift(A,dim=1,shift=1) and in j direction with eoshift

! -----------------------------------------------------------------------------
!  do k=1,km
!    CN = DUN*min(eoshift(DZU(:,:,k),dim=2,shift= 1),DZU(:,:,k))/DZU(:,:,k)
!    CS = DUS*min(eoshift(DZU(:,:,k),dim=2,shift=-1),DZU(:,:,k))/DZU(:,:,k)
!    CE = DUE*min( cshift(DZU(:,:,k),dim=1,shift= 1),DZU(:,:,k))/DZU(:,:,k)
!    CW = DUW*min( cshift(DZU(:,:,k),dim=1,shift=-1),DZU(:,:,k))/DZU(:,:,k)
!
!    D2UK = (CC *        UVEL(:,:,k)                 +                          &
!            CN *eoshift(UVEL(:,:,k),dim=2,shift= 1) +                          &
!            CS *eoshift(UVEL(:,:,k),dim=2,shift=-1) +                          &
!            CE * cshift(UVEL(:,:,k),dim=1,shift= 1) +                          &
!            CW * cshift(UVEL(:,:,k),dim=1,shift=-1))+                          &
!           (DMC*        VVEL(:,:,k)                 +                          &
!            DMN*eoshift(VVEL(:,:,k),dim=2,shift= 1) +                          &
!            DMS*eoshift(VVEL(:,:,k),dim=2,shift=-1) +                          &
!            DME* cshift(VVEL(:,:,k),dim=1,shift= 1) +                          &
!            DMW* cshift(VVEL(:,:,k),dim=1,shift=-1))
!    D2VK = (CC *        VVEL(:,:,k)                 +                          &
!            CN *eoshift(VVEL(:,:,k),dim=2,shift= 1) +                          &
!            CS *eoshift(VVEL(:,:,k),dim=2,shift=-1) +                          &
!            CE * cshift(VVEL(:,:,k),dim=1,shift= 1) +                          &
!            CW * cshift(VVEL(:,:,k),dim=1,shift=-1))-                          &
!           (DMC*        UVEL(:,:,k)                 +                          &
!            DMN*eoshift(UVEL(:,:,k),dim=2,shift= 1) +                          &
!            DMS*eoshift(UVEL(:,:,k),dim=2,shift=-1) +                          &
!            DME* cshift(UVEL(:,:,k),dim=1,shift= 1) +                          &
!            DMW* cshift(UVEL(:,:,k),dim=1,shift=-1))
!    if (k==1) then
!!      write(*,*) CC(ni,nj),CN(ni,nj),CS(ni,nj),CE(ni,nj),CW(ni,nj)      !0
!!      write(*,*) DUN(ni,nj),DUS(ni,nj),DUE(ni,nj),DUW(ni,nj)      !0
!!      write(*,*) DMC(ni,nj),DMN(ni,nj),DMS(ni,nj),DME(ni,nj),DMW(ni,nj)      !0
!!      write(*,*) DZU(ni,nj,k),DZU(ni+1,nj,k),DZU(ni-1,nj,k),DZU(ni,nj+1,k),DZU(ni,nj-1,k)      !0
!!      write(*,*) UVEL(ni,nj,k),UVEL(ni+1,nj,k),UVEL(ni-1,nj,k),UVEL(ni,nj+1,k),UVEL(ni,nj-1,k) !0
!!      write(*,*)  DUE(ni,nj)*min(DZU(ni+1,nj,k),DZU(ni,nj,k))/DZU(ni,nj,k)
!    endif
!
!    where (k <= KMU)
!      D2UK = AMF*D2UK
!      D2VK = AMF*D2VK
!    elsewhere
!      D2UK = c0
!      D2VK = c0
!    end where
!
!    HDUK(:,:,k) = am*((CC *        D2UK                 +                      &
!                       CN *eoshift(D2UK,dim=2,shift= 1) +                      &
!                       CS *eoshift(D2UK,dim=2,shift=-1) +                      &
!                       CE * cshift(D2UK,dim=1,shift= 1) +                      &
!                       CW * cshift(D2UK,dim=1,shift=-1))+                      &
!                      (DMC*        D2VK                 +                      &
!                       DMN*eoshift(D2VK,dim=2,shift= 1) +                      &
!                       DMS*eoshift(D2VK,dim=2,shift=-1) +                      &
!                       DME* cshift(D2VK,dim=1,shift= 1) +                      &
!                       DMW* cshift(D2VK,dim=1,shift=-1)))
!    HDVK(:,:,k) = am*((CC *        D2VK                 +                      &
!                       CN *eoshift(D2VK,dim=2,shift= 1) +                      &
!                       CS *eoshift(D2VK,dim=2,shift=-1) +                      &
!                       CE * cshift(D2VK,dim=1,shift= 1) +                      &
!                       CW * cshift(D2VK,dim=1,shift=-1))-                      &
!                      (DMC*        D2UK                 +                      &
!                       DMN*eoshift(D2UK,dim=2,shift= 1) +                      &
!                       DMS*eoshift(D2UK,dim=2,shift=-1) +                      &
!                       DME* cshift(D2UK,dim=1,shift= 1) +                      &
!                       DMW* cshift(D2UK,dim=1,shift=-1)))
!
!    where (k > KMU)
!      HDUK(:,:,k) = c0
!      HDVK(:,:,k) = c0
!    endwhere
!  enddo !k
!
!  ! factors: (1e-2)^2 for two velocities from cm/s to m/s
!  !          1e-8 because i.s.o. dividing by cm^4 we divide by m^4 in nabla
!  !          terms
!  write(*,100) 'u d2 A_m d2 u      :', uAu_SO30
!  uAu_SO30 = sum(sum(UVEL(:,1:n_30S,:)*HDUK(:,1:n_30S,:)*                      &
!                        DZU(:,1:n_30S,:),dim=3)*UAREA(:,1:n_30S))              &
!                 *rho0*1.0E-12
!  write(*,100) 'v d2 A_m d2 v      :', vAv_SO30
!  vAv_SO30 = sum(sum(VVEL(:,1:n_30S,:)*HDVK(:,1:n_30S,:)*                      &
!                        DZU(:,1:n_30S,:),dim=3)*UAREA(:,1:n_30S))              &
!                 *rho0*1.0E-12
 write(*,*) ' '

  close(1)
  close(2)
  
  !=============================================================================
  ! 8. integrated WVEL and PD
  !=============================================================================

  call vol_int_part(1,      1,1,imt,n_30S,km,  PD,TAREA,dz,FBC,   PD_SO30)
  call vol_int_part(1,      1,1,imt,n_30S,km,WVEL,TAREA,dz,FBC, WVEL_SO30)
  call vol_int_part(1,n_30S+1,1,imt,  jmt,km,  PD,TAREA,dz,FBC,   PD_NO30)
  call vol_int_part(1,n_30S+1,1,imt,  jmt,km,WVEL,TAREA,dz,FBC, WVEL_NO30)
  call surf_int_2D(1,      1,imt,n_30S,SSH,TAREA,DZT(:,:,1),SSH_SO30)
  call surf_int_2D(1,n_30S+1,imt,jmt  ,SSH,TAREA,DZT(:,:,1),SSH_NO30)
  TAREA_SO30  = sum(TAREA,DZT(:,      1:n_30S,1).ne.0.0)
  TAREA_NO30  = sum(TAREA,DZT(:,n_30S+1:  jmt,1).ne.0.0)
  DZT_SO30    = sum(DZT(:,      1:n_30S,:))
  DZT_NO30    = sum(DZT(:,n_30S+1:  jmt,:))
  !=============================================================================
  ! OUTPUT
  !=============================================================================



300 FORMAT (25(A,","),A)
301 FORMAT (25(E15.7,","),E15.7)
  if (ny==1) then 
    write(3,300) 'year','cPKm_SO30','UPx_SO30', 'PUx_SO30', 'VPy_SO30',        &
                 'PVy_SO30', 'WPz_SO30', 'PWz_SO30',                           &
                 'PD_SO30', 'PD_NO30', 'WVEL_SO30', 'WVEL_NO30',               &
                 'SSH_SO30', 'SSH_NO30', 'TAREA_SO30', 'TAREA_NO30',           &
                 'DZT_SO30', 'DZT_NO30',                                       &
                 'VPy_SA30', 'PVy_SA30', 'UPx_SA30', 'PUx_SA30',               &
                 'VPy_NA30', 'PVy_NA30', 'UPx_NA30', 'PUx_NA30'
  endif
  write(3,301) real(y), cPKm_SO30,                                             &
                  UPx_SO30, PUx_SO30, VPy_SO30, PVy_SO30, WPz_SO30, PWz_SO30,  &
                  PD_SO30, PD_NO30, WVEL_SO30, WVEL_NO30,                      &
                  SSH_SO30, SSH_NO30, TAREA_SO30, TAREA_NO30,DZT_SO30,DZT_NO30,&
                  VPy_SA30, PVy_SA30, UPx_SA30, PUx_SA30,                      &
                  VPy_NA30, PVy_NA30, UPx_NA30, PUx_NA30

  if (y==290 .or. y==315) then
200 FORMAT (3(A,","),A)
201 FORMAT (3(E15.7,","),E15.7)
    open (4,file=trim('../../results/rKm_budget/lat_dependence_'//yr),         &
          form='formatted')
    write(4,200) 'j','cPKm_j','VPy_j','PVy_j'!,'cPKm_jA','VPy_jA','PVy_jA'
    cPKm_j  = 0.0
    VPy_j   = 0.0
    PVy_j   = 0.0
!    cPKm_jA = 0.0
!    VPy_jA  = 0.0
!    PVy_jA  = 0.0
    do j=1,jmt
      cPKm_j  = cPKm_j                                                         & 
              + sum(TAREA(:,j)*sum(cPKm(:,j,:)*FBC(:,j,:),2,DZT(:,j,:).ne.0.0))
      VPy_j   = VPy_j                                                          &
              + sum(TAREA(:,j)*sum(VPy4(:,j,:)*FBC(:,j,:),2,DZT(:,j,:).ne.0.0))
      PVy_j   = PVY_j                                                          &
              + sum(TAREA(:,j)*sum( PVy(:,j,:)*FBC(:,j,:),2,DZT(:,j,:).ne.0.0))
       
!      cPKm_jA = cPKm_jA                                                        & 
!              + sum(NA_mask(:,j)*TAREA(:,j)*sum(cPKm(:,j,:)*FBC(:,j,:),2,DZT(:,j,:).ne.0.0))
!      VPy_jA  = VPy_jA                                                         &
!              + sum(NA_mask(:,j)*TAREA(:,j)*sum(VPy4(:,j,:)*FBC(:,j,:),2,DZT(:,j,:).ne.0.0))
!      PVy_jA  = PVY_jA                                                         &
!              + sum(NA_mask(:,j)*TAREA(:,j)*sum( PVy(:,j,:)*FBC(:,j,:),2,DZT(:,j,:).ne.0.0))

      write(4,201) real(j), cPKm_j, VPy_j, PVy_j
    enddo
    close(4)
  endif


  ny = ny+1
enddo ! y



!===============================================================================
!  SUBROUTINES
!===============================================================================

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


subroutine read_yrly(x,y,ntavg,nfile,FIELD)
implicit none
!
!  reads yearly data that was previoously calculated and written into a binary
!  file
!

integer,                       intent(in)  :: x,y,ntavg,nfile
integer                                    :: t
real,    dimension(x,y,ntavg), intent(out) :: FIELD

do t = 1,ntavg
  read (nfile,rec=t) FIELD(:,:,t)
enddo

end subroutine read_yrly


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine average(x,y,ntavg,FIELD_yrly,AVG_FIELD)
!
!  meridional advection of a energy reservoir term
!
implicit none

integer,                       intent(in)  :: x,y,ntavg
integer                                    :: t
real,    dimension(x,y,ntavg), intent(in)  :: FIELD_yrly
real,    dimension(x,y),       intent(out) :: AVG_FIELD

AVG_FIELD(:,:) = 0.0
do t=1,ntavg
  AVG_FIELD(:,:) = AVG_FIELD(:,:) + 1.0/ntavg*FIELD_yrly(:,:,t)
enddo

end subroutine


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine mer_advection(DXT,DZT,VVEL,Eres,Eres_adv)
!  meridional advection of a energy reservoir term
implicit none

integer, parameter                          :: imt=3600,jmt=2400,km=42
integer                                     :: i,j,k
real,    dimension(imt,jmt,km), intent(in)  :: DZT, VVEL,Eres
real,    dimension(imt,jmt),    intent(in)  :: DXT
real,    dimension(imt,jmt),    intent(out) :: Eres_adv

!Eres_adv(:,:) = sum(Eres(:,:,:)*TTT_VVEL(:,:,:)*DZT(:,:,:),3)/1.e2 
do i=2,imt
  do j=2,jmt-1
    Eres_adv(i,j)=0.
    do k=1,km
    Eres_adv(i,j) = Eres_adv(i,j) +&
0.25*(Eres(i,j,k)+Eres(i,j+1,k))*(VVEL(i,j,k)*DZT(i,j,k)+VVEL(i-1,j,k)*DZT(i-1,j,k))/1.e2
    enddo
  enddo
enddo
i=1
  do j=2,jmt-1
    Eres_adv(1,j)=0.
    do k=1,km
    Eres_adv(i,j) = Eres_adv(i,j) +&
0.25*(Eres(i,j,k)+Eres(i,j+1,k))*(VVEL(i,j,k)*DZT(i,j,k)+VVEL(imt,j,k)*DZT(imt,j,k))/1.e2
    enddo
  enddo
end subroutine mer_advection


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end program
