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
  input_folder,geometry1_file,geometry2_file,ref_state_file,tavg_file

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
double precision, parameter                     ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,rho0 = 4.1/3.996*1000

imt            = 3600
jmt            = 2400
km             =   42

input_folder   = '/home/dijkbio/andre/LEC/input/'
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

!===============================================================================
!  geometry
!===============================================================================

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
allocate( geometry2(6,km),dz(km),tdepth(km),area(km),p_z(km),vol(km) )
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
call load_3D_field(1,nrec_RHO, RHO ) ! [g/cm^3]
call load_3D_field(1,nrec_PD,  PD  ) ! [g/cm^3]
call load_3D_field(1,nrec_Q,   Q   ) ! [g/cm^4]
call load_3D_field(1,nrec_SALT,SALT) ! [g/kg]
call load_3D_field(1,nrec_TEMP,TEMP) ! [degC]

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
  write(*,*) k,T_avg(k)/T_wavg(k),S_avg(k)/S_wavg(k),RHO_avg(k)/RHO_wavg(k),PD_avg(k)/PD_wavg(k),Q_avg(k)/Q_wavg(k),D0_avg(k),D0_wavg(k),Dp_avg(k),Dp_wavg(k)
!D0_avg(k)/D0_wavg(k),Dp_avg(k)/Dp_wavg(k)
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
end program ref_state
