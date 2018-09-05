!===============================================================================
program OHC_adv
!===============================================================================
implicit none

!===============================================================================
!
!  this script calculates heat advection from the original POP output files
!
!  calculated quantities:
!  1. meridional advection through section at specific latitude
!  2. zonal advection through section at certain longitude
!
!  generated output
!
!  for curvilinear grid:
!  longitude (i) 0=250E, 1100=0E, 1500=40E
!  latitude  (j) 0=78.5S, 518=55S, 677=45S, 866=30S, 1181=0S
!  depth     (k) 1=5m, 9=96.9m, 10=112m, 14=198m, 17=318m, 35=4125, 42=5875m
!
!===============================================================================

!===============================================================================
!  variables
!===============================================================================

character*120 ::                                                               &
  LEC_folder,file_name,file_base,                                              &
  ADV_N_file,ADV_W_file,ADV_E_file,                                            &
  ADV_N_k_file,ADV_W_k_file,ADV_E_k_file,                                      &
  ADV_N_i_file,ADV_W_j_file,ADV_E_j_file,                                      &
  geometry1_file, geometry2_file
character*4   ::                                                               &
  year, yr
character*2   ::                                                               &
  month
integer                             ::                                         &
  imt,jmt,km,y,start_year,end_year,n,rec_length,i,j,k,m,t,                     &
  imin,imax,jmin,jmax,nrec_UVEL,nrec_VVEL,nrec_TEMP
real, dimension(:),     allocatable ::                                         &
  dz,tdepth,area,                                                              &
  ADV_N_k, ADV_W_k, ADV_E_k, ADV_N_i, ADV_W_j, ADV_E_j

real, dimension(:,:),   allocatable ::                                         &
  DXT,DYT,TAREA,DXU,DYU,UAREA,geometry2,                                       &
  ADV_N, ADV_W, ADV_E,                                                         &
  ADV_N_k_year, ADV_W_k_year, ADV_E_k_year, ADV_N_i_year, ADV_W_j_year,        &
  ADV_E_j_year
real, dimension(:,:,:), allocatable ::                                         &
  DZT,DZU,UVEL,VVEL,TEMP,                                                      &
  ADV_N_year, ADV_W_year, ADV_E_year
double precision, parameter         ::                                         &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

imt        = 3600
jmt        = 2400
km         =   42

start_year =   75
end_year   =  326

imin       =  750 ! 35W
imax       = 1900 ! 80E
jmin       =    2 ! ca. 78S
jmax       =  600 ! ca. 60S

nrec_UVEL  =    1
nrec_VVEL  =   43

write (*,*) ''
write (*,*) '--- OHC advection along transects ---'
write (*,*) ''

!===============================================================================
!  FILES
!===============================================================================

LEC_folder      = '/home/dijkbio/andre/LEC/'
geometry1_file  = trim(LEC_folder)//'input/geometry1'
geometry2_file  = trim(LEC_folder)//'input/geometry2'

ADV_N_file   = trim(LEC_folder)//'results/OHC_adv/OHC_adv_N'
ADV_W_file   = trim(LEC_folder)//'results/OHC_adv/OHC_adv_W'
ADV_E_file   = trim(LEC_folder)//'results/OHC_adv/OHC_adv_E'
ADV_N_k_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_N_vint'
ADV_W_k_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_W_vint'
ADV_E_k_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_E_vint'
ADV_N_i_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_N_zint'
ADV_W_j_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_W_mint'
ADV_E_j_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_E_mint'

open( 2,file=ADV_N_file,access='direct',form='unformatted',                    &
       recl=(imax-imin+1)*km,status='unknown')
open( 3,file=ADV_W_file,access='direct',form='unformatted',                    &
       recl=(jmax-jmin+1)*km,status='unknown')
open( 4,file=ADV_E_file,access='direct',form='unformatted',                    &
       recl=(jmax-jmin+1)*km,status='unknown')
open( 5,file=ADV_N_k_file,access='direct',form='unformatted',                  &
       recl=imax-imin+1,status='unknown')
open( 6,file=ADV_W_k_file,access='direct',form='unformatted',                  &
       recl=jmax-jmin+1,status='unknown')
open( 7,file=ADV_E_k_file,access='direct',form='unformatted',                  &
       recl=jmax-jmin+1,status='unknown')
open( 8,file=ADV_N_i_file,access='direct',form='unformatted',                  &
       recl=km,status='unknown')
open( 9,file=ADV_W_j_file,access='direct',form='unformatted',                  &
       recl=km,status='unknown')
open(10,file=ADV_E_j_file,access='direct',form='unformatted',                  &
       recl=km,status='unknown')
!===============================================================================
!  GEOMETRY 
!===============================================================================

! read 2D geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km) )
open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,    &
       status='old')
read(1,rec=1) DXT ! [m]
read(1,rec=2) DXU ! [m]
read(1,rec=3) DYT ! [m]read(1,rec=4) DYU ! [m]
read(1,rec=4) DYU ! [m]
read(1,rec=5) TAREA ! [m^2]
read(1,rec=6) UAREA ! [m^2]
do k=1,km
  read(1,rec=6+   k) DZT(:,:,k) ! [m]
  read(1,rec=6+km+k) DZU(:,:,k) ! [m]
enddo
close(1)

! read 1D geometry fields
allocate( geometry2(6,km),dz(km),tdepth(km),area(km) )
open(1,file=geometry2_file,access='direct',form='unformatted',recl=6,     &
       status='old')
do k=1,km
  read(1,rec=k) geometry2(:,k) ! all real: k,dz[m],tdepth[m],area[m^2],p[bar]
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
close(1)

!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'


allocate( UVEL(imt,jmt,km),      VVEL(imt,jmt,km),     TEMP(imt,jmt,km),       &
          ADV_N(imax-imin+1,km), ADV_N_k(imax-imin+1), ADV_N_i(km),            &
          ADV_W(jmax-jmin+1,km), ADV_W_k(jmax-jmin+1), ADV_W_j(km),            &
          ADV_E(jmax-jmin+1,km), ADV_E_k(jmax-jmin+1), ADV_E_j(km),            &
          ADV_N_year(imax-imin+1,km,12), ADV_N_k_year(imax-imin+1,12), ADV_N_i_year(km,12), &
          ADV_W_year(jmax-jmin+1,km,12), ADV_W_k_year(jmax-jmin+1,12), ADV_W_j_year(km,12), &
          ADV_E_year(jmax-jmin+1,km,12), ADV_E_k_year(jmax-jmin+1,12), ADV_E_j_year(km,12) )

! loop over years
do y=start_year,end_year

  if ( y<276 ) then
    file_base = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc/tavg/t.t0.1_42l_nccs01.'
    nrec_TEMP  =  127 ! for 12/14 GB files
  elseif ( y>=276 ) then
    file_base = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_viebahn/tavg/t.t0.1_42l_nccs01.'
    nrec_TEMP  =  222 ! in ev case
  endif

  ! loop over months
  do m=1,12
    if ( m<12 ) then
      write(year,"(I0.4)") y
      write(month,"(I0.2)") m+1
    elseif ( m==12 ) then
      write(year,"(I0.4)") y+1
      write(month,"(I0.2)") 1
    endif

    ! file 27601 and 28101 missing!
    if ( y==275 .and. m==12 ) then
      write(year,"(I0.4)") y
      write(month,"(I0.2)") m
    elseif ( y==280 .and. m==12 ) then
      write(year,"(I0.4)") y
      write(month,"(I0.2)") m
    endif

    file_name = trim(file_base)//year//month
    open(1,file=file_name,access='direct',form='unformatted',                  &
           recl=imt*jmt,status='unknown')

    write(*,*) y, file_name

    do k=1,km
      read(1,rec=nrec_UVEL+k-1) UVEL(:,:,k)
      read(1,rec=nrec_VVEL+k-1) VVEL(:,:,k)
      read(1,rec=nrec_TEMP+k-1) TEMP(:,:,k)
    enddo

    !===========================================================================
    ! OHC advection
    !===========================================================================


    ! N
    call OHC_merid_adv(imin,imax,jmax,TEMP,VVEL,DZT,DXT,ADV_N,ADV_N_k,ADV_N_i)
    ! W
    call OHC_zonal_adv(jmin,jmax,imin,TEMP,UVEL,DZT,DYT,ADV_W,ADV_W_k,ADV_W_j)
    ! E
    call OHC_zonal_adv(jmin,jmax,imax,TEMP,UVEL,DZT,DYT,ADV_E,ADV_E_k,ADV_E_j)


    !===========================================================================
    ! temporary storage
    ADV_N_year(:,:,m) = ADV_N(:,:)
    ADV_W_year(:,:,m) = ADV_W(:,:)
    ADV_E_year(:,:,m) = ADV_E(:,:)

    ADV_N_k_year(:,m) = ADV_N_k(:)
    ADV_W_k_year(:,m) = ADV_W_k(:)
    ADV_E_k_year(:,m) = ADV_E_k(:)

    ADV_N_i_year(:,m) = ADV_N_i(:)
    ADV_W_j_year(:,m) = ADV_W_j(:)
    ADV_E_j_year(:,m) = ADV_E_j(:)

    ! yearly averages
    if ( m==12 ) then
      ADV_N(:,:) = 0.0
      ADV_W(:,:) = 0.0
      ADV_E(:,:) = 0.0
      ADV_N_k(:) = 0.0
      ADV_W_k(:) = 0.0
      ADV_E_k(:) = 0.0
      ADV_N_i(:) = 0.0
      ADV_W_j(:) = 0.0
      ADV_E_j(:) = 0.0
      do t=1,12
        ADV_N(:,:) = ADV_N(:,:) + ADV_N_year(:,:,t)/12.0
        ADV_W(:,:) = ADV_W(:,:) + ADV_W_year(:,:,t)/12.0
        ADV_E(:,:) = ADV_E(:,:) + ADV_E_year(:,:,t)/12.0
        ADV_N_k(:) = ADV_N_k(:) + ADV_N_k_year(:,t)/12.0
        ADV_W_k(:) = ADV_W_k(:) + ADV_W_k_year(:,t)/12.0
        ADV_E_k(:) = ADV_E_k(:) + ADV_E_k_year(:,t)/12.0
        ADV_N_i(:) = ADV_N_i(:) + ADV_N_i_year(:,t)/12.0
        ADV_W_j(:) = ADV_W_j(:) + ADV_W_j_year(:,t)/12.0
        ADV_E_j(:) = ADV_E_j(:) + ADV_E_j_year(:,t)/12.0
      enddo !t
      !=========================================================================
      !  OUTPUT
      !=========================================================================
      write( 2,rec=y-start_year+1) ADV_N(:,:)
      write( 3,rec=y-start_year+1) ADV_W(:,:)
      write( 4,rec=y-start_year+1) ADV_E(:,:)
      write( 5,rec=y-start_year+1) ADV_N_k(:)
      write( 6,rec=y-start_year+1) ADV_W_k(:)
      write( 7,rec=y-start_year+1) ADV_E_k(:)
      write( 8,rec=y-start_year+1) ADV_N_i(:)
      write( 9,rec=y-start_year+1) ADV_W_j(:)
      write(10,rec=y-start_year+1) ADV_E_j(:)
    endif

    close (1)
  enddo !m
enddo !y



contains

subroutine OHC_merid_adv(imin,imax,lat,TEMP,VVEL,DZT,DXT,ADV,ADV_k,ADV_i)
! merdional advection at latitude lat
! northward flux through Northern face of T-cell(i,lat)
! ADV       in [J/m^2/s]
! ADV_k/i   in [J/m/s]

integer                                         :: i,k
integer,                            intent(in)  :: imin,imax,lat
integer, parameter                              :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt),        intent(in)  :: DXT
real,    dimension(imt,jmt,km),     intent(in)  :: TEMP, VVEL, DZT
real,    dimension(imax-imin+1,km), intent(out) :: ADV
real,    dimension(km),             intent(out) :: ADV_i ! zonal-int
real,    dimension(imax-imin+1),    intent(out) :: ADV_k ! depth-int
real,    parameter                              :: rho0=4.1/3.996*1000, c=3996.0

do k=1,km
  do i=imin,imax
    ADV(i-imin+1,k) = (TEMP(i,lat+1,k)-TEMP(i,lat,k))                          &
                    * (VVEL(i-1,lat,k)+VVEL(i,lat,k))/2                        &  
                    * rho0 * c
  enddo !i
enddo !k

ADV_k(:) = 0.0
do k=1,km
  ADV_k(:) = ADV_k(:) +  ADV(:,k)*DZT(imin:imax,lat,k)
enddo !k

ADV_i(:) = 0.0
do i=imin,imax
  ADV_i(:) = ADV_i(:) + ADV(i-imin+1,:)*DXT(i,lat)
enddo !i

end subroutine OHC_merid_adv

subroutine OHC_zonal_adv(jmin,jmax,lon,TEMP,UVEL,DZT,DYT,ADV,ADV_k,ADV_j)
! zonal advection at longitude imin/imax
! eastward flux through E face of T-cell(lon,i)
! ADV       in [J/m^2/s]
! ADV_k/j   in [J/m/s]

integer                                         :: j,k
integer,                            intent(in)  :: jmin,jmax,lon
integer, parameter                              :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt),        intent(in)  :: DYT
real,    dimension(imt,jmt,km),     intent(in)  :: TEMP, UVEL, DZT
real,    dimension(jmax-jmin+1,km), intent(out) :: ADV
real,    dimension(km),             intent(out) :: ADV_j ! merid-int
real,    dimension(jmax-jmin+1),    intent(out) :: ADV_k ! depth-int
real,    parameter                              :: rho0=4.1/3.996*1000, c=3996.0

do k=1,km
  do j=jmin,jmax
    ADV(j-jmin+1,k) = (TEMP(lon+1,j,k)-TEMP(lon,j,k))                          &
                    * (UVEL(lon,j-1,k)+UVEL(lon,j,k))/2                        &  
                    * rho0 * c
  enddo !j
enddo !k

ADV_k(:) = 0.0
do k=1,km
  ADV_k(:) = ADV_k(:) + ADV(:,k)*DZT(lon,jmin:jmax,k)
enddo !k

ADV_j(:) = 0.0
do j=jmin,jmax
  ADV_j(:) = ADV_j(:) + ADV(j-jmin+1,:)*DYT(lon,j)
enddo !i

end subroutine OHC_zonal_adv
!===============================================================================
end program OHC_adv
!===============================================================================
