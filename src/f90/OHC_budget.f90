!===============================================================================
program OHC_adv
!===============================================================================
implicit none

!===============================================================================
!
!  This script calculates OHC budget from the original POP output files.
!  The area is a box around the  Weddel Sea deep water formation spot:
!  [78,60S]x[35W,80E] 
!
!  components:
!  (1) advection:
!      a. meridional advection through section at specific latitude
!      b. zonal advection through section at certain longitude
!  (2) d(OHC)/dt
!  (3) surface heat fluxes
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
  geometry1_file, geometry2_file
character*4   ::                                                               &
  year, yr
character*2   ::                                                               &
  month
integer                             ::                                         &
  imt,jmt,km,y,start_year,end_year,n,rec_length,i,j,k,m,t,                     &
  imin,imax,jmin,jmax,                                                         &
  nrec_UVEL,nrec_VVEL,nrec_TEMP
real, dimension(:),     allocatable ::                                         &
  dz,tdepth,area!,                                                              &

real, dimension(:,:),   allocatable ::                                         &
  DXT,DYT,TAREA,DXU,DYU,UAREA,geometry2!,                                       &

real, dimension(:,:,:), allocatable ::                                         &
  DZT,DZU,UVEL,VVEL,TEMP!,                                                      &

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

ADV_E_j_file = trim(LEC_folder)//'results/OHC_adv/OHC_adv_E_mint'

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


    !===========================================================================
    ! d(OHC)/dt
    !===========================================================================


    !===========================================================================
    ! surface heat flux
    !===========================================================================




    !===========================================================================
    ! temporary storage
    ADV_N_year(:,:,m) = ADV_N(:,:)

    ! yearly averages
    if ( m==12 ) then
      ADV_E_j(:) = 0.0
      do t=1,12
        ADV_E_j(:) = ADV_E_j(:) + ADV_E_j_year(:,t)/12.0
      enddo !t
      !=========================================================================
      !  OUTPUT
      !=========================================================================
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
