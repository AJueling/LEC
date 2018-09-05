!===============================================================================
program MOC_34S
!===============================================================================
implicit none

character*120 ::                                                               &
  grid_file,kmt_file,in_depths,pbc_file,input_folder,file_name,file_base,      &
  geometry1_file, geometry2_file
character*4   ::                                                               &
  year, yr
character*2   ::                                                               &
  month
integer                             ::                                         &
  imt,jmt,km,y,start_year,end_year,nt,nrec_UVEL,nrec_VVEL,n,rec_length,i,j,k,m,t
real, dimension(:), allocatable ::                                             &
  AMOC, GMOC, AMOC_avg, GMOC_avg, dz, tdepth, area
real, dimension(:,:), allocatable ::                                           &
  DXT, DYT, TAREA, DXU, DYU, UAREA, AMOC_year, GMOC_year, geometry2
real, dimension(:,:,:), allocatable ::                                         &
  DZT, DZU, VVEL
double precision, parameter         ::                                         &
  c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

imt        = 3600
jmt        = 2400
km         =   42

start_year =   75
end_year   =  326
nt         = end_year-start_year+1
n          =   51

nrec_UVEL  =    1
nrec_VVEL  =   43

input_folder    = '/home/dijkbio/andre/LEC/input/'
geometry1_file  = trim(input_folder)//'geometry1'
geometry2_file  = trim(input_folder)//'geometry2'

open(2,file='/home/dijkbio/andre/LEC/results/AMOC/AMOC_34.5S.bin',             &
       access='direct',form='unformatted',recl=km,status='unknown')
open(3,file='/home/dijkbio/andre/LEC/results/AMOC/GMOC_34.5S.bin',             &
       access='direct',form='unformatted',recl=km,status='unknown')

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

allocate( VVEL(imt,jmt,km), AMOC(km), GMOC(km),                                &
          AMOC_year(km,12), GMOC_year(km,12), AMOC_avg(km), GMOC_avg(km) )

! loop over years
do y=start_year,end_year

  if ( y<276 ) then
    file_base = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc/tavg/t.t0.1_42l_nccs01.'
  elseif ( y>=276 ) then
    file_base = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_viebahn/tavg/t.t0.1_42l_nccs01.'
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
      read(1,rec=nrec_VVEL+k-1) VVEL(:,:,k)
    enddo

    !===========================================================================
    ! OSF
    !===========================================================================
    ! preparing integral dz/dy + executing integral dx
    do k = 1,km
      AMOC(k) = -sum(VVEL(401:1400,815,k)*DZU(401:1400,815,k)*DXU(401:1400,815)/1.0E08,1)
      GMOC(k) = -sum(VVEL(:,815,k)*DZU(:,815,k)*DXU(:,815)/1.0E08,1)
    enddo

    ! OSF_v (depth integral dz bottom to depth z)
    do k = 1,km-1
      AMOC(km-k) = AMOC(km-k) + AMOC(km-k+1)
      GMOC(km-k) = GMOC(km-k) + GMOC(km-k+1)
    enddo

    AMOC_year(:,m) = AMOC(:)
    GMOC_year(:,m) = GMOC(:)

    ! yearly average OSF
    if ( m==12 ) then
      AMOC_avg(:) = 0.0
      GMOC_avg(:) = 0.0
      do t=1,12
        AMOC_avg(:) = AMOC_avg(:) + AMOC_year(:,t)/12.0
        GMOC_avg(:) = GMOC_avg(:) + GMOC_year(:,t)/12.0
      enddo !t
      ! write to output file
      write(2,rec=y-start_year+1) AMOC_avg(:)
      write(3,rec=y-start_year+1) GMOC_avg(:)
    endif

    close (1)
  enddo !m
enddo !y

!===============================================================================
end program MOC_34S
!===============================================================================
