!===============================================================================
program OSF_BSF
!===============================================================================
implicit none

!===============================================================================
!
!  this script calculates various quantities from the original POP output files
!  and the LEC files created with LEC.f90
!
!  calculated quantities:
!  2. Eulerian meridional overturning stream function
!
!  generated output
!  2. two binary file with Eul. stream functions as an output, one calculated
!     from v and one from w
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
  grid_file,kmt_file,in_depths,pbc_file,input_folder,file_name,file_base,      &
  OSF_v_file,OSF_v_mean_file,OSF_v_std_file,OSF_folder,                        &
  BSF_file,BSF_mean_file,BSF_std_file,BSF_folder,BSF_values_file,              &
  geometry1_file, geometry2_file
character*4   ::                                                               &
  year, yr
character*2   ::                                                               &
  month
integer                             ::                                         &
  imt,jmt,km,y,start_year,end_year,nt,nrec_UVEL,nrec_VVEL,n,rec_length,i,j,k,m,t
real                                ::                                         &
  BSF_drake, BSF_gyre
real, dimension(:),     allocatable ::                                         &
  dz,tdepth,area
real, dimension(:,:),   allocatable ::                                         &
  DXT,DYT,TAREA,DXU,DYU,UAREA,geometry2,                                       &
  OSF_v,OSF_v_avg,OSF_v_sum,OSF_v_std,OSF_v_mean,                              &
  BSF,BSF_avg,BSF_sum,BSF_std,BSF_mean
real, dimension(:,:,:), allocatable ::                                         &
  DZT,DZU,UVEL,VVEL,OSF_v_year,BSF_year
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

write (*,*) ''
write (*,*) '--- Overturning and Barotropic Stream Functions ---'
write (*,*) ''

!===============================================================================
!  FILES
!===============================================================================

input_folder    = '/home/dijkbio/andre/LEC/input/'
geometry1_file  = trim(input_folder)//'geometry1'
geometry2_file  = trim(input_folder)//'geometry2'

OSF_folder      = '/projects/0/samoc/jan/Andree/OSF/'
BSF_folder      = '/projects/0/samoc/jan/Andree/BSF/'

OSF_v_file      = trim(OSF_folder)//'OSF_v_yrly_avg'
OSF_v_mean_file = trim(OSF_folder)//'OSF_v_mean'
OSF_v_std_file  = trim(OSF_folder)//'OSF_v_std'
BSF_file        = trim(BSF_folder)//'BSF_yrly_avg'
BSF_mean_file   = trim(BSF_folder)//'BSF_mean'
BSF_std_file    = trim(BSF_folder)//'BSF_std'
BSF_values_file = trim(BSF_folder)//'BSF_drake_gyre'

open(2,file=OSF_v_file,access='direct',form='unformatted',                     &
       recl=jmt*km,status='unknown')
open(3,file=OSF_v_mean_file,access='direct',form='unformatted',                &
       recl=jmt*km,status='unknown')
open(4,file=OSF_v_std_file,access='direct',form='unformatted',                 &
       recl=jmt*km,status='unknown')
open(5,file=BSF_file,access='direct',form='unformatted',                       &
       recl=imt*jmt,status='unknown')
open(6,file=BSF_mean_file,access='direct',form='unformatted',                  &
       recl=imt*jmt,status='unknown')
open(7,file=BSF_std_file,access='direct',form='unformatted',                   &
       recl=imt*jmt,status='unknown')
open(8,file=BSF_values_file,access='direct',form='unformatted',                &
       recl=3,status='unknown')
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

write (*,*) ''
write (*,*) 'CALCULATIONS'


allocate( UVEL(imt,jmt,km), VVEL(imt,jmt,km),                                  &
          OSF_v(jmt,km),      OSF_v_avg(jmt,km),                               &
          OSF_v_mean(jmt,km), OSF_v_sum(jmt,km), OSF_v_year(jmt,km,12),        &
          BSF(imt,jmt),       BSF_avg(imt,jmt),                                &
          BSF_mean(imt,jmt),  BSF_sum(imt,jmt),   BSF_year(imt,jmt,12) )

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
      read(1,rec=nrec_UVEL+k-1) UVEL(:,:,k)
      read(1,rec=nrec_VVEL+k-1) VVEL(:,:,k)
    enddo

    !===========================================================================
    ! OSF
    !===========================================================================
    ! preparing integral dz/dy + executing integral dx
    do k = 1,km
      OSF_v(:,k) = -sum(VVEL(:,:,k)*DZU(:,:,k)*DXU(:,:)/1.0E06,1)
    enddo

    ! OSF_v (depth integral dz bottom to depth z)
    do k = 1,km-1
      OSF_v(:,km-k) = OSF_v(:,km-k) + OSF_v(:,km-k+1)
    enddo

    !===========================================================================
    ! BSF
    !===========================================================================

    ! depth integral
    BSF(:,:) = sum(UVEL(:,:,:)*DZU(:,:,:),3)*DYU(:,:)/1.0E02 ! [m^2/s]

    ! meridional integral
    do j = 2,jmt
      BSF(:,j) = BSF(:,j) + BSF(:,j-1)
    enddo


    !===========================================================================
    ! temporary storage
    OSF_v_year(:,:,m) = OSF_v(:,:)
    BSF_year(:,:,m)   = BSF(:,:)

    ! yearly average OSF
    if ( m==12 ) then
      OSF_v_avg(:,:) = 0.0
      BSF_avg(:,:)   = 0.0
      do t=1,12
        OSF_v_avg(:,:) = OSF_v_avg(:,:) + OSF_v_year(:,:,t)/12.0
        BSF_avg(:,:)   = BSF_avg(:,:)   + BSF_year(:,:,t)/12.0
      enddo !t
      ! write to output file
      write(2,rec=y-start_year+1) OSF_v_avg(:,:)
      write(5,rec=y-start_year+1) BSF_avg(:,:)
      BSF_drake = MAXVAL(BSF_avg(446,311:521))
      BSF_gyre  = MINVAL(BSF_avg(600:1300,100:500))
      write(8,rec=y-start_year+1) y, BSF_drake, BSF_gyre
    endif

    close (1)
  enddo !m
enddo !y


!===============================================================================
! mean of last n years
!===============================================================================

OSF_v_mean(:,:) = 0.0
BSF_mean(:,:)   = 0.0
do y=1,n
  read(2,rec=end_year-start_year-y+1) OSF_v
  OSF_v_mean(:,:) = OSF_v_mean(:,:) + 1.0/float(n)*OSF_v(:,:)
  read(5,rec=end_year-start_year-y+1) BSF
  BSF_mean(:,:)   = BSF_mean(:,:) + 1.0/float(n)*BSF(:,:)
enddo

write(3,rec=1) OSF_v_mean
write(6,rec=1) OSF_v_mean

!===============================================================================
! standard deviation of last n years
!===============================================================================

OSF_v_sum(:,:) = 0.0
BSF_sum(:,:)   = 0.0
do y=1,n
  read(2,rec=end_year-start_year-y+1) OSF_v
  OSF_v_sum = OSF_v_sum + ( OSF_v - OSF_v_mean )**2
  read(5,rec=end_year-start_year-y+1) BSF
  BSF_sum   = BSF_sum + ( BSF - BSF_mean )**2
enddo

write(4,rec=1) sqrt(OSF_v_sum/(n-1))
write(7,rec=1) sqrt(BSF_sum/(n-1))

!===============================================================================
end program OSF_BSF
!===============================================================================
