program Psi_analysis
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
!  latitude  (j) 0=78.5S, 518=55S, 677=45S, 866=30S, 1181=0S
!  depth     (k) 1=5m, 9=96.9m, 10=112m, 14=198m, 17=318m, 35=4125, 42=5875m
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!===============================================================================
!  variables
!===============================================================================

character*120 ::                                                               &
  grid_file,kmt_file,in_depths,pbc_file,input_folder,output_folder,file_base,  &
  Psi_v_file,Psi_v_mean_file,Psi_v_std_file,file_name
character*4   ::                                                               &
  year, yr
character*2   ::                                                               &
  month
integer                                         ::                             &
  imt,jmt,km,y,start_year,end_year,nt,nrec_VVEL,n,rec_length,i,j,k,m,t
integer,          dimension(:,:),   allocatable ::                             &
  kmT
double precision                                ::                             &
  volume 
double precision, dimension(:),     allocatable ::                             &
  dz,area
double precision, dimension(:,:),   allocatable ::                             &
  HTN,HTE,DXT,DYT,TAREA,DXU,DYU,UAREA,DZBC,HUW,HUS,WORK,WORK2,WORK3  
real,             dimension(:),     allocatable ::                             &
  tdepth
real,             dimension(:,:),   allocatable ::                             &
  Psi_v,Psi_v_avg,Psi_v_sum,Psi_v_std,Psi_v_mean
real,             dimension(:,:,:), allocatable ::                             &
  DZT,DZU,VVEL,Psi_v_year
double precision, parameter                     ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

imt        = 3600
jmt        = 2400
km         =   42

start_year =   75
end_year   =  326
nt         = end_year-start_year+1
n          =   51

nrec_VVEL  =   43

write (*,*) ''
write (*,*) '--- Psi ---'
write (*,*) ''

!===============================================================================
!  FILES
!===============================================================================
input_folder    = '/home/dijkbio/andre/LEC/input/'

grid_file       = trim(input_folder)//'grid.3600x2400.fob.da'
kmt_file        = trim(input_folder)//'kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'
in_depths       = trim(input_folder)//'in_depths.42.dat'
pbc_file        = trim(input_folder)//'dzbc_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'

output_folder   = '/projects/0/samoc/jan/Andree/Psi/'

Psi_v_file      = trim(output_folder)//'Psi_v_yrly_avg'
Psi_v_mean_file = trim(output_folder)//'Psi_v_mean'
Psi_v_std_file  = trim(output_folder)//'Psi_v_std'

open(2,file=Psi_v_file,access='direct',form='unformatted',                     &
       recl=jmt*km,status='unknown')
open(3,file=Psi_v_mean_file,access='direct',form='unformatted',                &
       recl=jmt*km,status='unknown')
open(4,file=Psi_v_std_file,access='direct',form='unformatted',                 &
       recl=jmt*km,status='unknown')
!===============================================================================
!  GEOMETRY 
!===============================================================================

allocate(                                                                      &
HTN(imt,jmt), HTE(imt,jmt), HUW(imt,jmt), HUS(imt,jmt), WORK(imt,jmt),         &
WORK2(imt,jmt), WORK3(imt,jmt), DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt),    &
DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), dz(km), tdepth(km), kmT(imt,jmt),  &
DZBC(imt,jmt), DZT(imt,jmt,km), DZU(imt,jmt,km) )

inquire (iolength = rec_length) HTN
open    (1,file=grid_file,access='direct',form='unformatted',recl=rec_length,  &
           status='unknown')
  read  (1,rec=3) HTN ! [cm]
  read  (1,rec=4) HTE ! [cm]
close   (1)

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

! kmT
inquire (iolength=rec_length) kmT
open    (1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,   &
           status='unknown')
  read  (1,rec=1) kmT
close   (1)

! dz
open    (1,file=in_depths,status='old')
  do k = 1, km
    read  (1,*) dz(k) ! [cm]
  enddo !k
close   (1)

! partial bottom cell depths
inquire (iolength=rec_length) DZBC
open    (1,file=pbc_file,access='direct',form='unformatted',recl=rec_length,   &
           status='unknown')
  read  (1,rec=1) DZBC ! [cm]
close   (1)
!write   (*,*) 'PBC file: ', pbc_file

! create DZT-field, also used as 3D TTT-mask
do k=1,km
  where(kmT>k)
    DZT(:,:,k) = dz(k)
  elsewhere(kmT==k)
    DZT(:,:,k) = DZBC
  elsewhere
    DZT(:,:,k) = c0
  endwhere
enddo !k
  
! DZU field
do k=1,km
  do j=1,jmt-1   ! eastern edge
      DZU(imt,j,k) = min(DZT(imt,j,k),DZT(1,j,k),DZT(imt,j+1,k),DZT(1,j+1,k))
    do i=1,imt-1  ! rest of array
      if ( k.le.kmT(i,j) ) then
        DZU(i,j,k) = min(DZT(i,j,k),DZT(i+1,j,k),DZT(i,j+1,k),DZT(i+1,j+1,k))
      endif
    enddo !j
  enddo !j
  ! northern edge
  do i=1,imt-1
    if ( i.lt.imt/2 ) then      ! first half of array, before wrapping around
      DZU(i,j,k) = min( DZT(i    ,jmt,k), DZT(i+1    ,jmt,k),                  &
                        DZT(imt-i,jmt,k), DZT(imt-i+1,jmt,k) )   
    elseif ( i.gt.imt/2 ) then  ! second half, after wrap
      DZU(i,jmt,k) = DZU(imt-i,jmt,k)
    endif
  enddo !i
enddo !k

!  create tdepth 1D array containing depths of TTT-points in [cm]
tdepth(1) = dz(1)/2.0
do k = 2,km
  tdepth(k) = tdepth(k-1)+p5*(dz(k-1)+dz(k))
enddo !k

deallocate( WORK, WORK2, WORK3, DZBC )

allocate( area(km) )

volume  = sum(TAREA*sum(DZT,3))*1.0E-06 ! volume between bottom and z=0 [m^3]
!write (*,"(A10, ES13.4E2, A4)") 'volume  =',  volume, 'm^3'

do k = 1, km
  area(k) = sum(TAREA,DZT(:,:,k).ne.0)*1.0E-04 ! [m^2]
!  write (*,*) area(k)
enddo !k

!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'


allocate( VVEL(imt,jmt,km), Psi_v(jmt,km), Psi_v_avg(jmt,km),                  &
          Psi_v_mean(jmt,km), Psi_v_sum(jmt,km), Psi_v_year(jmt,km,12) )

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

    ! preparing integral dz/dy + executing integral dx
    do k = 1,km
      Psi_v(:,k) = -sum(VVEL(:,:,k)*DZU(:,:,k)*DXU(:,:)/1.0E06,1)
    enddo

    ! Psi_v (depth integral dz bottom to depth z)
    do k = 1,km-1
      Psi_v(:,km-k) = Psi_v(:,km-k) + Psi_v(:,km-k+1)
    enddo

    ! temporary storage
    Psi_v_year(:,:,m) = Psi_v(:,:)

    ! yearly average Psi
    if ( m==12 ) then
      Psi_v_avg(:,:) = 0.0
      do t=1,12
        Psi_v_avg(:,:) = Psi_v_avg(:,:) + Psi_v_year(:,:,t)/12.0
      enddo !t
      ! write to output file
      write(2,rec=y-start_year+1) Psi_v_avg(:,:)
    endif


    close (1)
  enddo !m
enddo !y


!===============================================================================
! mean of last n years
!===============================================================================

Psi_v_mean(:,:) = 0.0
do y=1,n
  read(2,rec=end_year-start_year-y+1) Psi_v
  Psi_v_mean(:,:) = Psi_v_mean(:,:) + 1.0/float(n)*Psi_v(:,:)
enddo

write(3,rec=1) Psi_v_mean

!===============================================================================
! standard deviation of last n years
!===============================================================================

Psi_v_sum(:,:) = 0.0

do y=1,n
  read(2,rec=end_year-start_year-y+1) Psi_v
  Psi_v_sum = Psi_v_sum + ( Psi_v - Psi_v_mean )**2
enddo

write(4,rec=1) sqrt(Psi_v_sum/(n-1))












!===============================================================================
!  SUBROUTINES
!===============================================================================

contains
end program
