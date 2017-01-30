!
program MXL
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! 1. creates yearly averages of mixed layer fields
! 2. calculates mean and standard deviation for last n(=51) years
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
implicit none

!===============================================================================
!  variables
!===============================================================================

integer                              ::                                        &
  imt,jmt,start_year,end_year,y,m,nt,n,nrec_HMXL,nrec_XMXL,nrec_TMXL,nrec_TEMP,&
  counter,km,t,rec_length,k
integer, dimension(:,:), allocatable ::                                        &
  XMASK
character*160                        ::                                        & 
  cmd,file_base,file_name,output_folder,grid_file,                             &
  HMXL_file,XMXL_file,TMXL_file,HMXL_mean_file,XMXL_mean_file,TMXL_mean_file,  &
  HMXL_std_file,XMXL_std_file,TMXL_std_file,XMASK_file,                        &
  HMXL_ymask_file,XMXL_ymask_file,TMXL_ymask_file,TEMP_ymask_file,             &
  HMXL_mmask_file,XMXL_mmask_file,TMXL_mmask_file,TEMP_mmask_file,             &
  XMXL_ymask_max_file,XMXL_mmask_max_file
character*4                          ::                                        &
  year
character*2                          ::                                        &
  month
real                                 ::                                        &
  file_size,XMASK_area
real,    dimension(:),   allocatable ::                                        &
  TEMP_avg, TEMP_avg2
real,    dimension(:,:), allocatable ::                                        &
  TEMP,HMXL,XMXL,TMXL,HMXL_avg,XMXL_max,TMXL_min,                              &
  HMXL_mean,XMXL_mean,TMXL_mean,HMXL_sum,XMXL_sum,TMXL_sum
double precision, dimension(:,:),   allocatable ::                             &
  HTN,HTE,DXT,DYT,TAREA,WORK
double precision, parameter                     ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

imt        = 3600
jmt        = 2400
km         =   42

start_year =   75
end_year   =  326
nt         = end_year-start_year+1
n          =   51

!===============================================================================
grid_file = '/home/dijkbio/andre/LEC/input/grid.3600x2400.fob.da'
allocate( HTN(imt,jmt), HTE(imt,jmt), WORK(imt,jmt),                           &
          DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt) )
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
!===============================================================================

output_folder   = '/projects/0/samoc/jan/Andree/MXL/'

HMXL_file       = trim(output_folder)//'HMXL_yrly_avg'
XMXL_file       = trim(output_folder)//'XMXL_yrly_avg'
TMXL_file       = trim(output_folder)//'TMXL_yrly_avg'

HMXL_mean_file  = trim(output_folder)//'HMXL_mean'
XMXL_mean_file  = trim(output_folder)//'XMXL_mean'
TMXL_mean_file  = trim(output_folder)//'TMXL_mean'

HMXL_std_file   = trim(output_folder)//'HMXL_std'
XMXL_std_file   = trim(output_folder)//'XMXL_std'
TMXL_std_file   = trim(output_folder)//'TMXL_std'

XMASK_file      = trim(output_folder)//'XMASK'

HMXL_ymask_file = trim(output_folder)//'HMXL_XMASK_y'
XMXL_ymask_file = trim(output_folder)//'XMXL_XMASK_y'
TMXL_ymask_file = trim(output_folder)//'TMXL_XMASK_y'
TEMP_ymask_file = trim(output_folder)//'TEMP_XMASK_y'
XMXL_ymask_max_file = trim(output_folder)//'XMXL_XMASK_max_y'


HMXL_mmask_file = trim(output_folder)//'HMXL_XMASK_m'
XMXL_mmask_file = trim(output_folder)//'XMXL_XMASK_m'
TMXL_mmask_file = trim(output_folder)//'TMXL_XMASK_m'
TEMP_mmask_file = trim(output_folder)//'TEMP_XMASK_m'
XMXL_mmask_max_file = trim(output_folder)//'XMXL_XMASK_max_m'

open(2,file=HMXL_file,access='direct',form='unformatted',                      &
       recl=imt*jmt,status='unknown')
open(3,file=XMXL_file,access='direct',form='unformatted',                      &
       recl=imt*jmt,status='unknown')
open(4,file=TMXL_file,access='direct',form='unformatted',                      &
       recl=imt*jmt,status='unknown')

allocate( HMXL(imt,jmt),      XMXL(imt,jmt),      TMXL(imt,jmt),               &
          HMXL_avg(imt,jmt),  XMXL_max(imt,jmt),  TMXL_min(imt,jmt),           &
          HMXL_mean(imt,jmt), XMXL_mean(imt,jmt), TMXL_mean(imt,jmt),          &
          HMXL_sum(imt,jmt),  XMXL_sum(imt,jmt),  TMXL_sum(imt,jmt) )

!===============================================================================
!  calculating average, max and min per year
!===============================================================================

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

    ! reading the file size, because different binary files sizes with different
    ! fields were written out: up to 275: 12GB and 14GB mixed because some files
    ! were lost and had to be recreated with a newer output script, 
    ! and after 276: 56 GB with extra variables
    cmd = "stat --format='%s.0' "//trim(file_name)//" > size"
    call system(cmd) 
    open(11,file='size')
111 format (F14.1)
    read(11,*) file_size
    cmd = 'rm size'
    call system(cmd)
    close(11)

    write(*,*) y, file_name, file_size
    if     ( file_size==11923200000.0 ) then
      nrec_HMXL = 343
      nrec_XMXL = 344
      nrec_TMXL = 345
    elseif ( file_size==14826240000.0 ) then 
      nrec_HMXL = 427
      nrec_XMXL = 428
      nrec_TMXL = 429
    elseif ( file_size==59132160000.0 ) then 
      nrec_HMXL = 1709
      nrec_XMXL = 1710
      nrec_TMXL = 1711
    endif
    read(1,rec=nrec_HMXL) HMXL
    read(1,rec=nrec_XMXL) XMXL
    read(1,rec=nrec_TMXL) TMXL

    ! calculating average, maximum and minimum
    if ( m==1 ) then
      HMXL_avg(:,:) = 0.0
      XMXL_max(:,:) = XMXL(:,:) 
      TMXL_min(:,:) = TMXL(:,:) 
    endif

    HMXL_avg(:,:) = HMXL_avg(:,:) + 1.0/12.0*HMXL(:,:)
    where ( XMXL-XMXL_max>0.0 ) XMXL_max = XMXL
    where ( TMXL-TMXL_min<0.0 ) TMXL_min = TMXL

    close (1)
  enddo !m

  write(2,rec=y-start_year+1) HMXL_avg(:,:)
  write(3,rec=y-start_year+1) XMXL_max(:,:)
  write(4,rec=y-start_year+1) TMXL_min(:,:)
enddo !y

!===============================================================================
! mean of last n years
!===============================================================================
open( 5,file=HMXL_mean_file,access='direct',form='unformatted',                &
        recl=imt*jmt,status='unknown')
open( 6,file=XMXL_mean_file,access='direct',form='unformatted',                &
        recl=imt*jmt,status='unknown')
open( 7,file=TMXL_mean_file,access='direct',form='unformatted',                &
        recl=imt*jmt,status='unknown')

HMXL_mean(:,:) = 0.0
XMXL_mean(:,:) = 0.0
TMXL_mean(:,:) = 0.0
do y=1,n
  read(2,rec=end_year-start_year-y+2) HMXL
  read(3,rec=end_year-start_year-y+2) XMXL
  read(4,rec=end_year-start_year-y+2) TMXL
  HMXL_mean(:,:) = HMXL_mean(:,:) + 1.0/float(n)*HMXL(:,:)
  XMXL_mean(:,:) = XMXL_mean(:,:) + 1.0/float(n)*XMXL(:,:)
  TMXL_mean(:,:) = TMXL_mean(:,:) + 1.0/float(n)*TMXL(:,:)
enddo

write(5,rec=1) HMXL_mean
write(6,rec=1) XMXL_mean
write(7,rec=1) TMXL_mean

!===============================================================================
! standard deviation of last n years
!===============================================================================
write(*,*) 'test'
open( 8,file=HMXL_std_file,access='direct',form='unformatted',                 &
        recl=imt*jmt,status='unknown')
write(*,*) 'test'
open( 9,file=XMXL_std_file,access='direct',form='unformatted',                 &
        recl=imt*jmt,status='unknown')
open(10,file=TMXL_std_file,access='direct',form='unformatted',                 &
        recl=imt*jmt,status='unknown')

HMXL_sum(:,:) = 0.0
XMXL_sum(:,:) = 0.0
TMXL_sum(:,:) = 0.0

do y=1,n
  read(2,rec=end_year-start_year-y+2) HMXL
  read(3,rec=end_year-start_year-y+2) XMXL
  read(4,rec=end_year-start_year-y+2) TMXL
  HMXL_sum = HMXL_sum + ( HMXL - HMXL_mean )**2
  XMXL_sum = XMXL_sum + ( XMXL - XMXL_mean )**2
  TMXL_sum = TMXL_sum + ( TMXL - TMXL_mean )**2
enddo

write( 8,rec=1) sqrt(HMXL_sum/(n-1))
write(*,*) 'test'
write( 9,rec=1) sqrt(XMXL_sum/(n-1))
write(10,rec=1) sqrt(TMXL_sum/(n-1))

!close( 2) ! yrly_avg
!close( 3)
!close( 4)
!close( 5) ! mean
!close( 6)
!close( 7)
close( 8) ! std
close( 9)
close(10)

!===============================================================================
! create mask
!===============================================================================

allocate( XMASK(imt,jmt) )

XMASK(:,:) = 0
where ( XMXL_mean(750:1900,100:600).gt.1.0E05 ) XMASK(750:1900,100:600) = 1

XMASK_area  = sum(DXT*DYT*XMASK)                                                
write (*,*) 'XMASK area: ', XMASK_area/1.0E04, ' m^2'

open(11,file=XMASK_file,access='direct',form='unformatted',                    &
        recl=imt*jmt,status='unknown')
write(11,rec=1) XMASK

!===============================================================================
! calculating masked quantitites
!===============================================================================

allocate( TEMP(imt,jmt), TEMP_avg(km), TEMP_avg2(km) )

open(12,file=HMXL_ymask_file,access='direct',form='unformatted',               &
        recl=1,status='unknown')
open(13,file=XMXL_ymask_file,access='direct',form='unformatted',               &
        recl=1,status='unknown')
open(14,file=TMXL_ymask_file,access='direct',form='unformatted',               &
        recl=1,status='unknown')
open(15,file=XMXL_ymask_max_file,access='direct',form='unformatted',           &
        recl=1,status='unknown')

open(16,file=HMXL_mmask_file,access='direct',form='unformatted',               &
        recl=1,status='unknown')
open(17,file=XMXL_mmask_file,access='direct',form='unformatted',               &
        recl=1,status='unknown')
open(18,file=TMXL_mmask_file,access='direct',form='unformatted',               &
        recl=1,status='unknown')
open(19,file=XMXL_mmask_max_file,access='direct',form='unformatted',           &
        recl=1,status='unknown')

open(20,file=TEMP_mmask_file,access='direct',form='unformatted',           &
        recl=km,status='unknown')
open(21,file=TEMP_ymask_file,access='direct',form='unformatted',           &
        recl=km,status='unknown')

!same loops as above
! loop over years
counter = 0
do y=start_year,end_year

  if ( y<276 ) then
    file_base = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc/tavg/t.t0.1_42l_nccs01.'
  elseif ( y>=276 ) then
    file_base = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_viebahn/tavg/t.t0.1_42l_nccs01.'
  endif

  !=============================================================================
  ! calculating masked quantitites yearly
  !=============================================================================
  read(2,rec=y-start_year+1) HMXL
  read(3,rec=y-start_year+1) XMXL
  read(4,rec=y-start_year+1) TMXL

  ! average in masked area
  write(12,rec=y-start_year+1) real(sum(HMXL*XMASK*TAREA)/XMASK_area)
  ! average in masked area
  write(13,rec=y-start_year+1) real(sum(XMXL*XMASK*TAREA)/XMASK_area)
  ! average in masked area
  write(14,rec=y-start_year+1) real(sum(TMXL*XMASK*TAREA)/XMASK_area)
  ! maximum in masked area
  write(15,rec=y-start_year+1) real(maxval(XMXL*XMASK))

  !=============================================================================
  ! calculating masked quantitites monthly
  !=============================================================================
  ! loop over months
  do m=1,12
    counter = counter + 1 
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

    cmd = "stat --format='%s.0' "//trim(file_name)//" > size"
    call system(cmd) 
    open(11,file='size')
    read(11,*) file_size
    cmd = 'rm size'
    call system(cmd)
    close(11)

    write(*,*) y, file_name, file_size
    if     ( file_size==11923200000.0 ) then
      nrec_HMXL =  343
      nrec_XMXL =  344
      nrec_TMXL =  345
      nrec_TEMP =  127
    elseif ( file_size==14826240000.0 ) then 
      nrec_HMXL =  427
      nrec_XMXL =  428
      nrec_TMXL =  429
      nrec_TEMP =  127
    elseif ( file_size==59132160000.0 ) then 
      nrec_HMXL = 1709
      nrec_XMXL = 1710
      nrec_TMXL = 1711
      nrec_TEMP =  222
    endif

    read(1,rec=nrec_HMXL) HMXL
    read(1,rec=nrec_XMXL) XMXL
    read(1,rec=nrec_TMXL) TMXL

    ! as above
    write(16,rec=counter) real(sum(HMXL*XMASK*TAREA)/XMASK_area)
    write(17,rec=counter) real(sum(XMXL*XMASK*TAREA)/XMASK_area)
    write(18,rec=counter) real(sum(TMXL*XMASK*TAREA)/XMASK_area)
    write(19,rec=counter) real(maxval(XMXL*XMASK))

    !===========================================================================
    ! TEMP profiles
    !===========================================================================
    ! monthly
    do k=1,km
      read(1,rec=nrec_TEMP+k-1) TEMP
      call masked_avg(imt,jmt,DXT,DYT,XMASK,XMASK_area,TEMP,TEMP_avg(k))
    enddo
    write(20,rec=counter) TEMP_avg(:)

    ! yearly average TEMP profile
    if ( m==12 ) then
      TEMP_avg(:) = 0.0
      do t=1,12
        read(20,rec=counter-t+1) TEMP_avg2
        TEMP_avg(:) = TEMP_avg(:) + TEMP_avg2(:)/12.0
      enddo !t
      write(21,rec=y-start_year+1) TEMP_avg(:)
    endif
  enddo !m
enddo !y


!===============================================================================
!===============================================================================
end program
