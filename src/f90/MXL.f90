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
  counter,km,t,rec_length,k,n1000,nrec_PD,opt1,opt2,opt3,opt4,opt5,now,nrec_SALT
integer, dimension(:,:), allocatable ::                                        &
  XMASK, WMASK, SO30M
character*160                        ::                                        & 
  cmd,file_base,file_name,output_folder1,output_folder2,geometry1_file,        &
  HMXL_file,XMXL_file,TMXL_file,HMXL_mean_file,XMXL_mean_file,TMXL_mean_file,  &
  HMXL_std_file,XMXL_std_file,TMXL_std_file,XMASK_file,WMASK_file,             &
  HMXL_ymask_file,XMXL_ymask_file,TMXL_ymask_file,TEMP_ymask_file,             &
  HMXL_mmask_file,XMXL_mmask_file,TMXL_mmask_file,TEMP_mmask_file,             &
  SALT_ymask_file, SALT_mmask_file, PD_ymask_file, PD_mmask_file,              &
  XMXL_ymask_max_file,XMXL_mmask_max_file,LEC_folder,                          &
  XMXL_WMASK_file,TEMP_WMASK_s_file,TEMP_WMASK_v_file,PD_WMASK_s_file,         &
  PD_WMASK_v_file,HMXL_WMASK_file, TMXL_WMASK_file,                            &
  XMXL_SO30_file,TEMP_SO30_s_file,TEMP_SO30_v_file,PD_SO30_s_file,             &
  PD_SO30_v_file,HMXL_SO30_file, TMXL_SO30_file, geometry2_file,               &
  TEMP_WMASK_d_file,TEMP_SO30_d_file,PD_WMASK_d_file,PD_SO30_d_file
character*4                          ::                                        &
  year
character*2                          ::                                        &
  month
real                                 ::                                        &
  file_size,XMASK_area,XMASK_vol_1000,WMASK_area,WMASK_vol_1000,SO30M_vol_1000,&
  TEMP_vol_1000_s, PD_vol_1000_s,integral,TEMP_vol_1000_w, PD_vol_1000_w,      &
  SO30M_area,SO30M_area_l22,WMASK_area_l22
real, dimension(:),   allocatable    ::                                        &
  TEMP_avg, TEMP_avg2, dz, SALT_avg, SALT_avg2, PD_avg, PD_avg2
real, dimension(:,:), allocatable    ::                                        &
  TAREA,HMXL,XMXL,TMXL,HMXL_avg,XMXL_max,TMXL_min,                             &
  HMXL_mean,XMXL_mean,TMXL_mean,HMXL_sum,XMXL_sum,TMXL_sum,                    &
  DXT,DYT,ref_state
real,dimension(:,:,:), allocatable   ::                                        &
  DZT, TEMP, PD, ONES, SALT
double precision, parameter                     ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

opt1 = 0
opt2 = 0
opt3 = 0
opt4 = 1

imt        = 3600
jmt        = 2400
km         =   42

start_year =   75
end_year   =  326
nt         = end_year-start_year+1
n          =   51

n1000      =   22

!===============================================================================
LEC_folder = '/home/dijkbio/andre/LEC/'
geometry1_file  = trim(LEC_folder)//'input/geometry1'
geometry2_file  = trim(LEC_folder)//'input/geometry2'
! read 2D geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          ref_state(6,km), dz(km) )
open(2,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,    &
       status='old')
open(3,file=geometry2_file,access='direct',form='unformatted',recl=6,    &
       status='old')

read(2,rec=1) DXT ! [m]
read(2,rec=3) DYT ! [m]
read(2,rec=5) TAREA ! [m^2]
do k=1,km
  read(2,rec=k+6) DZT(:,:,k) ! [m]
  read(3,rec=k) ref_state(:,k) ! [varying]
  dz(k)       = ref_state(2,k) ! [m]
enddo !k
close(2)
close(3)


!===============================================================================

output_folder1      = '/projects/0/samoc/jan/Andree/MXL/'
output_folder2      = '/home/dijkbio/andre/LEC/results/MXL/'

! (3D) average of the fields for every model year
HMXL_file           = trim(output_folder1)//'HMXL_yrly_avg'
XMXL_file           = trim(output_folder1)//'XMXL_yrly_avg'
TMXL_file           = trim(output_folder1)//'TMXL_yrly_avg'

! (2D) map of mean and std of fields of the last n=50 years
HMXL_mean_file      = trim(output_folder2)//'HMXL_mean_map'
XMXL_mean_file      = trim(output_folder2)//'XMXL_mean_map'
TMXL_mean_file      = trim(output_folder2)//'TMXL_mean_map'
HMXL_std_file       = trim(output_folder2)//'HMXL_std_map'
XMXL_std_file       = trim(output_folder2)//'XMXL_std_map'
TMXL_std_file       = trim(output_folder2)//'TMXL_std_map'

! (2D) convection mask XMASK (where XMXL>1km) and WMASK (Weddel area)
XMASK_file          = trim(output_folder2)//'XMASK_map'
WMASK_file          = trim(output_folder2)//'WMASK_map'

! (1D) monthly or yearly averages over masked area of different quantities
TMXL_ymask_file     = trim(output_folder2)//'TMXL_XMASK_y.csv'
TMXL_mmask_file     = trim(output_folder2)//'TMXL_XMASK_m.csv'
HMXL_ymask_file     = trim(output_folder2)//'HMXL_XMASK_y.csv'
HMXL_mmask_file     = trim(output_folder2)//'HMXL_XMASK_m.csv'
XMXL_ymask_file     = trim(output_folder2)//'XMXL_XMASK_y.csv'
XMXL_mmask_file     = trim(output_folder2)//'XMXL_XMASK_m.csv'
XMXL_ymask_max_file = trim(output_folder2)//'XMXL_XMASK_max_y.csv'
XMXL_mmask_max_file = trim(output_folder2)//'XMXL_XMASK_max_m.csv'
TEMP_ymask_file     = trim(output_folder2)//'TEMP_XMASK_y'
TEMP_mmask_file     = trim(output_folder2)//'TEMP_XMASK_m'
SALT_ymask_file     = trim(output_folder2)//'SALT_XMASK_y'
SALT_mmask_file     = trim(output_folder2)//'SALT_XMASK_m'
PD_ymask_file     = trim(output_folder2)//'PD_XMASK_y'
PD_mmask_file     = trim(output_folder2)//'PD_XMASK_m'

HMXL_WMASK_file     = trim(output_folder2)//'HMXL_WMASK_m.csv'
TMXL_WMASK_file     = trim(output_folder2)//'TMXL_WMASK_m.csv'
XMXL_WMASK_file     = trim(output_folder2)//'XMXL_WMASK_m.csv'
TEMP_WMASK_s_file   = trim(output_folder2)//'TEMP_WMASK_surface_m.csv'
TEMP_WMASK_v_file   = trim(output_folder2)//'TEMP_WMASK_vol_m.csv'
PD_WMASK_s_file     = trim(output_folder2)//'PD_WMASK_surface_m.csv'
PD_WMASK_v_file     = trim(output_folder2)//'PD_WMASK_vol_m.csv'

HMXL_SO30_file      = trim(output_folder2)//'HMXL_SO30_m.csv'
TMXL_SO30_file      = trim(output_folder2)//'TMXL_SO30_m.csv'
XMXL_SO30_file      = trim(output_folder2)//'XMXL_SO30_m.csv'
TEMP_SO30_s_file    = trim(output_folder2)//'TEMP_SO30_surface_m.csv'
TEMP_SO30_v_file    = trim(output_folder2)//'TEMP_SO30_vol_m.csv'
PD_SO30_s_file      = trim(output_folder2)//'PD_SO30_surface_m.csv'
PD_SO30_v_file      = trim(output_folder2)//'PD_SO30_vol_m.csv'

TEMP_WMASK_d_file   = trim(output_folder2)//'TEMP_WMASK_l22_m.csv'
TEMP_SO30_d_file    = trim(output_folder2)//'TEMP_SO30_l22_m.csv'
PD_WMASK_d_file     = trim(output_folder2)//'PD_WMASK_l22_m.csv'
PD_SO30_d_file      = trim(output_folder2)//'PD_SO30_l22_m.csv'

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
if (opt1==1) then
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

    write(*,*) y!, file_name, file_size
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
endif !opt1==1
!===============================================================================
! mean of last n years
!===============================================================================
open( 5,file=HMXL_mean_file,access='direct',form='unformatted',                &
        recl=imt*jmt,status='unknown')
open( 6,file=XMXL_mean_file,access='direct',form='unformatted',                &
        recl=imt*jmt,status='unknown')
open( 7,file=TMXL_mean_file,access='direct',form='unformatted',                &
        recl=imt*jmt,status='unknown')
if (opt2==1) then
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
else
  read(5,rec=1) HMXL_mean
  read(6,rec=1) XMXL_mean
  read(7,rec=1) TMXL_mean
endif !opt2==1
!===============================================================================
! standard deviation of last n years
!===============================================================================
if (opt1==1) then
open( 8,file=HMXL_std_file,access='direct',form='unformatted',                 &
        recl=imt*jmt,status='unknown')
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
write( 9,rec=1) sqrt(XMXL_sum/(n-1))
write(10,rec=1) sqrt(TMXL_sum/(n-1))

close( 8) ! std
close( 9)
close(10)

endif !opt1==1
!===============================================================================
! create mask
!===============================================================================

allocate( XMASK(imt,jmt), WMASK(imt,jmt), SO30M(imt,jmt), ONES(imt,jmt,km) )
! creating masks
SO30M(:,:) = 0
SO30M(:,1:866) = 1
WMASK(:,:) = 0
WMASK(750:1900,100:600) = 1
XMASK(:,:) = 0
where ( XMXL_mean(750:1900,100:600).gt.1.0E05 ) XMASK(750:1900,100:600) = 1

! mask areas
SO30M_area     = sum(TAREA*float(SO30M),DZT(:,:, 1).ne.0.0)
WMASK_area     = sum(TAREA*float(WMASK),DZT(:,:, 1).ne.0.0)
XMASK_area     = sum(TAREA*float(XMASK),DZT(:,:, 1).ne.0.0)
SO30M_area_l22 = sum(TAREA*float(SO30M),DZT(:,:,22).ne.0.0)
WMASK_area_l22 = sum(TAREA*float(WMASK),DZT(:,:,22).ne.0.0)
ONES  = 1

call vol_int_full(  1,  1,1, imt,866,n1000,ONES,TAREA,dz,DZT,SO30M_vol_1000)
call vol_int_full(750,100,1,1900,600,n1000,ONES,TAREA,dz,DZT,WMASK_vol_1000)

XMASK_vol_1000 = 0.0
do k=1,n1000 ! n1000 = 22 => 1022 m
  XMASK_vol_1000 = XMASK_vol_1000 + sum(TAREA(:,:)*DZT(:,:,k)*XMASK(:,:))
enddo

write (*,*) 'XMASK: area = ', XMASK_area, ' m^2; vol = ',               &
             XMASK_vol_1000, ' m^3'
write (*,*) 'WMASK: area = ', WMASK_area, ' m^2; vol = ',               &
             WMASK_vol_1000, ' m^3; area_l22 = ', WMASK_area_l22
write (*,*) 'SO30M: area = ', SO30M_area, ' m^2; vol = ',               &
             SO30M_vol_1000, ' m^3; area_l22 = ', SO30M_area_l22


open(11,file=WMASK_file,access='direct',form='unformatted',                    &
        recl=imt*jmt,status='unknown')
write(11,rec=1) WMASK
close(11)
open(11,file=XMASK_file,access='direct',form='unformatted',                    &
        recl=imt*jmt,status='unknown')
write(11,rec=1) XMASK


!===============================================================================
! calculating masked quantitites
!===============================================================================

allocate( TEMP(imt,jmt,km), TEMP_avg(km), TEMP_avg2(km),                       &
          SALT(imt,jmt,km), SALT_avg(km), SALT_avg2(km),                       &
          PD(imt,jmt,km),   PD_avg(km),   PD_avg2(km) )

open(12,file=HMXL_ymask_file,    form='formatted')
open(13,file=XMXL_ymask_file,    form='formatted')
open(14,file=TMXL_ymask_file,    form='formatted')
open(15,file=XMXL_ymask_max_file,form='formatted')

open(16,file=HMXL_mmask_file,    form='formatted')
open(17,file=XMXL_mmask_file,    form='formatted')
open(18,file=TMXL_mmask_file,    form='formatted')
open(19,file=XMXL_mmask_max_file,form='formatted')

open(20,file=TEMP_mmask_file,access='direct',form='unformatted',           &
        recl=km,status='unknown')
open(21,file=TEMP_ymask_file,access='direct',form='unformatted',           &
        recl=km,status='unknown')

open(22,file=XMXL_WMASK_file,    form='formatted')
open(23,file=TEMP_WMASK_s_file,  form='formatted')
open(24,file=TEMP_WMASK_v_file,  form='formatted')
open(25,file=PD_WMASK_s_file,    form='formatted')
open(26,file=PD_WMASK_v_file,    form='formatted')
open(27,file=HMXL_WMASK_file,    form='formatted')
open(28,file=TMXL_WMASK_file,    form='formatted')

open(29,file=XMXL_SO30_file,     form='formatted')
open(30,file=TEMP_SO30_s_file,   form='formatted')
open(31,file=TEMP_SO30_v_file,   form='formatted')
open(32,file=PD_SO30_s_file,     form='formatted')
open(33,file=PD_SO30_v_file,     form='formatted')
open(34,file=HMXL_SO30_file,     form='formatted')
open(35,file=TMXL_SO30_file,     form='formatted')

open(36,file=TEMP_WMASK_d_file,  form='formatted')
open(37,file=TEMP_SO30_d_file,  form='formatted')
open(38,file=PD_WMASK_d_file,    form='formatted')
open(39,file=PD_SO30_d_file,    form='formatted')

open(40,file=SALT_mmask_file,access='direct',form='unformatted',               &
        recl=km,status='unknown')
open(41,file=SALT_ymask_file,access='direct',form='unformatted',               &
        recl=km,status='unknown')
open(42,file=PD_mmask_file,access='direct',form='unformatted',                 &
        recl=km,status='unknown')
open(43,file=PD_ymask_file,access='direct',form='unformatted',                 &
        recl=km,status='unknown')
!same loops as above
! loop over years
201 FORMAT (A,",",A)
202 FORMAT (I3,E15.7)
301 FORMAT (A,",",A,",",A)
302 FORMAT (I3","I2","E15.7)

write(12,201) 'year','HMXL XMASK yrly'
write(13,201) 'year','XMXL XMASK yrly'
write(14,201) 'year','TMXL XMASK yrly'
write(15,201) 'year','XMXL XMASK yrly max'

write(16,301) 'year','month','HMXL XMASK mnthly'
write(17,301) 'year','month','XMXL XMASK mnthly'
write(18,301) 'year','month','TMXL XMASK mnthly'
write(19,301) 'year','month','XMXL XMASK mnthly max'

write(22,301) 'year','month','XMXL WMASK mnthly'
write(27,301) 'year','month','HMXL WMASK mnthly'
write(28,301) 'year','month','TMXL WMASK mnthly'

write(23,301) 'year','month','TEMP WMASK surface mnthly'
write(24,301) 'year','month','TEMP WMASK top1000m mnthly'
write(25,301) 'year','month','PD WMASK surface mnthly'
write(26,301) 'year','month','PD WMASK top1000m mnthly'

write(29,301) 'year','month','XMXL SO30 mnthly'
write(34,301) 'year','month','HMXL SO30 mnthly'
write(35,301) 'year','month','TMXL SO30 mnthly'

write(30,301) 'year','month','TEMP SO30 surface mnthly'
write(31,301) 'year','month','TEMP SO30 top1000m mnthly'
write(32,301) 'year','month','PD SO30 surface mnthly'
write(33,301) 'year','month','PD SO30 top1000m mnthly'

write(36,301) 'year','month','TEMP WMASK layer22 mnthly'
write(37,301) 'year','month','TEMP SO30 layer22 mnthly'
write(38,301) 'year','month','PD WMASK layer22 mnthly'
write(39,301) 'year','month','PD SO30 layer22 mnthly'

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
  write(12,202) y,real(sum(HMXL*XMASK*TAREA)/XMASK_area)
  ! average in masked area
  write(13,202) y,real(sum(XMXL*XMASK*TAREA)/XMASK_area)
  ! average in masked area
  write(14,202) y,real(sum(TMXL*XMASK*TAREA)/XMASK_area)
  ! maximum in masked area
  write(15,202) y,real(maxval(XMXL*XMASK))

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

    write(*,*) y, m , file_name, file_size
    if     ( file_size==11923200000.0 ) then
      nrec_HMXL =  343
      nrec_XMXL =  344
      nrec_TMXL =  345
      nrec_TEMP =  127
      nrec_SALT =  169
      nrec_PD   =  301
    elseif ( file_size==14826240000.0 ) then  ! 14GB files before 277
      nrec_HMXL =  427
      nrec_XMXL =  428
      nrec_TMXL =  429
      nrec_TEMP =  127
      nrec_SALT =  169 
      nrec_PD   =  301
    elseif ( file_size==59132160000.0 ) then  ! 56GB extra_vars files after 276 
      nrec_HMXL = 1709
      nrec_XMXL = 1710
      nrec_TMXL = 1711
      nrec_TEMP =  222
      nrec_SALT =  264
      nrec_PD   =  385
    endif

    read(1,rec=nrec_HMXL) HMXL
    read(1,rec=nrec_XMXL) XMXL
    read(1,rec=nrec_TMXL) TMXL


    ! as above
    write(16,302) y,m,real(sum(HMXL*XMASK*TAREA)/XMASK_area)
    write(17,302) y,m,real(sum(XMXL*XMASK*TAREA)/XMASK_area)
    write(18,302) y,m,real(sum(TMXL*XMASK*TAREA)/XMASK_area)
    write(19,302) y,m,real(maxval(XMXL*XMASK))

    ! MXL in WMASK, SO30

    call surf_int_2D(750,100,1900,600,XMXL,TAREA,DZT(:,:,1),integral)
    write(22,302) y,m,integral/WMASK_area
    call surf_int_2D(750,100,1900,600,HMXL,TAREA,DZT(:,:,1),integral)
    write(27,302) y,m,integral/WMASK_area
    call surf_int_2D(750,100,1900,600,TMXL,TAREA,DZT(:,:,1),integral)
    write(28,302) y,m,integral/WMASK_area

    call surf_int_2D(1,1,imt,866,XMXL,TAREA,DZT(:,:,1),integral)
    write(29,302) y,m,integral/SO30M_area
    call surf_int_2D(1,1,imt,866,HMXL,TAREA,DZT(:,:,1),integral)
    write(34,302) y,m,integral/SO30M_area
    call surf_int_2D(1,1,imt,866,TMXL,TAREA,DZT(:,:,1),integral)
    write(35,302) y,m,integral/SO30M_area

    call load_3D_field(1,nrec_TEMP,TEMP)
    call load_3D_field(1,nrec_SALT,SALT)
    call load_3D_field(1,nrec_PD,PD)
    
    ! surface integrals
    call surf_int_2D(750,100,1900,600,TEMP(:,:,1),TAREA,DZT(:,:,1),integral)
    write(23,302) y,m,integral/WMASK_area
    call surf_int_2D(1,1,imt,866,TEMP(:,:,1),TAREA,DZT(:,:,1),integral)
    write(30,302) y,m,integral/SO30M_area

    call surf_int_2D(750,100,1900,600,TEMP(:,:,22),TAREA,DZT(:,:,22),integral)
    write(36,302) y,m,integral/WMASK_area_l22
    call surf_int_2D(1,1,imt,866,TEMP(:,:,22),TAREA,DZT(:,:,22),integral)
    write(37,302) y,m,integral/SO30M_area_l22

    call surf_int_2D(750,100,1900,600,PD(:,:,1),TAREA,DZT(:,:,1),integral)
    write(25,302) y,m,integral/WMASK_area
    call surf_int_2D(1,1,imt,866,PD(:,:,1),TAREA,DZT(:,:,1),integral)
    write(32,302) y,m,integral/SO30M_area

    call surf_int_2D(750,100,1900,600,PD(:,:,22),TAREA,DZT(:,:,22),integral)
    write(38,302) y,m,integral/WMASK_area_l22
    call surf_int_2D(1,1,imt,866,PD(:,:,22),TAREA,DZT(:,:,22),integral)
    write(39,302) y,m,integral/SO30M_area_l22

    ! volume integrals
    call vol_int_full(750,100,1,1900,600,n1000,TEMP,TAREA,dz,DZT,TEMP_vol_1000_w)
    call vol_int_full(750,100,1,1900,600,n1000,PD  ,TAREA,dz,DZT,PD_vol_1000_w)
    call vol_int_full(  1,  1,1, imt,866,n1000,TEMP,TAREA,dz,DZT,TEMP_vol_1000_s)
    call vol_int_full(  1,  1,1, imt,866,n1000,PD  ,TAREA,dz,DZT,PD_vol_1000_s)

    write(24,302) y,m,TEMP_vol_1000_w/WMASK_vol_1000
    write(26,302) y,m,PD_vol_1000_w  /WMASK_vol_1000
    write(31,302) y,m,TEMP_vol_1000_s/SO30M_vol_1000
    write(33,302) y,m,PD_vol_1000_s  /SO30M_vol_1000

    !===========================================================================
    ! TEMP profiles
    !===========================================================================
    ! monthly
    do k=1,km
      call masked_avg(DXT,DYT,XMASK,XMASK_area,TEMP(:,:,k),TEMP_avg(k))
      call masked_avg(DXT,DYT,XMASK,XMASK_area,SALT(:,:,k),SALT_avg(k))
      call masked_avg(DXT,DYT,XMASK,XMASK_area,PD(:,:,k)  ,PD_avg(k))
    enddo
    write(20,rec=counter) TEMP_avg(:)
    write(40,rec=counter) SALT_avg(:)
    write(42,rec=counter) PD_avg(:)

    ! yearly average TEMP profile
    if ( m==12 ) then
      TEMP_avg(:) = 0.0
      SALT_avg(:) = 0.0
      PD_avg(:) = 0.0
      do t=1,12
        read(20,rec=counter-t+1) TEMP_avg2
        read(40,rec=counter-t+1) SALT_avg2
        read(42,rec=counter-t+1) PD_avg2
        TEMP_avg(:) = TEMP_avg(:) + TEMP_avg2(:)/12.0
        SALT_avg(:) = SALT_avg(:) + SALT_avg2(:)/12.0
        PD_avg(:)   = PD_avg(:)   + PD_avg2(:)/12.0
      enddo !t
      write(21,rec=y-start_year+1) TEMP_avg(:)
      write(41,rec=y-start_year+1) SALT_avg(:)
      write(43,rec=y-start_year+1) PD_avg(:)
    endif
  enddo !m
enddo !y


!===============================================================================
!===============================================================================
end program
