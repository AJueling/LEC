program extra_analysis
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

! user input 
integer ::                                                                     &
  opt1,opt2,opt3,opt4,                                                         &
  imt,jmt,km,ntavg,                                                            &
  y, start_year, end_year, nt,                                                 &
  nrec_rPm, nrec_rPe, nrec_rKm, nrec_rKe,                                      &
  nrec_gPm, nrec_gPe, nrec_gKm, nrec_gKe,                                      &
  nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                                  &
  nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                                  &
  nrec_TEMP, nrec_VVEL, nrec_WVEL,                                             &
  nrec_HMXL, nrec_XMXL, nrec_TMXL

character*120 :: grid_file,kmt_file,in_depths,pbc_file
character*38  :: LEC_file
character*58  :: bin_file
character*3   :: year, yr
character*62  :: filename1
character*42  :: filename2

! internal model variables 
integer                                         :: rec_length,i,j,k
integer,          dimension(:,:),   allocatable :: kmT

double precision                                :: volume 
double precision, dimension(:),     allocatable :: dz,area
double precision, dimension(:,:),   allocatable ::                             &
  HTN, HTE, DXT, DYT, TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK, WORK2, WORK3  

real                                            ::                             &
  rPm_int, rPe_int, rKm_int, rKe_int
real,             dimension(:),     allocatable ::                             &
  tdepth, cPKm_zint, cPKe_zint,                                                &
  vrPm_zint, vrPe_zint, vrKm_zint, vrKe_zint,                                  &
  cPKm_rsum, cPKe_rsum,                                                        &
  vrPm_rsum, vrPe_rsum, vrKm_rsum, vrKe_rsum

real,             dimension(:,:),   allocatable ::                             &
  cPKm_vint, cPKe_vint, cPKm_avg, cPKe_avg, Psi_v, Psi_w, Psi_v_avg, Psi_w_avg,&
  vrPm, vrPe, vrKm, vrKe, vrPm_vint, vrPe_vint, vrKm_vint, vrKe_vint,          &
  vrPm_avg, vrPe_avg, vrKm_avg, vrKe_avg,                                      &
  ref_state

real,             dimension(:,:,:), allocatable ::                             &
  DZT, DZU, cPKm, cPKe, rPm, rPe, rKm, rKe, VVEL, WVEL, TTT_VVEL,              &
  cPKm_yrly, cPKe_yrly, vrPm_yrly, vrPe_yrly, vrKm_yrly, vrKe_yrly,            &
  Psi_v_yrly, Psi_w_yrly, HMXL_yrly, XMXL_yrly, TMXL_yrly

double precision, parameter ::                                                 &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,                          &
S0ref = 0.035, rho0 = 4.1/3.996*1000, c = 3996,                                &
g = 9.806, salinity_factor = -34.7*1.0E-04, RE = 6.37E06

write (*,*) ''
write (*,*) '--- EXTRA ANALYSIS ---'
write (*,*) ''

!===============================================================================
!  INPUT
!===============================================================================

opt1=1
opt3=1
ntavg=5

imt            = 3600
jmt            = 2400
km             =   42

input_folder   = '/home/dijkbio/andre/LEC/input/'
in_depths      = trim(input_folder)//'in_depths.42.dat'
geometry_file  = trim(input_folder)//'geometry'

! read geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km),         &
          ref_state(15,km) )
open(1,file=geometry_file,access='direct',form='unformatted',recl=imt*jmt,     &
       status='old')
open(2,file=ref_state_file,access='direct',form='unformatted',recl=15,         &
       status='old')
! k, dz, tdepth, area, p,
! <T>, <S>, <RHO>, <PD>, <D(D,T,0)>, <D(D,T,p)>, <Q>,
! dz(<PD>), dz(<D(S,T,0)>), dz(<D(S,T,p)>)
read(1,rec=1) DXT ! [m]
read(1,rec=2) DXU ! [m]
read(1,rec=3) DYT ! [m]
read(1,rec=4) DYU ! [m]
read(1,rec=5) TAREA ! [m^2]
read(1,rec=6) UAREA ! [m^2]
do k=1,km
  read(1,rec=6+   k) DZT(:,:,k) ! [m]
  read(1,rec=6+km+k) DZU(:,:,k) ! [m]
  read(2,rec=k)      ref_state(:,k) ! [varying]
  dz(k)            = ref_state(2,k) ! [m]
  tdepth(k)        = ref_state(3,k) ! [m]
  area(k)          = ref_state(4,k) ! [m^2]
enddo
close(1)
close(2)

read  (*,'(a58)')  bin_file
read  (*,'(a38)')  LEC_file

nrec_rPm  =    1
nrec_rPe  =   43
nrec_rKm  =   85
nrec_rKe  =  127
nrec_gPm  =  169
nrec_gPe  =  170
nrec_gKm  =  171
nrec_gKe  =  172
nrec_gPmh =  173
nrec_gPms =  174
nrec_gPeh =  175
nrec_gPes =  176
nrec_cPem =  177
nrec_cKem =  219
nrec_cPKm =  261
nrec_cPKe =  303

nrec_VVEL =   43

!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'

! years
if ( ntavg==1 ) then
  start_year = 276
  end_year   = 326
elseif ( ntavg==5 ) then
  start_year = 278
  end_year   = 322
elseif ( ntavg==11 ) then
  start_year = 281
  end_year   = 315 ! check this one again
endif
nt = end_year-start_year+1

! naming the output files
write(yr,"(I3)") ntavg
yr = adjustl(yr)

if ( opt1==1 ) then
  write(*,*) 'option 1: cPKm/cPKe terms'
  open( 3,file='/projects/0/samoc/jan/Andree/cPKm_'//trim(yr),                 &
         form='unformatted',access='direct',recl=imt*jmt)
  open( 4,file='/projects/0/samoc/jan/Andree/cPKe_'//trim(yr),                 &
         form='unformatted',access='direct',recl=imt*jmt)
  open(11,file='/projects/0/samoc/jan/Andree/cPKm_zint_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(12,file='/projects/0/samoc/jan/Andree/cPKe_zint_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(25,file='/projects/0/samoc/jan/Andree/cPKm_rsum_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(26,file='/projects/0/samoc/jan/Andree/cPKe_rsum_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
endif
if ( opt3==1 ) then
  write(*,*) 'option 3: meridional advection of energy reservoirs'
  open( 7,file='/projects/0/samoc/jan/Andree/vrPm_'//trim(yr),                 &
         form='unformatted',access='direct',recl=imt*jmt)
  open( 8,file='/projects/0/samoc/jan/Andree/vrPe_'//trim(yr),                 &
         form='unformatted',access='direct',recl=imt*jmt)
  open( 9,file='/projects/0/samoc/jan/Andree/vrKm_'//trim(yr),                 &
         form='unformatted',access='direct',recl=imt*jmt)
  open(10,file='/projects/0/samoc/jan/Andree/vrKe_'//trim(yr),                 &
         form='unformatted',access='direct',recl=imt*jmt)
  open(13,file='/projects/0/samoc/jan/Andree/vrPm_zint_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(14,file='/projects/0/samoc/jan/Andree/vrPe_zint_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(15,file='/projects/0/samoc/jan/Andree/vrKm_zint_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(16,file='/projects/0/samoc/jan/Andree/vrKe_zint_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(27,file='/projects/0/samoc/jan/Andree/vrPm_rsum_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(28,file='/projects/0/samoc/jan/Andree/vrPe_rsum_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(29,file='/projects/0/samoc/jan/Andree/vrKm_rsum_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
  open(30,file='/projects/0/samoc/jan/Andree/vrKe_rsum_'//trim(yr),            &
         form='unformatted',access='direct',recl=jmt)
endif


do y = start_year, end_year

  ! open file
  write(year,"(I3)") y
  filename1 = bin_file//year       ! original binary POP output file
  filename2 = LEC_file//'_'//year  ! LEC.f90 output file
  !write (*,*) filename1
  !write (*,*) filename2
  open (1,file=filename1,access='direct',form='unformatted',recl=imt*jmt,      &
        status='unknown')
  open (2,file=filename2,access='direct',form='unformatted',recl=imt*jmt,      &
        status='unknown')

  ! loading VVEL, WVEL, TTT_VVEL
  allocate( VVEL(imt,jmt,km),   WVEL(imt,jmt,km), TTT_VVEL(imt,jmt,km) )
  call load_3D_field(imt,jmt,km,1,nrec_VVEL,VVEL)
  call load_3D_field(imt,jmt,km,1,nrec_WVEL,WVEL)
  call uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,VVEL,TTT_VVEL)
  !=============================================================================
  !  1. cPKm/cPKe
  !=============================================================================

  if ( opt1==1 ) then
    allocate( cPKm(imt,jmt,km),   cPKe(imt,jmt,km),                            &
              cPKm_vint(imt,jmt), cPKe_vint(imt,jmt),                          &
              cPKm_zint(jmt),     cPKe_zint(jmt),                              &
              cPKm_rsum(jmt),     cPKe_rsum(jmt) )
  
    ! load fields
    call load_3D_field(imt,jmt,km,2,nrec_cPKm,cPKm)
    call load_3D_field(imt,jmt,km,2,nrec_cPKe,cPKe)
  
    !vertical integral
    call vert_int(cPKm,DZT,cPKm_vint)
    call vert_int(cPKe,DZT,cPKe_vint)
  
    ! zonal-depth integral
    call zonal_int(imt,jmt,DXT,cPKm_vint,cPKm_zint)
    call zonal_int(imt,jmt,DXT,cPKe_vint,cPKe_zint)
  
    ! write into interal fields cPKm_vint
    write ( 3,rec=y-start_year+1) cPKm_vint(:,:)
    write ( 4,rec=y-start_year+1) cPKe_vint(:,:)
    write(*,*) 'option 1e'
    write (11,rec=y-start_year+1) cPKm_zint(:)
    write (12,rec=y-start_year+1) cPKe_zint(:)

    deallocate( cPKm, cPKe, cPKm_vint, cPKe_vint,                              &
                cPKm_zint, cPKe_zint, cPKm_rsum, cPKe_rsum )
  endif
  
  !=============================================================================
  !  3. Meridional Advection Of Energy
  !=============================================================================
  
  if ( opt3==1 ) then
    allocate( rPm(imt,jmt,km),    rPe(imt,jmt,km),                             &
              rKm(imt,jmt,km),    rKe(imt,jmt,km),                             &
              vrPm_vint(imt,jmt), vrPe_vint(imt,jmt),                          &
              vrKm_vint(imt,jmt), vrKe_vint(imt,jmt),                          &
              vrPm_zint(jmt),     vrPe_zint(jmt),                              &
              vrKm_zint(jmt),     vrKe_zint(jmt),                              &
              vrPm_rsum(jmt),     vrPe_rsum(jmt),                              &
              vrKm_rsum(jmt),     vrKe_rsum(jmt) )
  
    ! load fields
    call load_3D_field(imt,jmt,km,2,nrec_rPm, rPm )
    call load_3D_field(imt,jmt,km,2,nrec_rPe, rPe )
    call load_3D_field(imt,jmt,km,2,nrec_rKm, rKm )
    call load_3D_field(imt,jmt,km,2,nrec_rKe, rKE )

    call vol_int(1,1,1,imt,jmt,km,rPm,TAREA,DZT,rPm_int)
    call vol_int(1,1,1,imt,jmt,km,rPe,TAREA,DZT,rPe_int)
    call vol_int(1,1,1,imt,jmt,km,rKm,TAREA,DZT,rKm_int)
    call vol_int(1,1,1,imt,jmt,km,rKe,TAREA,DZT,rKe_int)
    write (*,*) rPm_int, rPe_int, rKm_int, rKe_int
  
    ! calculate vertically integrated, meridional advection
    call mer_advection(imt,jmt,km,DXT,DZT,TTT_VVEL,rPm,vrPm_vint)
    call mer_advection(imt,jmt,km,DXT,DZT,TTT_VVEL,rPe,vrPe_vint)
    call mer_advection(imt,jmt,km,DXT,DZT,TTT_VVEL,rKm,vrKm_vint)
    call mer_advection(imt,jmt,km,DXT,DZT,TTT_VVEL,rKe,vrKe_vint)

    ! calculate zonal-depth integral
    call zonal_int(imt,jmt,DXT,vrPm_vint,vrPm_zint)
    call zonal_int(imt,jmt,DXT,vrPe_vint,vrPe_zint)
    call zonal_int(imt,jmt,DXT,vrKm_vint,vrKm_zint)
    call zonal_int(imt,jmt,DXT,vrKe_vint,vrKe_zint)
  
    ! write to output file
    write ( 7,rec=y-start_year+1) vrPm_vint(:,:)
    write ( 8,rec=y-start_year+1) vrPe_vint(:,:)
    write ( 9,rec=y-start_year+1) vrKm_vint(:,:)
    write (10,rec=y-start_year+1) vrKe_vint(:,:)
    write (13,rec=y-start_year+1) vrPm_zint(:)
    write (14,rec=y-start_year+1) vrPe_zint(:)
    write (15,rec=y-start_year+1) vrKm_zint(:)
    write (16,rec=y-start_year+1) vrKe_zint(:)
  
    deallocate( rPm,       rPe,       rKm,       rKe,                          &
                vrPm_vint, vrPe_vint, vrKm_vint, vrKe_vint,                    &
                vrPm_zint, vrPe_zint, vrKm_zint, vrKe_zint,                    &
                vrPm_rsum, vrPe_rsum, vrKm_rsum, vrKe_rsum )

  endif
  
  close(1)
  close(2)

enddo ! y

!=============================================================================
!  Averages
!=============================================================================

if ( opt1==1 ) then
  allocate( cPKm_yrly(imt,jmt,ntavg), cPKe_yrly(imt,jmt,ntavg),                &
            cPKm_avg(imt,jmt),        cPKe_avg(imt,jmt),                       &
            cPKm_zint(jmt),           cPKe_zint(jmt),                          &
            cPKm_rsum(jmt),           cPKe_rsum(jmt) )

  call read_yrly(imt,jmt,ntavg,3,cPKm_yrly)
  call read_yrly(imt,jmt,ntavg,4,cPKe_yrly)

  call average(imt,jmt,ntavg,cPKm_yrly,cPKm_avg)
  call average(imt,jmt,ntavg,cPKe_yrly,cPKe_avg)
  write ( 3,rec=nt+1) cPKm_avg(:,:)
  write ( 4,rec=nt+1) cPKe_avg(:,:)

  call zonal_int(imt,jmt,DXT,cPKm_avg,cPKm_zint)
  call zonal_int(imt,jmt,DXT,cPKe_avg,cPKe_zint)
  write (11,rec=nt+1) cPKm_zint(:)
  write (12,rec=nt+1) cPKe_zint(:)

  deallocate( cPKm_yrly, cPKe_yrly, cPKm_avg, cPKe_avg,                        &
              cPKm_zint, cPKe_zint, cPKm_rsum, cPKe_rsum)

  close( 3)
  close( 4)
  close(11)
  close(12)
  close(25)
  close(26)
endif


if (opt3==1 ) then
  allocate( vrPm_yrly(imt,jmt,ntavg), vrPe_yrly(imt,jmt,ntavg),                &
            vrKm_yrly(imt,jmt,ntavg), vrKe_yrly(imt,jmt,ntavg),                &
            vrPm_avg(imt,jmt),        vrPe_avg(imt,jmt),                       &
            vrKm_avg(imt,jmt),        vrKe_avg(imt,jmt),                       &
            vrPm_zint(jmt),           vrPe_zint(jmt),                          &
            vrKm_zint(jmt),           vrKe_zint(jmt),                          &
            vrPm_rsum(jmt),           vrPe_rsum(jmt),                          &
            vrKm_rsum(jmt),           vrKe_rsum(jmt) )

  call read_yrly(imt,jmt,ntavg, 7,vrPm_yrly)
  call read_yrly(imt,jmt,ntavg, 8,vrPe_yrly)
  call read_yrly(imt,jmt,ntavg, 9,vrKm_yrly)
  call read_yrly(imt,jmt,ntavg,10,vrKe_yrly)

  call average(imt,jmt,ntavg,vrPm_yrly,vrPm_avg)
  call average(imt,jmt,ntavg,vrPe_yrly,vrPe_avg)
  call average(imt,jmt,ntavg,vrKm_yrly,vrKm_avg)
  call average(imt,jmt,ntavg,vrKe_yrly,vrKe_avg)
  write ( 7,rec=nt+1) vrPm_avg(:,:)
  write ( 8,rec=nt+1) vrPe_avg(:,:)
  write ( 9,rec=nt+1) vrKm_avg(:,:)
  write (10,rec=nt+1) vrKe_avg(:,:)

  call zonal_int(imt,jmt,DXT,vrPm_avg,vrPm_zint)
  call zonal_int(imt,jmt,DXT,vrPe_avg,vrPe_zint)
  call zonal_int(imt,jmt,DXT,vrKm_avg,vrKm_zint)
  call zonal_int(imt,jmt,DXT,vrKe_avg,vrKe_zint)
  write (13,rec=nt+1) vrPm_zint(:)
  write (14,rec=nt+1) vrPe_zint(:)
  write (15,rec=nt+1) vrKm_zint(:)
  write (16,rec=nt+1) vrKe_zint(:)

  deallocate( vrPm_yrly, vrPe_yrly, vrKm_yrly, vrKe_yrly,                      &
              vrPm_avg,  vrPe_avg,  vrKm_avg,  vrKe_avg,                       &
              vrPm_zint, vrPe_zint, vrKm_zint, vrKe_zint,                      &
              vrPm_rsum, vrPe_rsum, vrKm_rsum, vrKe_rsum )  

  close( 7)
  close( 8)
  close( 9)
  close(10)
  close(13)
  close(14)
  close(15)
  close(16)
  close(27)
  close(28)
  close(29)
  close(30)
endif


!===============================================================================
!  SUBROUTINES
!===============================================================================

contains


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine rsum(imt,jmt,km,DXT,DYT,DZT,FIELD,rsum)
implicit none
!
!  calculates running sum
!

integer,                                 intent(in)  :: imt,jmt,km,opt
integer                                              :: j
real,             dimension(imt,jmt,km), intent(in)  :: DZT
if ( opt==1) then
  real,           dimension(imt,jmt,km), intent(in)  :: FIELD
elseif ( opt==2 ) then
  real,           dimension(imt,jmt),    intent(in)  :: FIELD
endif
real,             dimension(jmt),        intent(out) :: rsum
double precision, dimension(imt,jmt),    intent(in)  :: DXT,DYT

if ( opt==1 ) then
  rsum(:) = sum(sum(FIELD(:,:,:)*DZT(:,:,:),3)*DXT(:,:)*DYT(:,:),1)
  do j=2,jmt
    rsum(j) = rsum(j) + rsum(j-1)
  enddo
elseif ( opt==2 ) then
  rsum(:) = sum(FIELD(:,:)*DXT(:,:)*DYT(:,:),1)
  do j=2,jmt
    rsum(j) = rsum(j) + rsum(j-1)
  enddo
endif

end subroutine rsum

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine read_yrly(x,y,ntavg,nfile,FIELD)
implicit none
!
!  reads yearly data that was previoously calculated and written into a binary
!  file
!

integer                     , intent(in)  :: x,y,ntavg,nfile
integer                                   :: t
real,   dimension(x,y,ntavg), intent(out) :: FIELD

do t = 1,ntavg
  read (nfile,rec=t) FIELD(:,:,t)
enddo

end subroutine read_yrly

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine uu2tt_3D(imt,jmt,km,DZT,DZU,TAREA,UAREA,TUU_FIELD,TTT_FIELD)
implicit none
!
!  interpolates a 3D field from UU to TT
!

integer                                  :: imt,jmt,km
real, dimension(imt,jmt)                 :: WORK
double precision, dimension(imt,jmt),    intent(in)  :: TAREA, UAREA
real, dimension(imt,jmt,km), intent(in)  :: TUU_FIELD, DZU, DZT
real, dimension(imt,jmt,km), intent(out) :: TTT_FIELD

do k = 1,km
  call uu2tt(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,TUU_FIELD(:,:,k),WORK)
  TTT_FIELD(:,:,k) = WORK
enddo

end subroutine uu2tt_3D


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine zonal_int(imt,jmt,DXT,FIELD,FIELD_zint)
!
!  "zonal" integral, true for Southern Hemisphere, Northern one is distorted
!
implicit none

integer,                                intent(in)  :: imt,jmt
double precision,   dimension(imt,jmt), intent(in)  :: DXT
real,               dimension(imt,jmt), intent(in)  :: FIELD
real,               dimension(jmt),     intent(out) :: FIELD_zint

FIELD_zint(:) = sum(FIELD(:,:)*DYT(:,:),1)

end subroutine

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

subroutine mer_advection(imt,jmt,km,DXT,DZT,TTT_VVEL,Eres,Eres_adv)
!
!  meridional advection of a energy reservoir term
!
implicit none

integer,          intent(in)                         :: imt, jmt, km
real,             dimension(imt,jmt,km), intent(in)  :: DZT, TTT_VVEL, Eres
double precision, dimension(imt,jmt),    intent(in)  :: DXT
real,             dimension(imt,jmt),    intent(out) :: Eres_adv

Eres_adv(:,:) = sum(Eres(:,:,:)*TTT_VVEL(:,:,:)*DZT(:,:,:),3) * DXT(:,:)

end subroutine mer_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine osf(imt,jmt,km,DXU,DZU,VVEL,WVEL,Psi_v,Psi_w)
!
!  overturning stream function calculated both via VVEL and WVEL
!
implicit none

integer,          intent(in)                        :: imt, jmt, km
integer                                             :: l
real,             dimension(imt,jmt,km)             :: DZU, VVEL, WVEL
double precision, dimension(imt,jmt),   intent(in)  :: DXU
real,             dimension(jmt,km),    intent(out) :: Psi_v, Psi_w 


! preparing integral dz/dy + executing integral dx
do k = 1,km
  Psi_v(:,k) = -sum(VVEL(:,:,k)*DZU(:,:,k)*DXU(:,:)/1.0E06,1)
  Psi_w(:,k) =  sum(WVEL(:,:,k)*DYU(:,:)  *DXU(:,:)/1.0E06,1)
enddo

! Psi_v (depth integral dz bottom to depth z)
do k = 1,km-1
  Psi_v(:,km-k) = Psi_v(:,km-k) + Psi_v(:,km-k+1)
enddo

! Psi_w (meridional integral dy South Pole to latitude j)
do j = 2,jmt
  Psi_w(j,:) = Psi_w(j,:) + Psi_w(j-1,:)
enddo

end subroutine osf

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine load_3D_field(imt,jmt,km,nfile,nrec,FIELD)
!
!  loads 3D field (starting with record number nrec) from file with number nfile
!
implicit none

integer, intent(in)                         :: imt, jmt, km, nfile,  nrec
real,    dimension(imt,jmt)                 :: WORK
real,    dimension(imt,jmt,km), intent(out) :: FIELD

do k = 1,km
  read (nfile,rec=nrec+k-1) WORK
  FIELD(:,:,k) = WORK
enddo

end subroutine load_3d_field

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine load_2D_field(imt,jmt,nfile,nrec,FIELD)
!
!  loads 2D field (starting with record number nrec) from file with number nfile
!
implicit none

integer, intent(in)                      :: imt, jmt, nfile,  nrec
real,    dimension(imt,jmt), intent(out) :: FIELD

read (nfile,rec=nrec) FIELD

end subroutine load_2d_field

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine surf_int(imin,jmin,imax,jmax,FIELD,TAREA,MASK,INTEGRAL)
!
!     calculates surface integral (*[m^2]) for TT-surface, using Kahan Summation
!
implicit none

! input/output variables
integer,                                 intent(in)  :: imin,jmin,imax,jmax
real,                dimension(imt,jmt), intent(in)  :: FIELD, MASK!=DZT(:,:,1)
double precision,    dimension(imt,jmt), intent(in)  :: TAREA 
real,                                    intent(out) :: INTEGRAL
real                                                 :: c, y, t

INTEGRAL = 0.0
c = 0.0
do j = jmin,jmax
 do i = imin,imax
  if ( MASK(i,j).ne.0.0 ) then
   y = TAREA(i,j) * FIELD(i,j) - c
   t = INTEGRAL + y
   c = ( t - INTEGRAL ) - y
   INTEGRAL = t
  endif
 enddo !i
enddo !j

INTEGRAL = INTEGRAL/1.0E04

!INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)*FIELD(imin:imax,jmin:jmax),          &
!               MASK(imin:imax,jmin:jmax).ne.0.0)/1.0E04
!write (*,*) '2D: ', INTEGRAL
end subroutine surf_int

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine vert_int(FIELD,DZT,INTEGRAL)
!
!     calculates vertical integral (*[m]) for TT-surface, using Kahan Summation
!
implicit none

! input/output variables
real,    dimension(imt,jmt,km), intent(in)  :: FIELD
real,    dimension(imt,jmt,km), intent(in)  :: DZT 
real,    dimension(imt,jmt),    intent(out) :: INTEGRAL

INTEGRAL = sum(DZT*FIELD,3)/1.0e02
! DZT in [cm], factor 1e-2 used to achieve [m]

end subroutine vert_int

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
subroutine vol_int(imin,jmin,kmin,imax,jmax,kmax,FIELD,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field, using Kahan Summation
!
implicit none

! input/output variables
integer,             intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,                dimension(imt,jmt,km), intent(in)  :: DZT, FIELD
double precision,    dimension(imt,jmt   ), intent(in)  :: TAREA
real,                                       intent(out) :: INTEGRAL
real                                                    :: c, y, t

INTEGRAL = 0.0
c = 0.0
do k = kmin,kmax
  do j = jmin,jmax
    do i = imin,imax
      if ( DZT(i,j,k).ne.0.0 ) then
        y = TAREA(i,j) * DZT(i,j,k) * FIELD(i,j,k) - c
        t = INTEGRAL + y
        c = ( t - INTEGRAL ) - y
        INTEGRAL = t
      endif
    enddo !i
  enddo !j
enddo !k

INTEGRAL = INTEGRAL/1.0E06
! write (*,*) '3D: ', INTEGRAL

!INTEGRAL = sum(TAREA(imin:imax,jmin:jmax)                                      &
!         * sum(DZT(imin:imax,jmin:jmax,kmin:kmax)                            &
!         * FIELD(imin:imax,jmin:jmax,kmin:kmax),3))/1.0E06
write (*,*) '3D: ', INTEGRAL

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
subroutine wtt2ttt(imt,jmt,km,W_FIELD,DZT,T_FIELD)
!
!     interpolated real field from WTT-grid to TTT-grid
!     interpolating between two W-levels
!     T-point is located exactly in the middle of two W-points
!     T_FIELD(i,j,k) = W_FIELD(i,j,k) + ( W_FIELD(i,j,k+1) - W_FIELD(i,j,k) ) &
!                                       / dz(k) * dz(k) / 2
!
implicit none

! input/output variables
integer,                         intent(in)  :: imt,jmt,km
real,     dimension(imt,jmt,km), intent(in)  :: DZT
real,     dimension(imt,jmt,km), intent(in)  :: W_FIELD ! WTT
real,     dimension(imt,jmt,km), intent(out) :: T_FIELD ! TTT
! internal variables
double precision,  parameter                 :: p5=0.5, c0=0.0, c3=3.0

do k = 1,km-1
 do j = 1,jmt
  do i = 1,imt
   if ( DZT(i,j,k).ne.c0 ) then
    T_FIELD(i,j,k) = p5 * ( W_FIELD(i,j,k) + W_FIELD(i,j,k+1) )
   endif
  enddo !i
 enddo !j
enddo !k

! bottom level
do j = 1,jmt
  do i = 1,imt
    if ( DZT(i,j,km).ne.c0 ) then
      ! extrapolating to last level with gradient between last two W-levels
      ! T_FIELD(i,j,km) = W_FIELD(i,j,km) 
      !                 + (W_FIELD(i,j,km)-W_FIELD(i,j,km-1)) / dz(k) * dz(k)/2
      T_FIELD(i,j,km) = p5 * (c3*W_FIELD(i,j,km) - W_FIELD(i,j,km-1) )
    endif
  enddo !i
enddo !j

end subroutine

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine uu2tt(imt,jmt,k,DZT_k,DZU_k,TAREA,UAREA,FIELD_k,NEW_FIELD_k)
!
!     interpolates real 2D field from UU-grid to TT-grid on T-depth levels
!
implicit none

! input/output variables
integer,                              intent(in)  :: imt,jmt,k
real,             dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
double precision, dimension(imt,jmt), intent(in)  :: TAREA,UAREA
real,             dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
real,             dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT
real,             dimension(imt,jmt)              :: NEW_FIELD_k_old ! (T)TT
real                                              :: A, B 
! surface integral, used for comparing old and new method od interpolating 

! southernmost gridpoints (j=1) on TTT-grid must be 0, as there would not be 
! UU gridpoints south of it otherwise

do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    NEW_FIELD_k_old(1,j) = 0.25 * ( FIELD_k(1  ,j  )                       &
                                  + FIELD_k(imt,j  )                       &
                                  + FIELD_k(1  ,j-1)                       &
                                  + FIELD_k(imt,j-1) )
    NEW_FIELD_k(1,j) = 0.25                                                & 
                  * ( FIELD_k(1  ,j  ) * UAREA(1  ,j  ) * DZU_k(1  ,j  )   &
                    + FIELD_k(imt,j  ) * UAREA(imt,j  ) * DZU_k(imt,j  )   &
                    + FIELD_k(1  ,j-1) * UAREA(1  ,j-1) * DZU_k(1  ,j-1)   &
                    + FIELD_k(imt,j-1) * UAREA(imt,j-1) * DZU_k(imt,j-1) ) &
                  / TAREA(1,j) / DZT_k(1,j)
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      NEW_FIELD_k_old(i,j) = 0.25 * ( FIELD_k(i  ,j  )                      &
                                    + FIELD_k(i-1,j  )                      &
                                    + FIELD_k(i  ,j-1)                      &
                                    + FIELD_k(i-1,j-1) )
      NEW_FIELD_k(i,j) = 0.25                                               &
                   * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  ) * DZU_k(i  ,j  )   &
                     + FIELD_k(i-1,j  ) * UAREA(i-1,j  ) * DZU_k(i-1,j  )   &
                     + FIELD_k(i  ,j-1) * UAREA(i  ,j-1) * DZU_k(i  ,j-1)   &
                     + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) * DZU_k(i-1,j-1) ) &
                   / TAREA(i,j) / DZT_k(i,j)
    endif
  enddo !i
enddo !j

call surf_int(1,1,imt,jmt,NEW_FIELD_k_old,TAREA,DZT_k,A)
call surf_int(1,1,imt,jmt,NEW_FIELD_k    ,TAREA,DZT_k,B)
!write (*,*) 'uu2tt method error: ', k, (A-B)/B

end subroutine uu2tt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

subroutine uu2tt_scalar(imt,jmt,k,DZT_k,DZU_k,TAREA,UAREA,FIELD_k,NEW_FIELD_k)
!
!     interpolates real scalar 2D field (as opposed to velocity fields which  
!     can be related to volume transports) from UU-grid to TT-grid on T-depth 
!     levels
!
implicit none

! input/output variables
integer,                              intent(in)  :: imt,jmt,k
real,             dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
double precision, dimension(imt,jmt), intent(in)  :: TAREA,UAREA
real,             dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
real,             dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT
real :: A, B ! surface integral, used for comparing old/new interp. method

! southernmost gridpoints (j=1) on TTT-grid must be 0, as there would not
! be UU gridpoints south of it otherwise
do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    NEW_FIELD_k(1,j) = 0.25 * ( FIELD_k(1  ,j  ) * UAREA(1  ,j  )          &
                              + FIELD_k(imt,j  ) * UAREA(imt,j  )          &
                              + FIELD_k(1  ,j-1) * UAREA(1  ,j-1)          &
                              + FIELD_k(imt,j-1) * UAREA(imt,j-1) )        &
                            / TAREA(1,j)
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      NEW_FIELD_k(i,j) = 0.25 * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  )         &
                                + FIELD_k(i-1,j  ) * UAREA(i-1,j  )         &
                                + FIELD_k(i  ,j-1) * UAREA(i  ,j-1)         &
                                + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) )       &
                              / TAREA(i,j)
    endif
  enddo !i
enddo !j

! results in very small correction O(-4)

end subroutine uu2tt_scalar

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
  enddo
enddo

end subroutine e_rshift

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
  enddo
enddo

end subroutine s_rshift

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end program
