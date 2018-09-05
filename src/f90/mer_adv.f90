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

character*120 ::&
grid_file,kmt_file,in_depths,pbc_file,geometry1_file,geometry2_file,input_folder
character*38  :: LEC_file
character*58  :: bin_file
character*3   :: year, yr
character*62  :: filename1
character*42  :: filename2

! internal model variables 
integer                                         :: rec_length,i,j,k
integer,          dimension(:,:),   allocatable :: kmT

real, dimension(:),     allocatable :: dz,area
real, dimension(:,:),   allocatable :: DXT, DYT, TAREA, DXU, DYU, UAREA  

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
  Psi_v_yrly, Psi_w_yrly, HMXL_yrly, XMXL_yrly, TMXL_yrly, FBC

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
geometry1_file  = trim(input_folder)//'geometry1'
geometry2_file  = trim(input_folder)//'geometry2'

! read geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km),         &
          dz(km),       area(km),     ref_state(6,km), FBC(imt,jmt,km) )
open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,     &
       status='old')
open(2,file=geometry2_file,access='direct',form='unformatted',recl=6,         &
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
  write (*,*) k
  read(1,rec=6+   k) DZT(:,:,k) ! [m]
  read(1,rec=6+km+k) DZU(:,:,k) ! [m]
  read(2,rec=k)      ref_state(:,k) ! [varying]
  dz(k)            = ref_state(2,k) ! [m]
  area(k)          = ref_state(4,k) ! [m^2]
  FBC(:,:,k)       = dz(k)
enddo
write (*,*) 'done'

close(1)
close(2)

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
nrec_WVEL =  405

!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'
write (*,*) opt1,opt3


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
!  open(25,file='/projects/0/samoc/jan/Andree/cPKm_rsum_'//trim(yr),            &
!         form='unformatted',access='direct',recl=jmt)
!  open(26,file='/projects/0/samoc/jan/Andree/cPKe_rsum_'//trim(yr),            &
!         form='unformatted',access='direct',recl=jmt)
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
!  open(27,file='/projects/0/samoc/jan/Andree/vrPm_rsum_'//trim(yr),            &
!         form='unformatted',access='direct',recl=jmt)
!  open(28,file='/projects/0/samoc/jan/Andree/vrPe_rsum_'//trim(yr),            &
!         form='unformatted',access='direct',recl=jmt)
!  open(29,file='/projects/0/samoc/jan/Andree/vrKm_rsum_'//trim(yr),            &
!         form='unformatted',access='direct',recl=jmt)
!  open(30,file='/projects/0/samoc/jan/Andree/vrKe_rsum_'//trim(yr),            &
!         form='unformatted',access='direct',recl=jmt)
endif


do y = start_year, end_year

  ! open file
  write(year,"(I3)") y
  filename1 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.5year.'//year       ! original binary POP output file
  filename2 = '/projects/0/samoc/jan/Andree/LEC_bin_5_'//year  ! LEC.f90 output file
  !write (*,*) filename1
  !write (*,*) filename2
  write (*,*) filename1, filename2
  open (1,file=trim(filename1),access='direct',form='unformatted',recl=imt*jmt,      &
        status='unknown')
  open (2,file=trim(filename2),access='direct',form='unformatted',recl=imt*jmt,      &
        status='unknown')

  ! loading VVEL, WVEL, TTT_VVEL
  allocate( VVEL(imt,jmt,km), TTT_VVEL(imt,jmt,km) )
  
  call load_3D_field(1,nrec_VVEL,VVEL)
  !call load_3D_field(1,nrec_WVEL,WVEL)
  call uu2tt_3D(DZT,DZU,TAREA,UAREA,VVEL,TTT_VVEL)
  !=============================================================================
  !  1. cPKm/cPKe
  !=============================================================================

  if ( opt1==1 ) then
    allocate( cPKm(imt,jmt,km),   cPKe(imt,jmt,km),                            &
              cPKm_vint(imt,jmt), cPKe_vint(imt,jmt),                          &
              cPKm_zint(jmt),     cPKe_zint(jmt),                              &
              cPKm_rsum(jmt),     cPKe_rsum(jmt) )
  
    ! load fields
    call load_3D_field(2,nrec_cPKm,cPKm)
    call load_3D_field(2,nrec_cPKe,cPKe)
  
    !vertical integral
    call vert_int(cPKm,FBC,cPKm_vint)
    call vert_int(cPKe,FBC,cPKe_vint)
  
    ! zonal-depth integral
    call zonal_int(DXT,cPKm_vint,cPKm_zint)
    call zonal_int(DXT,cPKe_vint,cPKe_zint)
  
    ! write into interal fields cPKm_vint
    write ( 3,rec=y-start_year+1) cPKm_vint(:,:)
    write ( 4,rec=y-start_year+1) cPKe_vint(:,:)
    write(*,*) 'option 1'
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
    call load_3D_field(2,nrec_rPm, rPm )
    call load_3D_field(2,nrec_rPe, rPe )
    call load_3D_field(2,nrec_rKm, rKm )
    call load_3D_field(2,nrec_rKe, rKE )

    call vol_int_full(1,1,1,imt,jmt,km,rPm,TAREA,dz,DZT,rPm_int)
    call vol_int_full(1,1,1,imt,jmt,km,rPe,TAREA,dz,DZT,rPe_int)
    call vol_int_part(1,1,1,imt,jmt,km,rKm,TAREA,dz,DZT,rKm_int)
    call vol_int_part(1,1,1,imt,jmt,km,rKe,TAREA,dz,DZT,rKe_int)
    write (*,*) rPm_int, rPe_int, rKm_int, rKe_int
  
    ! calculate vertically integrated, meridional advection
    call mer_advection(DXT,FBC,TTT_VVEL,rPm,vrPm_vint)
    call mer_advection(DXT,FBC,TTT_VVEL,rPe,vrPe_vint)
    call mer_advection(DXT,DZT,TTT_VVEL,rKm,vrKm_vint)
    call mer_advection(DXT,DZT,TTT_VVEL,rKe,vrKe_vint)

    ! calculate zonal-depth integral
    call zonal_int(DXT,vrPm_vint,vrPm_zint)
    call zonal_int(DXT,vrPe_vint,vrPe_zint)
    call zonal_int(DXT,vrKm_vint,vrKm_zint)
    call zonal_int(DXT,vrKe_vint,vrKe_zint)
  
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
                vrPm_rsum, vrPe_rsum, vrKm_rsum, vrKe_rsum,                    &
                VVEL, TTT_VVEL )

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

  call zonal_int(DXT,cPKm_avg,cPKm_zint)
  call zonal_int(DXT,cPKe_avg,cPKe_zint)
  write (11,rec=nt+1) cPKm_zint(:)
  write (12,rec=nt+1) cPKe_zint(:)

  deallocate( cPKm_yrly, cPKe_yrly, cPKm_avg, cPKe_avg,                        &
              cPKm_zint, cPKe_zint, cPKm_rsum, cPKe_rsum)

  close( 3)
  close( 4)
  close(11)
  close(12)
!  close(25)
!  close(26)
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

  call zonal_int(DXT,vrPm_avg,vrPm_zint)
  call zonal_int(DXT,vrPe_avg,vrPe_zint)
  call zonal_int(DXT,vrKm_avg,vrKm_zint)
  call zonal_int(DXT,vrKe_avg,vrKe_zint)
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
!  close(27)
!  close(28)
!  close(29)
!  close(30)
endif


!===============================================================================
!  SUBROUTINES
!===============================================================================

contains


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!subroutine rsum(imt,jmt,km,DXT,DYT,DZT,FIELD,rsum)
!implicit none
!
!  calculates running sum
!

!integer, intent(in)  :: imt,jmt,km,opt
!integer                                     :: j
!real,    dimension(imt,jmt,km), intent(in)  :: DZT
!if ( opt==1) then
!  real,  dimension(imt,jmt,km), intent(in)  :: FIELD
!elseif ( opt==2 ) then
!  real,  dimension(imt,jmt),    intent(in)  :: FIELD
!endif
!real,    dimension(jmt),        intent(out) :: rsum
!real,    dimension(imt,jmt),    intent(in)  :: DXT,DYT
!
!if ( opt==1 ) then
!  rsum(:) = sum(sum(FIELD(:,:,:)*DZT(:,:,:),3)*DXT(:,:)*DYT(:,:),1)
!  do j=2,jmt
!    rsum(j) = rsum(j) + rsum(j-1)
!  enddo
!elseif ( opt==2 ) then
!  rsum(:) = sum(FIELD(:,:)*DXT(:,:)*DYT(:,:),1)
!  do j=2,jmt
!    rsum(j) = rsum(j) + rsum(j-1)
!  enddo
!endif

!end subroutine rsum

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
end program
