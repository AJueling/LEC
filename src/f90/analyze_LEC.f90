program analyze_LEC
implicit none

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  this script calculates the integrals of the LEC output file created with
!  LEC.f90
!  
!  for curvilinear grid:
!  longitude (i) 0=250E, 1100=0E, 1500=40E, 1700=60E
!  latitude  (j) 0=78.5S, 518=55S, 677=45S, 866=30S, 1181=0S
!  depth     (k) 1=5m, 9=96.9m, 10=112m, 14=198m, 17=318m, 35=4125, 42=5875m
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!===============================================================================
!  variables
!===============================================================================

! user input 
integer ::                                                                     &
  imt, jmt, km, ntavg, k,                                                      &
  nrec_rPm, nrec_rPe, nrec_rKm, nrec_rKe,                                      &
  nrec_gPm, nrec_gPe, nrec_gKm, nrec_gKe,                                      &
  nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                                  &
  nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                                  &
  nrec_TEMP, nrec_PU_x, nrec_PV_y

character*120 :: LEC_folder,geometry1_file,geometry2_file,LEC_file,            &
  output_file, output_folder

character*3   :: year

real,    dimension(:),     allocatable :: dz,area,tdepth

real,    dimension(:,:),   allocatable ::                                      &
  DXT, DYT, TAREA, DXU, DYU, UAREA, geometry2,                                 &
  gPm, gPe, gKm, gKe, gPmh, gPeh, gPms, gPes

real,    dimension(:,:,:), allocatable ::                                      &
  rPm, rPe, rKm, rKe, cPem, cKem, cPKm, cPKe, TEMP, DZT, DZU, PU_x, PV_y, FBC

imt = 3600
jmt = 2400
km  =   42

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
nrec_TEMP =  345
nrec_PU_x =  387
nrec_PV_y =  429

write (*,*) ''
write (*,*) '--- ANALYZING ENERGY BUDGET ---'
write (*,*) ''

!===============================================================================
!  INPUT
!===============================================================================

read  (*,'(a120)') LEC_file
read  (*,'(a3)')   year
read  (*,*)        ntavg

output_folder = '/home/dijkbio/andre/LEC/results/analyze_LEC/'
if      ( ntavg==1  ) then
  output_file = trim(output_folder)//'analysis_LEC_1_'//year
else if ( ntavg==5  ) then
  output_file = trim(output_folder)//'analysis_LEC_5_'//year
else if ( ntavg==11 ) then
  output_file = trim(output_folder)//'analysis_LEC_11_'//year
endif
write(*,*) ntavg, year
write(*,*) output_file
write(*,*) output_folder
!===============================================================================
!  GEOMETRY 
!===============================================================================

LEC_folder      = '/home/dijkbio/andre/LEC/'
geometry1_file  = trim(LEC_folder)//'input/geometry1'
geometry2_file  = trim(LEC_folder)//'input/geometry2'

! read 2D geometry fields
allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km) )
open(2,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,    &
       status='old')
read(2,rec=1) DXT ! [m]
read(2,rec=2) DXU ! [m]
read(2,rec=3) DYT ! [m]
read(2,rec=4) DYU ! [m]
read(2,rec=5) TAREA ! [m^2]
read(2,rec=6) UAREA ! [m^2]
do k=1,km
  read(2,rec=6+   k) DZT(:,:,k) ! [m]
  read(2,rec=6+km+k) DZU(:,:,k) ! [m]
enddo
close(2)

! read 1D geometry fields
! all real: k,dz[m],tdepth[m],area[m^2],p[bar]
allocate( geometry2(6,km),dz(km),tdepth(km),area(km) )
open(2,file=geometry2_file,access='direct',form='unformatted',recl=6,     &
       status='old')
do k=1,km
  read(2,rec=k) geometry2(:,k)
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
close(2)

!===============================================================================
!  LOAD LEC FIELDS
!===============================================================================

allocate(                                                                      &
rPm(imt,jmt,km),  rPe(imt,jmt,km),  rKm(imt,jmt,km),  rKe(imt,jmt,km),         &
gPm(imt,jmt),     gPe(imt,jmt),     gKm(imt,jmt),     gKe(imt,jmt),            &
gPmh(imt,jmt),    gPms(imt,jmt),    gPeh(imt,jmt),    gPes(imt,jmt),           &
cPem(imt,jmt,km), cKem(imt,jmt,km), cPKm(imt,jmt,km), cPKe(imt,jmt,km),        & 
TEMP(imt,jmt,km), PU_x(imt,jmt,km), PV_y(imt,jmt,km) ) 


!open file
open (1,file=LEC_file,access='direct',form='unformatted',recl=imt*jmt,  &
        status='unknown')
write (*,*) 'file: ', LEC_file

call load_3D_field(1,nrec_rPm ,rPm )
call load_3D_field(1,nrec_rPe ,rPe )
call load_3D_field(1,nrec_rKm ,rKm )
call load_3D_field(1,nrec_rKe ,rKe )

call load_2D_field(1,nrec_gPm ,gPm )
call load_2D_field(1,nrec_gPe ,gPe )
call load_2D_field(1,nrec_gKm ,gKm )
call load_2D_field(1,nrec_gKe ,gKe )
call load_2D_field(1,nrec_gPmh,gPmh)
call load_2D_field(1,nrec_gPms,gPms)
call load_2D_field(1,nrec_gPeh,gPeh)
call load_2D_field(1,nrec_gPes,gPes)

call load_3D_field(1,nrec_cPem,cPem)
call load_3D_field(1,nrec_cKem,cKem)
call load_3D_field(1,nrec_cPKm,cPKm)
call load_3D_field(1,nrec_cPKe,cPKe)
call load_3D_field(1,nrec_TEMP,TEMP)
call load_3D_field(1,nrec_PU_x,PU_x)
call load_3D_field(1,nrec_PV_y,PV_y)

close(1)

!===============================================================================
!  CALCULATIONS
!===============================================================================
write (*,*) ''
write (*,*) 'CALCULATIONS'

! Southern Ocean integrals [90S,30S]
call regional_integrals(1,1,imt,866,'SO30')
! Southern Ocean integrals [90S,30S]x[60W,30E]
call regional_integrals(500,1,1400,866,'SA30')
! Southern Ocean/Weddell Gyre integrals [90S,30S]x[50W,60E]
call regional_integrals(600,1,1700,866,'WG30')
! Southern Ocean integrals [90S,45S]
call regional_integrals(1,1,imt,677,'SO45')
! Southern Ocean integrals [90S,55S]
call regional_integrals(1,1,imt,518,'SO55')
! Weddell Gyre to Kerguelen Plateau  integrals [90S,55S]x[35W,80E]
call regional_integrals(750,1,1900,600,'WGKP')
! North of Southern Ocean [30S,90N]
call regional_integrals(1,867,imt,jmt,'NO30P')
! North of Southern Ocean [30S,90N]
call regional_integrals(800,1,1300,518,'SOMSR')
write (*,*) '>>> done writing output'


!===============================================================================
!===============================================================================
contains

!===============================================================================
subroutine regional_integrals(imin,jmin,imax,jmax,area_name)

character*4         :: area_name
integer, intent(in) :: imin,jmin,imax,jmax
integer, parameter  :: imt=3600,jmt=2400,km=42
real                ::                                                         & 
  gPm_int ,gPe_int ,gKm_int ,gKe_int ,                                         &
  gPmh_int,gPeh_int,gPms_int,gPes_int,                                         &
  rPm_int ,rPe_int ,rKm_int ,rKe_int ,                                         &
  cPem_int,cKem_int,cPKm_int,cPKe_int,                                         &
  rPm_int_top ,rPe_int_top ,rKm_int_top ,rKe_int_top ,                         &
  cPem_int_top,cKem_int_top,cPKm_int_top,cPKe_int_top,                         &
  rPm_int_bot ,rPe_int_bot ,rKm_int_bot ,rKe_int_bot ,                         &
  cPem_int_bot,cKem_int_bot,cPKm_int_bot,cPKe_int_bot,                         &
  ohc_int,ohc_int_top,ohc_int_bot,                                             &
  PU_x_int,PU_x_int_top,PU_x_int_bot,                                          &
  PV_y_int,PV_y_int_top,PV_y_int_bot

! generation
call surf_int_2D(imin,jmin,imax,jmax,gPm,TAREA,DZT(:,:,1),gPm_int)
call surf_int_2D(imin,jmin,imax,jmax,gPe,TAREA,DZT(:,:,1),gPe_int)
call surf_int_2D(imin,jmin,imax,jmax,gKm,TAREA,DZT(:,:,1),gKm_int)
call surf_int_2D(imin,jmin,imax,jmax,gKe,TAREA,DZT(:,:,1),gKe_int)

call surf_int_2D(imin,jmin,imax,jmax,gPmh,TAREA,DZT(:,:,1),gPmh_int)
call surf_int_2D(imin,jmin,imax,jmax,gPms,TAREA,DZT(:,:,1),gPms_int)
call surf_int_2D(imin,jmin,imax,jmax,gPeh,TAREA,DZT(:,:,1),gPeh_int)
call surf_int_2D(imin,jmin,imax,jmax,gPes,TAREA,DZT(:,:,1),gPes_int)

! full depth
call r_c_integrals(imin,jmin,1,imax,jmax,km,                                   &
                   rPm_int,rPe_int,rKm_int,rKe_int,                            &
                   cPem_int,cKem_int,cPKm_int,cPKe_int,                        &
                   ohc_int, PU_x_int, PV_y_int)

! top
call r_c_integrals(imin,jmin,1,imt,jmt,17,                                     &
                   rPm_int_top,rPe_int_top,rKm_int_top,rKe_int_top,            &
                   cPem_int_top,cKem_int_top,cPKm_int_top,cPKe_int_top,        &
                   ohc_int_top, PU_x_int_top, PV_y_int_top)

! bottom
call r_c_integrals(imin,jmin,35,imt,jmt,km,                                    &
                   rPm_int_bot,rPe_int_bot,rKm_int_bot,rKe_int_bot,            &
                   cPem_int_bot,cKem_int_bot,cPKm_int_bot,cPKe_int_bot,        &
                   ohc_int_bot, PU_x_int_bot, PV_y_int_bot)


! output
open(3,file=trim(output_file)//'_'//area_name//'.out',form='formatted')

200 FORMAT (40(A,","),A)
201 FORMAT (40(E15.7,","),E15.7)
write (3,200)                                                                  &
  'gPm' ,'gPe' ,'gKm' ,'gKe' ,                                                 &
  'gPmh','gPeh','gPms','gPes',                                                 &
  'rPm' ,'rPe' ,'rKm' ,'rKe' ,                                                 &
  'cPem','cKem','cPKm','cPKe',                                                 &
  'rPm_top' ,'rPe_top' ,'rKm_top' ,'rKe_top' ,                                 &
  'cPem_top','cKem_top','cPKm_top','cPKe_top',                                 &
  'rPm_bot' ,'rPe_bot' ,'rKm_bot' ,'rKe_bot' ,                                 &
  'cPem_bot','cKem_bot','cPKm_bot','cPKe_bot',                                 &
  'temp_int','temp_int_top','temp_int_bot',                                    &
  'PU_x_int','PU_x_int_top','PU_x_int_bot',                                    &
  'PV_y_int','PV_y_int_top','PV_y_int_bot'
write (3,201)                                                                  &
  gPm_int ,gPe_int ,gKm_int ,gKe_int ,                                         &
  gPmh_int,gPeh_int,gPms_int,gPes_int,                                         &
  rPm_int ,rPe_int ,rKm_int ,rKe_int ,                                         &
  cPem_int,cKem_int,cPKm_int,cPKe_int,                                         &
  rPm_int_top ,rPe_int_top ,rKm_int_top ,rKe_int_top ,                         &
  cPem_int_top,cKem_int_top,cPKm_int_top,cPKe_int_top,                         &
  rPm_int_bot ,rPe_int_bot ,rKm_int_bot ,rKe_int_bot ,                         &
  cPem_int_bot,cKem_int_bot,cPKm_int_bot,cPKe_int_bot,                         &
  ohc_int,ohc_int_top,ohc_int_bot,                                             &
  PU_x_int,PU_x_int_top,PU_x_int_bot,                                          &
  PV_y_int,PV_y_int_top,PV_y_int_bot
close(3)
end subroutine regional_integrals
!===============================================================================

!===============================================================================
subroutine r_c_integrals(imin,jmin,kmin,imax,jmax,kmax,                        &
                         rPm_int,rPe_int,rKm_int,rKe_int,                      &
                         cPem_int,cKem_int,cPKm_int,cPKe_int,                  &
                         ohc_int, PU_x_int, PV_y_int )

integer, intent(in) :: imin,jmin,kmin,imax,jmax,kmax
real, intent(out)   :: rPm_int,rPe_int,rKm_int,rKe_int,                        &
                       cPem_int,cKem_int,cPKm_int,cPKe_int,                    &
                       ohc_int, PU_x_int, PV_y_int

! reservoirs
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,rPm,TAREA,dz,DZT,rPm_int)
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,rPe,TAREA,dz,DZT,rPe_int)
call vol_int_part(imin,jmin,kmin,imax,jmax,kmax,rKm,TAREA,dz,DZT,rKm_int)
call vol_int_part(imin,jmin,kmin,imax,jmax,kmax,rKe,TAREA,dz,DZT,rKe_int)

! exchange terms
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,cPem,TAREA,dz,DZT,cPem_int)
call vol_int_part(imin,jmin,kmin,imax,jmax,kmax,cKem,TAREA,dz,DZT,cKem_int)
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,cPKm,TAREA,dz,DZT,cPKm_int)
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,cPKe,TAREA,dz,DZT,cPKe_int)

! OHC
call vol_int_part(imin,jmin,kmin,imax,jmax,kmax,TEMP,TAREA,dz,DZT,ohc_int)

! pressure work terms
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,PU_x,TAREA,dz,DZT,PU_x_int)
call vol_int_full(imin,jmin,kmin,imax,jmax,kmax,PV_y,TAREA,dz,DZT,PV_y_int)

end subroutine r_c_integrals
!===============================================================================

!===============================================================================
!===============================================================================
end program
