program int_test
implicit none

!===============================================================================
!
!  compares surface integrals un-/weighted byt pbc-depth
!  
!===============================================================================

!===============================================================================
!  variables
!===============================================================================


character*120 ::                                                               &
  input_folder,in_depths,geometry1_file,geometry2_file,int_test_file,tavg_file,&
  LEC_file,projects_folder

integer                                         ::                             &
  imt,jmt,km,rec_length,i,j,k,                                                 &
  nrec_rPm, nrec_rPe, nrec_rKm, nrec_rKe,                                      &
  nrec_gPm, nrec_gPe, nrec_gKm, nrec_gKe,                                      &
  nrec_gPmh, nrec_gPms, nrec_gPeh, nrec_gPes,                                  &
  nrec_cPem, nrec_cKem, nrec_cPKm, nrec_cPKe,                                  &
  nrec_TEMP

real                                            ::                             &
  rPm_int,rPe_int,rKm_int,rKe_int,cPem_int,cKem_int,cPKm_int,cPKe_int,TEMP_int,&
  gPm_int,gPe_int,gKm_int,gKe_int,gPmh_int,gPms_int,gPeh_int,gPes_int

real,             dimension(:),     allocatable ::                             &
  dz,tdepth,area,p_z,vol,                                                      &
  rPm_avg, rPe_avg, rKm_avg, rKe_avg, cPem_avg, cKem_avg, cPKm_avg, cPKe_avg,  &
  rPm_wavg,rPe_wavg,rKm_wavg,rKe_wavg,cPem_wavg,cKem_wavg,cPKm_wavg,cPKe_wavg
real,             dimension(:,:),   allocatable ::                             &
  DXT, DYT, TAREA, DXU, DYU, UAREA, geometry2,                                 &
  gPm,gPe,gKm,gKe,gPmh,gPms,gPeh,gPes
real,             dimension(:,:,:), allocatable ::                             &
  rPm,rPe,rKm,rKe,cPem,cKem,cPKm,cPKe,TEMP,&
  DZT,DZU
real, parameter                                 ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.,rho0 = 4.1/3.996*1000

imt            = 3600
jmt            = 2400
km             =   42

input_folder   = '/home/dijkbio/andre/LEC/input/'
projects_folder = '/projects/0/samoc/jan/Andree/'
in_depths      = trim(input_folder)//'in_depths.42.dat'
geometry1_file = trim(input_folder)//'geometry1'
geometry2_file = trim(input_folder)//'geometry2'

int_test_file  = '/home/dijkbio/andre/LEC/results/int_test/int_test'

tavg_file      = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.51.year.301'
LEC_file       = trim(projects_folder)//'LEC_bin_5_301'


write (*,*) ''
write (*,*) '--- INTEGRAL TESTING ---'
write (*,*) ''

!===============================================================================
!  LOAD GEOMETRY
!===============================================================================

allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km),         &
          DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt), DZU(imt,jmt,km) )
open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,     &
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

allocate( geometry2(6,km),dz(km),tdepth(km),area(km),p_z(km),vol(km) )
open(1,file=geometry2_file,access='direct',form='unformatted',recl=6,     &
       status='old')
do k=1,km
  read(1,rec=k) geometry2(:,k) ! all real: k, dz[m], tdepth[m], area[m^2], p[bar]
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
p_z(:)    = geometry2(5,:)
vol(:)    = geometry2(6,:)
close(1)

!===============================================================================
!  LOAD LEC FIELDS
!===============================================================================

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

allocate(                                                                      &
rPm(imt,jmt,km),  rPe(imt,jmt,km),  rKm(imt,jmt,km),  rKe(imt,jmt,km),         &
gPm(imt,jmt),     gPe(imt,jmt),     gKm(imt,jmt),     gKe(imt,jmt),            &
gPmh(imt,jmt),    gPms(imt,jmt),    gPeh(imt,jmt),    gPes(imt,jmt),           &
cPem(imt,jmt,km), cKem(imt,jmt,km), cPKm(imt,jmt,km), cPKe(imt,jmt,km),        &
TEMP(imt,jmt,km) )

!open file
open (1,file=LEC_file,access='direct',form='unformatted',recl=imt*jmt,  &
        status='unknown')
write (*,*) 'file: ', LEC_file

call load_3D_field(1,nrec_rPm ,rPm )
call load_3D_field(1,nrec_rPe ,rPe )
call load_3D_field(1,nrec_rKm ,rKm )
call load_3D_field(1,nrec_rKe ,rKe )

!call load_2D_field(   1,nrec_gPm ,gPm )
!call load_2D_field(   1,nrec_gPe ,gPe )
!call load_2D_field(   1,nrec_gKm ,gKm )
!call load_2D_field(   1,nrec_gKe ,gKe )
!call load_2D_field(   1,nrec_gPmh,gPmh)
!call load_2D_field(   1,nrec_gPms,gPms)
!call load_2D_field(   1,nrec_gPeh,gPeh)
!call load_2D_field(   1,nrec_gPes,gPes)

call load_3D_field(1,nrec_cPem,cPem)
call load_3D_field(1,nrec_cKem,cKem)
call load_3D_field(1,nrec_cPKm,cPKm)
call load_3D_field(1,nrec_cPKe,cPKe)

!call load_3D_field(1,nrec_TEMP,TEMP)

close(1)

!===============================================================================
!  TESTING
!===============================================================================

write(*,*) '  TESTING starts'

allocate( rPm_avg(km),   rPe_avg(km),   rKm_avg(km),   rKe_avg(km), & 
          cPKm_avg(km),  cPKe_avg(km),  cPem_avg(km),  cKem_avg(km) )
allocate( rPm_wavg(km),  rPe_wavg(km),  rKm_wavg(km),  rKe_wavg(km), & 
          cPKm_wavg(km), cPKe_wavg(km), cPem_wavg(km), cKem_wavg(km) )

call area_avg( rPm,DZT,TAREA,area, rPm_avg)
call area_avg( rPe,DZT,TAREA,area, rPe_avg)
call area_avg( rKm,DZT,TAREA,area, rKm_avg)
call area_avg( rKe,DZT,TAREA,area, rKe_avg)
call area_avg(cPem,DZT,TAREA,area,cPem_avg)
call area_avg(cKem,DZT,TAREA,area,cKem_avg)
call area_avg(cPKm,DZT,TAREA,area,cPKm_avg)
call area_avg(cPKe,DZT,TAREA,area,cPKe_avg)

call area_avg_weighted( rPm,DZT,TAREA,vol, rPm_wavg)
call area_avg_weighted( rPe,DZT,TAREA,vol, rPe_wavg)
call area_avg_weighted( rKm,DZT,TAREA,vol, rKm_wavg)
call area_avg_weighted( rKe,DZT,TAREA,vol, rKe_wavg)
call area_avg_weighted(cPem,DZT,TAREA,vol,cPem_wavg)
call area_avg_weighted(cKem,DZT,TAREA,vol,cKem_wavg)
call area_avg_weighted(cPKm,DZT,TAREA,vol,cPKm_wavg)
call area_avg_weighted(cPKe,DZT,TAREA,vol,cPKe_wavg)

do k=1,km
  write(*,*) rPm_avg(k)/rPm_wavg(k),&
             rPe_avg(k)/rPe_wavg(k),&
             rKm_avg(k)/rKm_wavg(k),&
             rKe_avg(k)/rKe_wavg(k)
enddo

write(*,*) ''
do k=1,km
  write(*,*) cPem_avg(k)/cPem_wavg(k),&
             cKem_avg(k)/cKem_wavg(k),&
             cPKm_avg(k)/cPKm_wavg(k),&
             cPKe_avg(k)/cPKe_wavg(k)
enddo

!write(*,*) 'avg*vol (inconsistent)'
!write(*,*) sum(cPem_avg*vol),sum(cKem_avg*vol),&
!           sum(cPKm_avg*vol),sum(cPKe_avg*vol)
write(*,*) 'avg*area*dz (consistent, unweighted)'
write(*,*) sum(cPem_avg*area*dz),sum(cKem_avg*area*dz),&
           sum(cPKm_avg*area*dz),sum(cPKe_avg*area*dz)
write(*,*) 'wavg*vol (consistent, weighted)'
write(*,*) sum(cPem_wavg*vol),sum(cKem_wavg*vol),&
           sum(cPKm_wavg*vol),sum(cPKe_wavg*vol)
!write(*,*) 'wavg*area*dz (inconsistent)'
!write(*,*) sum(cPem_wavg*area*dz),sum(cKem_wavg*area*dz),&
!           sum(cPKm_wavg*area*dz),sum(cPKe_wavg*area*dz)
write(*,*) 'vol_int'
call vol_int(1,1,1,imt,jmt,km,cPem,TAREA,DZT,cPem_int)
call vol_int(1,1,1,imt,jmt,km,cKem,TAREA,DZT,cKem_int)
call vol_int(1,1,1,imt,jmt,km,cPKm,TAREA,DZT,cPKm_int)
call vol_int(1,1,1,imt,jmt,km,cPKe,TAREA,DZT,cPKe_int)
write(*,*) cPem_int,cKem_int,cPKm_int,cPKe_int
write(*,*) 'vol_int2'
call vol_int2(1,1,1,imt,jmt,km,cPem,TAREA,DZT,cPem_int)
call vol_int2(1,1,1,imt,jmt,km,cKem,TAREA,DZT,cKem_int)
call vol_int2(1,1,1,imt,jmt,km,cPKm,TAREA,DZT,cPKm_int)
call vol_int2(1,1,1,imt,jmt,km,cPKe,TAREA,DZT,cPKe_int)
write(*,*) cPem_int,cKem_int,cPKm_int,cPKe_int


!===============================================================================
!  OUTPUT
!===============================================================================

! create readme
! >>>


open(1,file=int_test_file,access='direct',form='unformatted',recl=16,          &
       status='unknown')
do k=1,km
  write(1,rec=k) rPm_avg(k),   rPe_avg(k),   rKm_avg(k),   rKe_avg(k),         &
                 cPKm_avg(k),  cPKe_avg(k),  cPem_avg(k),  cKem_avg(k),        &
                 rPm_wavg(k),  rPe_wavg(k),  rKm_wavg(k),  rKe_wavg(k),        & 
                 cPKm_wavg(k), cPKe_wavg(k), cPem_wavg(k), cKem_wavg(k)
enddo
close(1)

!===============================================================================
!===============================================================================
end program int_test
