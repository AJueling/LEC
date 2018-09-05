program pbc_test
implicit none

!===============================================================================
!
!  testing calculation on cPKm/e and rPm at vertical center of full bottom cell, 
!  not at center of partial bottom cell as before
!
!  terminology: fbc: full bottom cell, pbs: partial bottom cell
!  
!===============================================================================




!===============================================================================
!  variables
!===============================================================================


character*120 ::                                                               &
  input_folder,in_depths,geometry1_file,geometry2_file,tavg51_file,            &
  LEC_file,projects_folder,kmt_file,pbc_test_file

integer                                         ::                             &
  imt,jmt,km,rec_length,i,j,k,opt,                                             &
  nrec_PD,nrec_WVEL,nrec_PDW,nrec_TEMP,nrec_SALT,nrec_WPD

integer,          dimension(:,:),   allocatable ::                             &
  kmT

real                                            ::                             &
  cPKm_o_pint,cPKe_o_pint,cPKt_o_pint,cPKm_o_fint,cPKe_o_fint,cPKt_o_fint,     &
  cPKm_a_pint,cPKe_a_pint,cPKt_a_pint,cPKm_a_fint,cPKe_a_fint,cPKt_a_fint,     &
  cPKm_b_pint,cPKe_b_pint,cPKt_b_pint,cPKm_b_fint,cPKe_b_fint,cPKt_b_fint,     &
  pressure,pbar

real,             dimension(:),     allocatable ::                             &
  dz,tdepth,area,p_z,vol

real,             dimension(:,:),   allocatable ::                             &
  DXT,DYT,TAREA,DXU,DYU,UAREA,geometry2,fac,dz_k,dz_k1

real,             dimension(:,:,:), allocatable ::                             &
  DZT,DZU,PD,WVEL,PDW,WPD,TTT_WVEL,PDW_new,cPKm,cPKe,cPKt,TEMP,SALT,           &
  PD_a,PD_b,TEMP_b,SALT_b

real, parameter                                 ::                             &
  c0=0.,p125=0.125,p25=0.25,p5=0.5,c1=1.,rho0 = 4.1/3.996*1000,g=9.806

double precision                                ::                             &
  PD_dble

!read (*,*) opt
opt            = 1 ! density option: 1 - PD, 2 - RHO
!opt            = 2 ! density option: 1 - PD, 2 - RHO
!write(*,*) opt

imt            = 3600
jmt            = 2400
km             =   42

input_folder   = '/home/dijkbio/andre/LEC/input/'
projects_folder = '/projects/0/samoc/jan/Andree/'
in_depths      = trim(input_folder)//'in_depths.42.dat'
kmt_file       = trim(input_folder)//'kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'
geometry1_file = trim(input_folder)//'geometry1'
geometry2_file = trim(input_folder)//'geometry2'

pbc_test_file  = '/home/dijkbio/andre/LEC/results/pbc_test/pbc_test'


tavg51_file    = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.year.301'
!tavg51_file    = trim(projects_folder)//'t.t0.1_42l_nccs01.tavg.51.year.301'
LEC_file       = trim(projects_folder)//'LEC_bin_5_301'

write(*,*) ''
write(*,*) '  Testing effects of pbc '
if (opt==1) then
  write(*,*) '  using PD, PDW, and WPD'
elseif (opt==2) then
  write(*,*) '  using RHO, RHOW, and WRHO'
endif
write(*,*) ''

!===============================================================================
!  LOAD GEOMETRY
!===============================================================================

allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km) )
open(1,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,     &
       status='old')
read(1,rec=1) DXT ! [m]
read(1,rec=3) DYT ! [m]
read(1,rec=5) TAREA ! [m^2]
do k=1,km
  read(1,rec=6+   k) DZT(:,:,k) ! [m]
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

allocate( kmT(imt,jmt) )
inquire (iolength=rec_length) kmT
open(1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,   &
           status='old')
  read(1,rec=1) kmT
close(1)
!===============================================================================
!  LOAD LEC FIELDS
!===============================================================================

if (opt==1) then
  nrec_PD   =  995 ! PD
  nrec_PDW  = 1457 ! PDW
  nrec_WPD  = 1667 ! WPD
elseif (opt==2) then
  nrec_PD   =  306 ! RHO
  nrec_PDW  = 1289 ! RHOW
  nrec_WPD  = 1415 ! WRHO
endif

nrec_WVEL =  405
nrec_SALT =  264
nrec_TEMP =  222

open (1,file=tavg51_file,access='direct',form='unformatted',recl=imt*jmt,  &
        status='unknown')
write (*,*) 'file: ', tavg51_file

allocate( PD(imt,jmt,km), WVEL(imt,jmt,km), PDW(imt,jmt,km), WPD(imt,jmt,km),&
          TEMP(imt,jmt,km), SALT(imt,jmt,km) )
call load_3D_field(1,nrec_PD  ,PD  )
call load_3D_field(1,nrec_WVEL,WVEL)
call load_3D_field(1,nrec_PDW ,PDW )
call load_3D_field(1,nrec_WPD ,WPD )
call load_3D_field(1,nrec_TEMP,TEMP)
call load_3D_field(1,nrec_SALT,SALT)

close(1)

!===============================================================================
!   testing pressure level of PD
!===============================================================================

pbar = pressure(dz(1)/2.0)
!write(*,*) dz(1)/2.0,'m ;',pressure(dz(1)/2.0), 'bar'
!call state(dble(SALT(1000,1000,10)),dble(TEMP(1000,1000,10)),dble(pbar),PD_dble)
!write(*,*) 'z=dz(1)/2',PD(1000,1000,10),PD_dble/1.0E03
!call state(dble(SALT(1000,1000,10)),dble(TEMP(1000,1000,10)),dble(0),PD_dble)
!write(*,*) 'z=0      ',PD(1000,1000,10),PD_dble/1.0E03
! >>> cannot exactly reproduce PD error in O(10^-6 kg/m^3)

!===============================================================================
!   OLD
!   calculate original cPK terms for comparison
!===============================================================================

allocate( TTT_WVEL(imt,jmt,km) )
call wtt2ttt(WVEL,DZT,TTT_WVEL) ! @ bottom cell: 0.5 * WVEL(top of bottom cell)

allocate( cPKm(imt,jmt,km), cPKe(imt,jmt,km), cPKt(imt,jmt,km)  )
do k = 1, km
  cPKt(:,:,k) = -g * PDW(:,:,k)*1.0E01
  cPKm(:,:,k) = -g * PD(:,:,k) * TTT_WVEL(:,:,k) * 1.0E01
  cPKe(:,:,k) = cPKt(:,:,k) - cPKm(:,:,k)
enddo !k

write(*,*) 'o - original, a - density gradient, b - temp/salinity gradient'
write(*,*) ''
write(*,*) '               partial integral                full integral'
call vol_int(1,1,1,imt,jmt,km,cPKm,TAREA,DZT,cPKm_o_pint)
call vol_int(1,1,1,imt,jmt,km,cPKe,TAREA,DZT,cPKe_o_pint)
call vol_int(1,1,1,imt,jmt,km,cPKt,TAREA,DZT,cPKt_o_pint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKm,dz,TAREA,DZT,cPKm_o_fint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKe,dz,TAREA,DZT,cPKe_o_fint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKt,dz,TAREA,DZT,cPKt_o_fint)
write(*,*) 'original'
write(*,*) 'cPKm_o =    ',&
            cPKm_o_pint, 'W', cPKm_o_pint/-4.9E11,&
            cPKm_o_fint, 'W', cPKm_o_fint/-4.9E11
write(*,*) 'cPKe_o =    ', &
            cPKe_o_pint, 'W', cPKe_o_pint/ 7.3E11,&
            cPKe_o_fint, 'W', cPKe_o_fint/ 7.3E11
write(*,*) 'cPKt_o =    ',&
            cPKt_o_pint, 'W', cPKt_o_pint/ 2.4E11,&
            cPKt_o_fint, 'W', cPKt_o_fint/ 2.4E11

!===============================================================================
!   NEW
!===============================================================================

allocate( PD_a(imt,jmt,km),   PD_b(imt,jmt,km), PDW_new(imt,jmt,km), &
          TEMP_b(imt,jmt,km), SALT_b(imt,jmt,km) )
do k=1,km
  PD_a(:,:,k)    = PD(:,:,k)
  PD_b(:,:,k)    = PD(:,:,k)
  TEMP_b(:,:,k)  = TEMP(:,:,k)
  SALT_b(:,:,k)  = SALT(:,:,k)
  PDW_new(:,:,k) = PDW(:,:,k)
enddo

allocate( fac(imt,jmt), dz_k(imt,jmt), dz_k1(imt,jmt) )



do k = 2,km
  !=============================================================================
  ! 1.  cPKm  
  !     calculate new density at center of fbc
  ! (a) extrapolation of density based on vert. density gradient
  ! (b) extrapolation via vert. TEMP/SALT gradient, then eof
  !=============================================================================
  dz_k(:,:)  = dz(k)
  dz_k1(:,:) = dz(k-1)
  if (opt==1) then
    pbar = pressure(dz(1)/2.0) ! for PD case
  elseif (opt==2) then  
    pbar = pressure(tdepth(k)) ! for RHO case
  endif

  where(kmT==k)
    fac(:,:) = ( dz_k1(:,:)+dz_k(:,:) ) / ( dz_k1(:,:)+DZT(:,:,k) )
  ! (a)
    PD_a(:,:,k)   =   PD(:,:,k-1) + (   PD(:,:,k)-  PD(:,:,k-1) ) * fac(:,:)
  ! (b)
    TEMP_b(:,:,k) = TEMP(:,:,k-1) + ( TEMP(:,:,k)-TEMP(:,:,k-1) ) * fac(:,:)
    SALT_b(:,:,k) = SALT(:,:,k-1) + ( SALT(:,:,k)-SALT(:,:,k-1) ) * fac(:,:)
  endwhere
  !write(*,*) k, sum(fac),&
  !              sum(TEMP(:,:,k))/sum(TEMP_b(:,:,k)),&
  !              sum(SALT(:,:,k))/sum(SALT_b(:,:,k))

  do j=1,jmt
    do i=1,imt
      if ( DZT(i,j,k).ne.0.0 .and. DZT(i,j,k).ne.dz(k) ) then
        call state(dble(SALT_b(i,j,k)),dble(TEMP_b(i,j,k)),dble(pbar),PD_dble)
        PD_b(i,j,k) = real(PD_dble/1.0E03)
  !      if(j==1117 .and. k==5) then
  !        write(*,*) 'orig', TEMP(i,j,k), SALT(i,j,k), PD(i,j,k)
  !        write(*,*) 'a   ', TEMP(i,j,k), SALT(i,j,k), PD_a(i,j,k)
  !        write(*,*) 'b   ', TEMP_b(i,j,k), SALT_b(i,j,k), PD_b(i,j,k)
  !      endif
  ! >> newly calculated densities are higher which is sensible
      endif
    enddo !i
  enddo !j

  !=============================================================================
  ! 2.  cPKe
  !     recalculate PDW at center of fbc 
  !=============================================================================

  fac(:,:) = c0
  where(kmT==k)
    fac(:,:) = max(( dz_k(:,:)-DZT(:,:,k) ) / ( DZT(:,:,k)+DZT(:,:,k-1) ),c0)
    PDW_new(:,:,k) = PDW(:,:,k) + (2.0*PDW(:,:,k)-dz_k(:,:)*WPD(:,:,k))*fac(:,:)
  endwhere

  !do j=1,imt
  !  do i=1,imt
  !    if (fac(i,j)<c0) then
  !      write(*,*) i,j,k,fac(i,j),dz(k),DZT(i,j,k)
  !    endif
  !  enddo
  !enddo

enddo !k




do k = 1, km
  cPKt(:,:,k) = -g * PDW_new(:,:,k)*1.0E01
  cPKm(:,:,k) = -g * PD_a(:,:,k) * TTT_WVEL(:,:,k) * 1.0E01
  cPKe(:,:,k) = cPKt(:,:,k) - cPKm(:,:,k)
enddo !k

!write (*,*) 'Interpolated PD, pbc integral'
call vol_int(1,1,1,imt,jmt,km,cPKm,TAREA,DZT,cPKm_a_pint)
call vol_int(1,1,1,imt,jmt,km,cPKe,TAREA,DZT,cPKe_a_pint)
call vol_int(1,1,1,imt,jmt,km,cPKt,TAREA,DZT,cPKt_a_pint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKm,dz,TAREA,DZT,cPKm_a_fint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKe,dz,TAREA,DZT,cPKe_a_fint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKt,dz,TAREA,DZT,cPKt_a_fint)
write(*,*) 'density extrapolation, new PDW'
write(*,*) 'cPKm_a =    ',&
            cPKm_a_pint, 'W', cPKm_a_pint/-4.9E11,&
            cPKm_a_fint, 'W', cPKm_a_fint/-4.9E11
write(*,*) 'cPKe_a =    ', &
            cPKe_a_pint, 'W', cPKe_a_pint/ 7.3E11,&
            cPKe_a_fint, 'W', cPKe_a_fint/ 7.3E11
write(*,*) 'cPKt_a =    ',&
            cPKt_a_pint, 'W', cPKt_a_pint/ 2.4E11,&
            cPKt_a_fint, 'W', cPKt_a_fint/ 2.4E11

do k = 1, km
  cPKt(:,:,k) = -g * PDW_new(:,:,k)*1.0E01
  cPKm(:,:,k) = -g * PD_b(:,:,k) * TTT_WVEL(:,:,k) * 1.0E01
  cPKe(:,:,k) = cPKt(:,:,k) - cPKm(:,:,k)
enddo !k

!write (*,*) 'Interpolated TEMP/SALT, pbc integral'
call vol_int(1,1,1,imt,jmt,km,cPKm,TAREA,DZT,cPKm_b_pint)
call vol_int(1,1,1,imt,jmt,km,cPKe,TAREA,DZT,cPKe_b_pint)
call vol_int(1,1,1,imt,jmt,km,cPKt,TAREA,DZT,cPKt_b_pint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKm,dz,TAREA,DZT,cPKm_b_fint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKe,dz,TAREA,DZT,cPKe_b_fint)
call vol_int_fbc(1,1,1,imt,jmt,km,cPKt,dz,TAREA,DZT,cPKt_b_fint)
write(*,*) 'temp./salt extrapolation, new PDW'
write(*,*) 'cPKm_b =    ',&
            cPKm_b_pint, 'W', cPKm_b_pint/-4.9E11,&
            cPKm_b_fint, 'W', cPKm_b_fint/-4.9E11
write(*,*) 'cPKe_b =    ', &
            cPKe_b_pint, 'W', cPKe_b_pint/ 7.3E11,&
            cPKe_b_fint, 'W', cPKe_b_fint/ 7.3E11
write(*,*) 'cPKt_b =    ',&
            cPKt_b_pint, 'W', cPKt_b_pint/ 2.4E11,&
            cPKt_b_fint, 'W', cPKt_b_fint/ 2.4E11


write(*,*) ''
!write(*,*) sum(PDW), sum(PDW_new), sum(WPD)

!===============================================================================
!  OUTPUT
!===============================================================================


contains

!===============================================================================
subroutine vol_int_fbc(imin,jmin,kmin,imax,jmax,kmax,FIELD,dz,TAREA,DZT,INTEGRAL)
!
!     calculates volume integral (*[m^3]) for TTT-field via Kahan summation
!
implicit none

! input/output variables
integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: imin,jmin,kmin,imax,jmax,kmax
real,    dimension(km),         intent(in)  :: dz
real,    dimension(imt,jmt,km), intent(in)  :: DZT,FIELD
real,    dimension(imt,jmt),    intent(in)  :: TAREA
real,                           intent(out) :: INTEGRAL
real                                        :: c, y, t

INTEGRAL = 0.0
c = 0.0
do k = kmin,kmax
  do j = jmin,jmax
    do i = imin,imax
      if ( DZT(i,j,k).ne.0.0 ) then
        y = TAREA(i,j) * dz(k) * FIELD(i,j,k) - c
        t = INTEGRAL + y
        c = ( t - INTEGRAL ) - y
        INTEGRAL = t
      endif
    enddo !i
  enddo !j
enddo !k

INTEGRAL = INTEGRAL

end subroutine vol_int_fbc

!===============================================================================
!===============================================================================
end program pbc_test
