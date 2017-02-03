program WVEL_program
implicit none

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! this is to test whether the surface integral of  WVEL is 0
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!===============================================================================
!  variables
!===============================================================================

character*120 ::                                                               &
  geometry1_file,geometry2_file,input_folder,output_folder,file_base,          &
  file1, file2, file3, file1_out, file2_out, file3_out
integer                                         ::                             &
  imt,jmt,km,y,start_year,end_year,nt,nrec_WVEL,n,i,j,k,m,t,r
real                                            ::                             &
  WVEL_sint, WVEL_sint_abs, TTT_WVEL_sint, TTT_WVEL_sint_abs
real,             dimension(:),     allocatable ::                             &
  tdepth,dz,area
real,             dimension(:,:),   allocatable ::                             &
  Psi_v,Psi_v_avg,Psi_v_sum,Psi_v_std,Psi_v_mean, DXT,DYT,TAREA,DXU,DYU,UAREA
real,             dimension(:,:,:), allocatable ::                             &
  DZT,DZU,WVEL,TTT_WVEL
double precision, parameter                     ::                             &
c0 = 0., p125 = 0.125, p25 = 0.25, p5 = 0.5, c1 = 1.

imt        = 3600
jmt        = 2400
km         =   42

write (*,*) ''
write (*,*) '--- WVEL test ---'
write (*,*) ''

!===============================================================================
!  FILES
!===============================================================================
input_folder    = '/home/dijkbio/andre/LEC/input/'
geometry1_file  = trim(LEC_folder)//'input/geometry1'
geometry2_file  = trim(LEC_folder)//'input/geometry2'
!===============================================================================
!  GEOMETRY 
!===============================================================================

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
allocate( geometry2(6,km),tdz(km),dz(km),tdepth(km),area(km),p_z(km),vol(km) )
open(2,file=geometry2_file,access='direct',form='unformatted',recl=6,     &
       status='old')
do k=1,km
  read(2,rec=k) geometry2(:,k)
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
close(2)a


!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'


file1 = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_viebahn/tavg/t.t0.1_42l_nccs01.030102'
file2 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.5year.301'
file3 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.51.year.301'

file1_out = 'WVEL_output/WVEL_m'
file2_out = 'WVEL_output/WVEL_5'
file3_out = 'WVEL_output/WVEL_51'

open(1,file=file1,access='direct',form='unformatted',recl=imt*jmt,status='unknown')
open(2,file=file2,access='direct',form='unformatted',recl=imt*jmt,status='unknown')
open(3,file=file3,access='direct',form='unformatted',recl=imt*jmt,status='unknown')

open(11,file=file1_out,access='direct',form='unformatted',recl=5,status='unknown')
open(12,file=file2_out,access='direct',form='unformatted',recl=5,status='unknown')
open(13,file=file3_out,access='direct',form='unformatted',recl=5,status='unknown')

nrec_WVEL = 405

do r=1,3
  allocate( WVEL(imt,jmt,km), TTT_WVEL(imt,jmt,km) )

  !load WVEL
  do k=1,km
    read(r,rec=nrec_WVEL+k-1) WVEL(:,:,k)
  enddo

  ! create interpolated file
  call wtt2ttt(WVEL,DZT,TTT_WVEL)

  ! calculate and write out quantitites
  do k = 1,km
    call surf_int(1,1,imt,jmt,WVEL(:,:,k),TAREA,DZT(:,:,k),WVEL_sint)
    WVEL_sint = WVEL_sint/area(k)
    call surf_int(1,1,imt,jmt,abs(WVEL(:,:,k)),TAREA,DZT(:,:,k),WVEL_sint_abs)
    WVEL_sint_abs = WVEL_sint_abs/area(k)

    call surf_int(1,1,imt,jmt,TTT_WVEL(:,:,k),TAREA,DZT(:,:,k),TTT_WVEL_sint)
    TTT_WVEL_sint = TTT_WVEL_sint/area(k)
    call surf_int(1,1,imt,jmt,abs(TTT_WVEL(:,:,k)),TAREA,DZT(:,:,k),TTT_WVEL_sint_abs)
    TTT_WVEL_sint_abs = TTT_WVEL_sint_abs/area(k)

    write(r+10,rec=k) real(k), WVEL_sint, WVEL_sint_abs, TTT_WVEL_sint, TTT_WVEL_sint_abs
  enddo !k

  deallocate( WVEL, TTT_WVEL)
enddo !r


!===============================================================================
!===============================================================================
end program WVEL_program
