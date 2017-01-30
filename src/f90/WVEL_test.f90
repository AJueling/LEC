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
  grid_file,kmt_file,in_depths,pbc_file,input_folder,output_folder,file_base,  &
  file1, file2, file3, file1_out, file2_out, file3_out
integer                                         ::                             &
  imt,jmt,km,y,start_year,end_year,nt,nrec_WVEL,n,rec_length,i,j,k,m,t,r
integer,          dimension(:,:),   allocatable ::                             &
  kmT
double precision                                ::                             &
  volume 
double precision, dimension(:),     allocatable ::                             &
  dz,area
double precision, dimension(:,:),   allocatable ::                             &
  HTN,HTE,DXT,DYT,TAREA,DXU,DYU,UAREA,DZBC,HUW,HUS,WORK,WORK2,WORK3
real                                            ::                             &
  WVEL_sint, WVEL_sint_abs, TTT_WVEL_sint, TTT_WVEL_sint_abs
real,             dimension(:),     allocatable ::                             &
  tdepth
real,             dimension(:,:),   allocatable ::                             &
  Psi_v,Psi_v_avg,Psi_v_sum,Psi_v_std,Psi_v_mean
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

grid_file       = trim(input_folder)//'grid.3600x2400.fob.da'
kmt_file        = trim(input_folder)//'kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'
in_depths       = trim(input_folder)//'in_depths.42.dat'
pbc_file        = trim(input_folder)//'dzbc_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'
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
  call wtt2ttt(imt,jmt,km,WVEL,DZT,TTT_WVEL)

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
