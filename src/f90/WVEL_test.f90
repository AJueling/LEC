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
  geometry1_file,geometry2_file,LEC_folder,output_folder,file_base,kmt_file,   &
  file1, file2, file3, file1_out, file2_out, file3_out
integer                                         ::                             &
  imt,jmt,km,nrec_WVEL,nrec_RHO,nrec_PD,n,i,j,k,r,rec_length
integer,          dimension(:,:),   allocatable ::&
  kmT
real                                            ::                             &
  WVEL_usi,WVEL_abs_usi,WVEL_pos_usi,WVEL_neg_usi,&
  WVEL_RHO_usi,WVEL_pos_RHO_usi,WVEL_neg_RHO_usi,&
  WVEL_PD_usi,WVEL_pos_PD_usi,WVEL_neg_PD_usi,&
  WVEL_RHO_new_usi,WVEL_pos_RHO_new_usi,WVEL_neg_RHO_new_usi,&
  WVEL_PD_new_usi,WVEL_pos_PD_new_usi,WVEL_neg_PD_new_usi,&
  WVEL_wsi,WVEL_abs_wsi,WVEL_pos_wsi,WVEL_neg_wsi,&
  WVEL_RHO_wsi,WVEL_pos_RHO_wsi,WVEL_neg_RHO_wsi,&
  WVEL_PD_wsi,WVEL_pos_PD_wsi,WVEL_neg_PD_wsi,&
  WVEL_RHO_new_wsi,WVEL_pos_RHO_new_wsi,WVEL_neg_RHO_new_wsi,&
  WVEL_PD_new_wsi,WVEL_pos_PD_new_wsi,WVEL_neg_PD_new_wsi
real,             dimension(:),     allocatable ::                             &
  tdepth,dz,area
real,             dimension(:,:),   allocatable ::                             &
  DXT,DYT,TAREA,geometry2
real,             dimension(:,:,:), allocatable ::                             &
  DZT,WVEL,WTT_WVEL,RHO,PD,WVEL_pos,WVEL_neg,DZT_pos,DZT_neg,RHO_new,PD_new
real, parameter                                 ::                             &
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
LEC_folder    = '/home/dijkbio/andre/LEC/'
geometry1_file  = trim(LEC_folder)//'input/geometry1'
geometry2_file  = trim(LEC_folder)//'input/geometry2'
!===============================================================================
!  GEOMETRY 
!===============================================================================

! kmT
kmt_file ='/home/dijkbio/andre/LEC/input/kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black'
allocate( kmT(imt,jmt) )
inquire (iolength=rec_length) kmT
open    (1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,   &
           status='unknown')
  read  (1,rec=1) kmT
close   (1)

allocate( DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), DZT(imt,jmt,km) )
open(2,file=geometry1_file,access='direct',form='unformatted',recl=imt*jmt,    &
       status='old')
read(2,rec=1) DXT ! [m]
read(2,rec=3) DYT ! [m]
read(2,rec=5) TAREA ! [m^2]
do k=1,km
  read(2,rec=6+k) DZT(:,:,k) ! [m]
enddo
close(2)

! read 1D geometry fields
! all real: k,dz[m],tdepth[m],area[m^2],p[bar]
allocate( geometry2(6,km),dz(km),tdepth(km),area(km) )
open(2,file=geometry2_file,access='direct',form='unformatted',recl=6,          &
       status='old')
do k=1,km
  read(2,rec=k) geometry2(:,k)
enddo
dz(:)     = geometry2(2,:)
tdepth(:) = geometry2(3,:)
area(:)   = geometry2(4,:)
close(2)


!===============================================================================
!  CALCULATIONS
!===============================================================================

write (*,*) ''
write (*,*) 'CALCULATIONS'

file1 = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_viebahn/tavg/t.t0.1_42l_nccs01.030102'
file2 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.year.301'
file3 = '/projects/0/samoc/jan/Andree/t.t0.1_42l_nccs01.tavg.5year.301'

output_folder = '/home/dijkbio/andre/LEC/results/WVEL_test/'
file1_out     = trim(output_folder)//'WVEL_m'
file2_out     = trim(output_folder)//'WVEL_1'
file3_out     = trim(output_folder)//'WVEL_5'

open(1,file=file1,access='direct',form='unformatted',recl=imt*jmt,status='unknown')
open(2,file=file2,access='direct',form='unformatted',recl=imt*jmt,status='unknown')
open(3,file=file3,access='direct',form='unformatted',recl=imt*jmt,status='unknown')

open(11,file=file1_out,access='direct',form='unformatted',recl=33,status='unknown')
open(12,file=file2_out,access='direct',form='unformatted',recl=33,status='unknown')
open(13,file=file3_out,access='direct',form='unformatted',recl=33,status='unknown')

nrec_WVEL = 405
nrec_RHO  = 306
nrec_PD   = 995

do r=1,3
  write(*,*) r
  allocate( WTT_WVEL(imt,jmt,km), WVEL(imt,jmt,km), &
            RHO(imt,jmt,km), PD(imt,jmt,km) )

  !load WVEL
  do k=1,km
    read(r,rec=nrec_WVEL+k-1) WTT_WVEL(:,:,k)
    read(r,rec=nrec_RHO+k-1)  RHO(:,:,k)
    read(r,rec=nrec_PD+k-1)   PD(:,:,k)
  enddo

  ! create interpolated file
  call wtt2ttt(WTT_WVEL,DZT,WVEL)

  ! creating positive and negative WVEL fields
  allocate( WVEL_pos(imt,jmt,km), WVEL_neg(imt,jmt,km), &
            DZT_pos(imt,jmt,km), DZT_neg(imt,jmt,km),  &
            RHO_new(imt,jmt,km), PD_NEW(imt,jmt,km) )

  WVEL_pos = WVEL
  DZT_pos  = DZT
  WVEL_neg = WVEL
  DZT_neg  = DZT

  RHO_new  = RHO
  PD_new   = PD

  ! calculate and write out quantitites
  do k = 1,km

    where (WVEL(:,:,k)<c0)
      WVEL_pos(:,:,k) = c0
      DZT_pos(:,:,k)  = c0
    elsewhere (WVEL(:,:,k)>c0)
      WVEL_neg(:,:,k) = c0
      DZT_neg(:,:,k)  = c0
    endwhere

    if (k>1) then
      do i=1,imt
        do j=1,jmt
          !if (DZT(i,j,k)<dz(k) .and. DZT(i,j,k).ne.c0) then
          if (kmT(i,j)==k) then
            RHO_new(i,j,k) = (RHO(i,j,k)-RHO(i,j,k-1)) &
                            /(DZT(i,j,k)+dz(k-1)) &
                            *(dz(k)+dz(k-1))
            PD_new(i,j,k)  = (PD(i,j,k)-PD(i,j,k-1)) &
                            /(DZT(i,j,k)+dz(k-1)) &
                            *(dz(k)+dz(k-1))
          endif
        enddo
      enddo
    endif
    ! 1. unweighted surface integrals
    ! 1.1. WVEL
    call surf_int(1,1,imt,jmt,               WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_usi)
    ! 1.2. abs(WVEL)
    call surf_int(1,1,imt,jmt,          abs(WVEL(:,:,k)),TAREA,DZT(:,:,k),&
                  WVEL_abs_usi)
    ! 1.3. WVEL_pos
    call surf_int(1,1,imt,jmt,           WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_usi)
    ! 1.4. WVEL_neg
    call surf_int(1,1,imt,jmt,           WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_usi)
    ! 1.5. WVEL     * RHO
    call surf_int(1,1,imt,jmt,RHO(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_RHO_usi)
    ! 1.6. WVEL_pos * RHO
    call surf_int(1,1,imt,jmt,RHO(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_RHO_usi)
    ! 1.7. WVEL_neg * RHO
    call surf_int(1,1,imt,jmt,RHO(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_RHO_usi)
    ! 1.8. WVEL     * PD
    call surf_int(1,1,imt,jmt,PD(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_PD_usi)
    ! 1.9. WVEL_pos * PD
    call surf_int(1,1,imt,jmt,PD(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_PD_usi)
    ! 1.10. WVEL_neg * PD
    call surf_int(1,1,imt,jmt,PD(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_PD_usi)
    ! 1.11. WVEL     * RHO_new
    call surf_int(1,1,imt,jmt,RHO_new(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_RHO_new_usi)
    ! 1.12. WVEL_pos * RHO_new
    call surf_int(1,1,imt,jmt,RHO_new(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_RHO_new_usi)
    ! 1.13. WVEL_neg * RHO_new
    call surf_int(1,1,imt,jmt,RHO_new(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_RHO_new_usi)
    ! 1.14. WVEL     * PD_new
    call surf_int(1,1,imt,jmt,PD_new(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_PD_new_usi)
    ! 1.15. WVEL_pos * PD_new
    call surf_int(1,1,imt,jmt,PD_new(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_PD_new_usi)
    ! 1.16. WVEL_neg * PD_new
    call surf_int(1,1,imt,jmt,PD_new(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_PD_new_usi)

    ! 2. weighted surface integrals
    ! 2.1. WVEL
    call surf_int_3D(1,1,imt,jmt,               WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_wsi)
    ! 2.2. abs(WVEL)
    call surf_int_3D(1,1,imt,jmt,          abs(WVEL(:,:,k)),TAREA,DZT(:,:,k),&
                  WVEL_abs_wsi)
    ! 2.3. WVEL_pos
    call surf_int_3D(1,1,imt,jmt,           WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_wsi)
    ! 2.4. WVEL_neg
    call surf_int_3D(1,1,imt,jmt,           WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_wsi)
    ! 2.5. WVEL     * RHO
    call surf_int_3D(1,1,imt,jmt,RHO(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_RHO_wsi)
    ! 2.6. WVEL_pos * RHO
    call surf_int_3D(1,1,imt,jmt,RHO(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_RHO_wsi)
    ! 2.7. WVEL_neg * RHO
    call surf_int_3D(1,1,imt,jmt,RHO(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_RHO_wsi)
    ! 2.8. WVEL     * PD
    call surf_int_3D(1,1,imt,jmt,PD(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_PD_wsi)
    ! 2.9. WVEL_pos * PD
    call surf_int_3D(1,1,imt,jmt,PD(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_PD_wsi)
    ! 2.10. WVEL_neg * PD
    call surf_int_3D(1,1,imt,jmt,PD(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_PD_wsi)
    ! 2.11. WVEL     * RHO_new
    call surf_int_3D(1,1,imt,jmt,RHO_new(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_RHO_new_wsi)
    ! 2.12. WVEL_pos * RHO_new
    call surf_int_3D(1,1,imt,jmt,RHO_new(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_RHO_new_wsi)
    ! 2.13. WVEL_neg * RHO_new
    call surf_int_3D(1,1,imt,jmt,RHO_new(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_RHO_new_wsi)
    ! 2.14. WVEL     * PD_new
    call surf_int_3D(1,1,imt,jmt,PD_new(:,:,k)*    WVEL(:,:,k),TAREA,DZT(:,:,k),&
                  WVEL_PD_new_wsi)
    ! 2.15. WVEL_pos * PD_new
    call surf_int_3D(1,1,imt,jmt,PD_new(:,:,k)*WVEL_pos(:,:,k),TAREA,DZT_pos(:,:,k),&
                  WVEL_pos_PD_new_wsi)
    ! 2.16. WVEL_neg * PD_new
    call surf_int_3D(1,1,imt,jmt,PD_new(:,:,k)*WVEL_neg(:,:,k),TAREA,DZT_neg(:,:,k),&
                  WVEL_neg_PD_new_wsi)

    write(r+10,rec=k) real(k),&
                      WVEL_usi,WVEL_abs_usi,WVEL_pos_usi,WVEL_neg_usi,&
                      WVEL_RHO_usi,WVEL_pos_RHO_usi,WVEL_neg_RHO_usi,&
                      WVEL_PD_usi,WVEL_pos_PD_usi,WVEL_neg_PD_usi,&
                      WVEL_RHO_new_usi,WVEL_pos_RHO_new_usi,WVEL_neg_RHO_new_usi,&
                      WVEL_PD_new_usi,WVEL_pos_PD_new_usi,WVEL_neg_PD_new_usi,&
                      WVEL_wsi,WVEL_abs_wsi,WVEL_pos_wsi,WVEL_neg_wsi,&
                      WVEL_RHO_wsi,WVEL_pos_RHO_wsi,WVEL_neg_RHO_wsi,&
                      WVEL_PD_wsi,WVEL_pos_PD_wsi,WVEL_neg_PD_wsi,&
                      WVEL_RHO_new_wsi,WVEL_pos_RHO_new_wsi,WVEL_neg_RHO_new_wsi,&
                      WVEL_PD_new_wsi,WVEL_pos_PD_new_wsi,WVEL_neg_PD_new_wsi
  enddo !k

  deallocate( WTT_WVEL, WVEL, RHO, PD, RHO_new, PD_new,&
              WVEL_pos, WVEL_neg, DZT_pos, DZT_neg )
enddo !r


!===============================================================================
!===============================================================================
end program WVEL_program
