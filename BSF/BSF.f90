!  Compile on Huygens with
!    xlf90 -O3 -o SALT_BUDGET SALT_BUDGET.f90 -lnetcdf -lnetcdff
!
   program ENERGY_BUDGET
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Program based on SALT_BUDGET.f90 version from January 2014
!  by ( Michael  Kliphuis & Matthijs den Toom (IMAU, August 2011), with comments by Dewi (IMAU, January 2014):
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   use netcdf 
   implicit none

!=====================================================================================================================================================
!  variables
!=====================================================================================================================================================

   ! user input 
   integer :: ncase, imt,jmt,km

   integer :: nrec_UVEL

   character*120 :: grid_file,kmt_file,in_depths,pbc_file,tavg_file,outputfile
   character*3   :: year

   ! program output
   real :: BSF_drake, BSF_gyre

   ! internal model variables 
   integer :: rec_length, i, ip, j, k, ncid, maskid

   integer, dimension(:,:), allocatable :: kmT

   integer, dimension(:,:,:), allocatable :: TUU_mask

   real, dimension(:,:), allocatable ::    UVEL_k, BSF

   real, dimension(:,:,:), allocatable ::  UVEL

   double precision, dimension(:), allocatable ::  tdepth, dz

   double precision, dimension(:,:), allocatable :: HTN, HTE, WORK, DXT, DYT,&
   TAREA, DXU, DYU, UAREA, DZBC, HUW, HUS, WORK2, WORK3 

   real, dimension(:,:,:), allocatable :: DZT

   ! Parameters
   double precision, parameter ::  c0 = 0., p5 = 0.5, c1 = 1.
   ! salinity factor from http://www.cesm.ucar.edu/models/cesm1.0/cesm/
   ! cesmBbrowser/html_code/pop/constants.F90.html "flux (kg/m^2/s) to salt flux
   ! (msu*cm/s)"


   !netCDF
   include 'netcdf.inc'

   write (*,*) ''

!=====================================================================================================================================================
!  user input
!=====================================================================================================================================================

   read  (*,*) imt, jmt, km
   read  (*,'(a120)') grid_file
   read  (*,'(a120)') kmt_file
   read  (*,'(a120)') in_depths
   read  (*,'(a120)') pbc_file
   read  (*,'(a120)') tavg_file
   read  (*,'(a3)') year
   read  (*,*) nrec_UVEL
   read  (*,'(a120)') outputfile



   write(*,*) 'Year:', year

!=====================================================================================================================================================
!  read and create horizontal grid spacing, define TAREA
!=====================================================================================================================================================

   allocate(                                                                   &
   HTN(imt,jmt), HTE(imt,jmt), HUW(imt,jmt), HUS(imt,jmt), WORK(imt,jmt),      &
   WORK2(imt,jmt), WORK3(imt,jmt), DXT(imt,jmt), DYT(imt,jmt), TAREA(imt,jmt), &
   DXU(imt,jmt), DYU(imt,jmt), UAREA(imt,jmt) )

   inquire (iolength = rec_length) HTN
   open (1,file=grid_file,access='direct',form='unformatted',recl=rec_length,status='unknown')
   read (1,rec=3) HTN ! [cm]
   read (1,rec=4) HTE ! [cm]
   close(1)

   ! Why do we do this?
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

!   deallocate( WORK, WORK2, WORK3 )

!=====================================================================================================================================================
!  read bathymetry and depth level variables, define DZT
!=====================================================================================================================================================

   allocate( dz(km), tdepth(km), kmT(imt,jmt), DZBC(imt,jmt) )
   allocate( DZT(imt,jmt,km), TUU_mask(imt,jmt,km) )

   ! kmT
   inquire (iolength=rec_length) kmT
   open    (1,file=kmt_file,access='direct',form='unformatted',recl=rec_length,status='unknown')
   read    (1,rec=1) kmT
   close   (1)

   ! dz
   open    (1,file=in_depths,status='old')
   do k = 1, km
   read  (1,*) dz(k) ! [cm]
   enddo
   close   (1)

   ! partial bottom cell depths
   inquire (iolength=rec_length) DZBC
   open    (1,file=pbc_file,access='direct',form='unformatted',recl=rec_length,status='unknown')
   read    (1,rec=1) DZBC ! [cm]
   close   (1)

   ! create DZT-field, also used as 3D TTT-mask
   do k=1,km
     where(kmT>k)
       DZT(:,:,k) = dz(k)
     elsewhere(kmT==k)
       DZT(:,:,k) = DZBC
     elsewhere
       DZT(:,:,k) = c0
     endwhere
   enddo
   deallocate( kmT, DZBC )
   
!=====================================================================================================================================================
!  read fields
!=====================================================================================================================================================


   ! 3D
   allocate(UVEL_k(imt,jmt),UVEL(imt,jmt,km))

   !open file
   inquire (iolength=rec_length) UVEL_k
   open (1,file=tavg_file,access='direct',form='unformatted',recl=rec_length,status='unknown')

   !then read 3-D fields, storing full salinity field, and velocity field for sections
   do k = 1, km
     read (1,rec=nrec_UVEL+k-1)  UVEL_k     ! [cm/s]
     where(UVEL_k>100.0)
       UVEL_k = c0
     elsewhere(UVEL_k<-100.0)
       UVEL_k = c0 
     endwhere
     UVEL(:,:,k)               = UVEL_k
   enddo


   deallocate( UVEL_k )

   !write (*,*) 'UVEL(500,1000,1): ', UVEL(500,1000,1)
   !write (*,*) 'UVEL(1000,500,1): ', UVEL(1000,500,1)
   !write (*,*) 'UVEL(1,1,1): ', UVEL(1,1,1)

   write (*,*) 'file: ', tavg_file

!=====================================================================================================================================================
!  calculate BSF
!=====================================================================================================================================================

   allocate(BSF(imt,jmt))

   BSF(:,:) = c0
   do k = 1, km
     BSF(:,:) = BSF(:,:) + UVEL(:,:,k)*dz(k)/1.0E04 ! depth integral in [m^2/s]
   enddo


   BSF(:,1) = BSF(:,1)*DYU(:,1)/1.0E02 ! now volume transport in [m^3/s]
   do j = 2,jmt/2
   !only necessary up to equator; north of it the grid is distorted anyways
     BSF(:,j) = BSF(:,j)*DYU(:,j)/1.0E02 + BSF(:,j-1)
   enddo


   !BSF(:,:) = BSF(:,:) - BSF(446,521)


   BSF_drake = MAXVAL(BSF(446,311:521))
   BSF_gyre  = MINVAL(BSF(600:1300,100:500))


!=====================================================================================================================================================
!   FILE OUTPUT
!=====================================================================================================================================================

100 FORMAT (132(A,","),A)
101 FORMAT (132(E15.7,","),E15.7)
   open(1,file=outputfile,form='formatted',status='old',action='write',position='append')
   write (1,*) year, BSF_drake, BSF_gyre
   close (1)

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 

      subroutine check(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
        write(*,*) trim(nf90_strerror(status))
        stop
      endif
 
      end subroutine check
 
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   end program

!
!  Below are some useful (model-native?) subroutines
!

!***********************************************************************

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
        end do
      end do

      end subroutine e_rshift

!***********************************************************************

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
        end do
      end do

      end subroutine s_rshift

!***********************************************************************

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

!***********************************************************************

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

!***********************************************************************
