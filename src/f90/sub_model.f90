!===============================================================================
subroutine load_3D_field(nfile,nrec,FIELD)
!
!  loads 3D field (starting with record number nrec) from file with
!  number nfile
!
implicit none

integer                                     :: k
integer, parameter                          :: imt=3600,jmt=2400,km=42
integer,                        intent(in)  :: nfile,nrec
real,    dimension(imt,jmt)                 :: WORK
real,    dimension(imt,jmt,km), intent(out) :: FIELD

do k = 1,km
  read (nfile,rec=nrec+k-1) WORK
  FIELD(:,:,k) = WORK
enddo

end subroutine load_3d_field
!===============================================================================

!===============================================================================
subroutine load_2D_field(nfile,nrec,FIELD)
!
!  loads 2D field (starting with record number nrec) from file with
!  number nfile
!
implicit none

integer, parameter                       :: imt=3600,jmt=2400,km=42
integer,                     intent(in)  :: nfile,  nrec
real,    dimension(imt,jmt), intent(out) :: FIELD

read (nfile,rec=nrec) FIELD

end subroutine load_2D_field
!===============================================================================

!===============================================================================
subroutine e_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from east (double precision arrays)
!-----------------------------------------------------------------------
integer                                           :: i,j,n,ip
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result

do j=1,jmt
  do i=1,imt
    ip = i + 1
    if(i == imt) ip = 1
    XOUT(i,j) = X(ip,j)
  enddo
enddo

end subroutine e_rshift
!===============================================================================

!===============================================================================
subroutine s_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from south (real arrays)
!-----------------------------------------------------------------------
integer                                           :: i,j,n,jm
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result

do j=1,jmt
  jm = j - 1
  if(j == 1) jm = 1
  do i=1,imt
    XOUT(i,j) = X(i,jm)
  enddo
enddo

end subroutine s_rshift
!===============================================================================

!===============================================================================
subroutine w_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from west (real arrays)
!-----------------------------------------------------------------------
integer                                           :: i,j,n,im
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result

do j=1,jmt
  do i=1,imt
    im = i - 1
    if(i == 1) im = imt
    XOUT(i,j) = X(im,j)
  end do
end do

end subroutine w_rshift
!===============================================================================

!===============================================================================
subroutine n_rshift(XOUT,X,imt,jmt)
implicit none
!-----------------------------------------------------------------------
!     shift from north (double precision arrays)
!-----------------------------------------------------------------------
integer                                           :: i,j,n,jp
integer,                              intent(in)  :: imt,jmt
double precision, dimension(imt,jmt), intent(in)  :: X    ! array to be shifted
double precision, dimension(imt,jmt), intent(out) :: XOUT ! shifted result

do j=1,jmt
  jp = j+1
  if(j == jmt) jp = jmt
  do i=1,imt
    XOUT(i,j) = X(i,jp)
  end do
end do

end subroutine n_rshift
!===============================================================================

