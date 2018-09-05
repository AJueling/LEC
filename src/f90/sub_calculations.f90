
!===============================================================================
subroutine mer_advection(DXT,DZT,VVEL,Eres,Eres_adv)
!  meridional advection of a energy reservoir term
implicit none

integer, parameter                          :: imt=3600,jmt=2400,km=42
integer                                     :: i,j,k
real,    dimension(imt,jmt,km), intent(in)  :: DZT, VVEL,Eres
real,    dimension(imt,jmt),    intent(in)  :: DXT
real,    dimension(imt,jmt),    intent(out) :: Eres_adv

!Eres_adv(:,:) = sum(Eres(:,:,:)*TTT_VVEL(:,:,:)*DZT(:,:,:),3)/1.e2 
do i=2,imt
  do j=2,jmt-1
    Eres_adv(i,j)=0.
    do k=1,km
    Eres_adv(i,j) = Eres_adv(i,j) +&
0.25*(Eres(i,j,k)+Eres(i,j+1,k))*(VVEL(i,j,k)*DZT(i,j,k)+VVEL(i-1,j,k)*DZT(i-1,j,k))/1.e2
    enddo
  enddo
enddo
i=1
  do j=2,jmt-1
    Eres_adv(1,j)=0.
    do k=1,km
    Eres_adv(i,j) = Eres_adv(i,j) +&
0.25*(Eres(i,j,k)+Eres(i,j+1,k))*(VVEL(i,j,k)*DZT(i,j,k)+VVEL(imt,j,k)*DZT(imt,j,k))/1.e2
    enddo
  enddo
end subroutine mer_advection
!===============================================================================
