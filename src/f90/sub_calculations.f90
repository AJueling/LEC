
!===============================================================================
subroutine mer_advection(DXT,DZT,TTT_VVEL,Eres,Eres_adv)
!  meridional advection of a energy reservoir term
implicit none

integer, parameter                          :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt,km), intent(in)  :: DZT, TTT_VVEL,Eres
real,    dimension(imt,jmt),    intent(in)  :: DXT
real,    dimension(imt,jmt),    intent(out) :: Eres_adv

Eres_adv(:,:) = sum(Eres(:,:,:)*TTT_VVEL(:,:,:)*DZT(:,:,:),3) * DXT(:,:)

end subroutine mer_advection
!===============================================================================
