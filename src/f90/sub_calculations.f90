ubroutine mer_advection(imt,jmt,km,DXT,DZT,TTT_VVEL,Eres,Eres_adv)
!
!  meridional advection of a energy reservoir term
!
implicit none

integer,          intent(in)                         :: imt, jmt, km
real,             dimension(imt,jmt,km), intent(in)  :: DZT, TTT_VVEL,
Eres
double precision, dimension(imt,jmt),    intent(in)  :: DXT
real,             dimension(imt,jmt),    intent(out) :: Eres_adv

Eres_adv(:,:) = sum(Eres(:,:,:)*TTT_VVEL(:,:,:)*DZT(:,:,:),3) * DXT(:,:)

end subroutine mer_advection

subroutine osf(imt,jmt,km,DXU,DZU,VVEL,WVEL,Psi_v,Psi_w)
!
!  overturning stream function calculated both via VVEL and WVEL
!
implicit none

integer,          intent(in)                        :: imt, jmt, km
integer                                             :: l
real,             dimension(imt,jmt,km)             :: DZU, VVEL, WVEL
double precision, dimension(imt,jmt),   intent(in)  :: DXU
real,             dimension(jmt,km),    intent(out) :: Psi_v, Psi_w


! preparing integral dz/dy + executing integral dx
do k = 1,km
  Psi_v(:,k) = -sum(VVEL(:,:,k)*DZU(:,:,k)*DXU(:,:)/1.0E06,1)
  Psi_w(:,k) =  sum(WVEL(:,:,k)*DYU(:,:)  *DXU(:,:)/1.0E06,1)
enddo

! Psi_v (depth integral dz bottom to depth z)
do k = 1,km-1
  Psi_v(:,km-k) = Psi_v(:,km-k) + Psi_v(:,km-k+1)
enddo

! Psi_w (meridional integral dy South Pole to latitude j)
do j = 2,jmt
  Psi_w(j,:) = Psi_w(j,:) + Psi_w(j-1,:)
enddo

end subroutine osf

