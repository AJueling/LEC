!===============================================================================
!  VERTICAL
!===============================================================================

!===============================================================================
subroutine wtt2ttt(W_FIELD,DZT,T_FIELD)
!
!     interpolated real field from WTT-grid to TTT-grid
!     interpolating between two W-levels
!     T-point is located exactly in the middle of two W-points
!     T_FIELD(i,j,k) = W_FIELD(i,j,k) + ( W_FIELD(i,j,k+1) -
!     W_FIELD(i,j,k) ) &
!                                       / dz(k) * dz(k) / 2
!
implicit none

! input/output variables
integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt,km), intent(in)  :: DZT
real,    dimension(imt,jmt,km), intent(in)  :: W_FIELD ! WTT
real,    dimension(imt,jmt,km), intent(out) :: T_FIELD ! TTT
real,    parameter                          :: p5=0.5, c0=0.0, c3=3.0

do k = 1,km-1
  do j = 1,jmt
    do i = 1,imt
      if ( DZT(i,j,k).ne.c0 ) then
        T_FIELD(i,j,k) = p5 * ( W_FIELD(i,j,k) + W_FIELD(i,j,k+1) )
      endif
    enddo !i
  enddo !j
enddo !k

! bottom level
do j = 1,jmt
  do i = 1,imt
    if ( DZT(i,j,km).ne.c0 ) then
      ! extrapolating to last level with gradient between last two
      ! W-levels
      ! T_FIELD(i,j,km) = W_FIELD(i,j,km) 
      !                 + (W_FIELD(i,j,km)-W_FIELD(i,j,km-1)) / dz(k) *
      !                 dz(k)/2
      !T_FIELD(i,j,km) = p5 * (c3*W_FIELD(i,j,km) - W_FIELD(i,j,km-1) )
      T_FIELD(i,j,km) = p5 * W_FIELD(i,j,km)
    endif
  enddo !i
enddo !j

!do k = 1,km
!  call surf_int(1,1,imt,jmt,W_FIELD(:,:,:),TAREA,DZT(:,:,k),A)
!  call surf_int(1,1,imt,jmt,T_FIELD(:,:,k),TAREA,DZT(:,:,k),B)
!  write (*,*) k, ' wtt2ttt orig', A
!  write (*,*) k, ' wtt2ttt new ', B
!  write (*,*) k, ' wtt2ttt factor      ', B/A
!enddo !k
end subroutine wtt2ttt
!===============================================================================


!===============================================================================
!  HORIZONTAL
!===============================================================================

!===============================================================================
subroutine uu2tt(DZT_k,DZU_k,TAREA,UAREA,FIELD_k,NEW_FIELD_k)
!
!     interpolates real 2D field from UU-grid to TT-grid on T-depth
!     levels
!
implicit none

! input/output variables
integer                                  :: i,j,k
integer, parameter                       :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
real,    dimension(imt,jmt), intent(in)  :: TAREA,UAREA
real,    dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
real,    dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT
real,    parameter                       :: c0=0.0
! southernmost gridpoints (j=1) on TTT-grid must be 0,
! as there would not be UU gridpoints south of it otherwise
do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    NEW_FIELD_k(1,j) = 0.25                                                    &
                  * ( FIELD_k(1  ,j  ) * UAREA(1  ,j  ) * DZU_k(1  ,j  )       &
                    + FIELD_k(imt,j  ) * UAREA(imt,j  ) * DZU_k(imt,j  )       &
                    + FIELD_k(1  ,j-1) * UAREA(1  ,j-1) * DZU_k(1  ,j-1)       &
                    + FIELD_k(imt,j-1) * UAREA(imt,j-1) * DZU_k(imt,j-1))      &
                  / TAREA(1,j) / DZT_k(1,j)
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      NEW_FIELD_k(i,j) = 0.25                                                  &
                   * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  ) * DZU_k(i  ,j  )      &
                     + FIELD_k(i-1,j  ) * UAREA(i-1,j  ) * DZU_k(i-1,j  )      &
                     + FIELD_k(i  ,j-1) * UAREA(i  ,j-1) * DZU_k(i  ,j-1)      &
                     + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) * DZU_k(i-1,j-1) )    &
                   / TAREA(i,j) / DZT_k(i,j)
    endif
  enddo !i
enddo !j

end subroutine uu2tt
!===============================================================================

!===============================================================================
subroutine uu2tt_3D(DZT,DZU,TAREA,UAREA,TUU_FIELD,TTT_FIELD)
implicit none
!  interpolates a 3D field from UU to TT

integer                                  :: i,j,k
integer, parameter                       :: imt=3600,jmt=2400,km=42
real, dimension(imt,jmt)                 :: WORK
real, dimension(imt,jmt),    intent(in)  :: TAREA, UAREA
real, dimension(imt,jmt,km), intent(in)  :: TUU_FIELD, DZU, DZT
real, dimension(imt,jmt,km), intent(out) :: TTT_FIELD

do k = 1,km
  call uu2tt(DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,TUU_FIELD(:,:,k),WORK)
  TTT_FIELD(:,:,k) = WORK
enddo

end subroutine uu2tt_3D
!===============================================================================

!===============================================================================
subroutine uu2tt_scalar(DZT_k,DZU_k,TAREA,UAREA,FIELD_k,NEW_FIELD_k)
!     interpolates real scalar 2D field (as opposed to velocity fields
!     which  
!     can be related to volume transports) from UU-grid to TT-grid on
!     T-depth 
!     levels
implicit none

! input/output variables
integer                                  :: i,j,k
integer, parameter                       :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
real,    dimension(imt,jmt), intent(in)  :: TAREA,UAREA
real,    dimension(imt,jmt), intent(in)  :: FIELD_k     ! (T)UU
real,    dimension(imt,jmt), intent(out) :: NEW_FIELD_k ! (T)TT
real,    parameter                       :: c0=0.0

! southernmost gridpoints (j=1) on TTT-grid must be 0,
! as there would not be UU gridpoints south of it otherwise
do j = 2, jmt
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    NEW_FIELD_k(1,j) = 0.25 * ( FIELD_k(1  ,j  ) * UAREA(1  ,j  )              &
                              + FIELD_k(imt,j  ) * UAREA(imt,j  )              &
                              + FIELD_k(1  ,j-1) * UAREA(1  ,j-1)              &
                              + FIELD_k(imt,j-1) * UAREA(imt,j-1) )            &
                            / TAREA(1,j)
  endif
  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      NEW_FIELD_k(i,j) = 0.25 * ( FIELD_k(i  ,j  ) * UAREA(i  ,j  )            &
                                + FIELD_k(i-1,j  ) * UAREA(i-1,j  )            &
                                + FIELD_k(i  ,j-1) * UAREA(i  ,j-1)            &
                                + FIELD_k(i-1,j-1) * UAREA(i-1,j-1) )          &
                              / TAREA(i,j)
    endif
  enddo !i
enddo !j

! results in very small correction O(-4)

end subroutine uu2tt_scalar
!===============================================================================

!===============================================================================
subroutine uu2tt_scalar_3D(DZT,DZU,TAREA,UAREA,TUU_FIELD,TTT_FIELD)
implicit none
!
!  interpolates a 3D field from UU to TT
!

integer                                     :: i,j,k
integer, parameter                          :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt)                 :: WORK
real,    dimension(imt,jmt),    intent(in)  :: TAREA, UAREA
real,    dimension(imt,jmt,km), intent(in)  :: TUU_FIELD, DZU, DZT
real,    dimension(imt,jmt,km), intent(out) :: TTT_FIELD

do k = 1,km
  call uu2tt_scalar(DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,TUU_FIELD(:,:,k),WORK)
  TTT_FIELD(:,:,k) = WORK
enddo

end subroutine uu2tt_scalar_3D
!===============================================================================


!===============================================================================
!  MOMENTUM FLUX
!===============================================================================

!===============================================================================
subroutine interp_mom_fluxes(DZT_k,DZU_k,DXT,DXU,DYT,DYU,&
                                   UE_FIELD,VN_FIELD,UE_NEW_FIELD,VN_NEW_FIELD)
!
!     interpolates momentum fluxes from east side of U-cell to momentum
!     fluxes
!     at center of T-cells, for this a weighted average of the flux
!     North/South
!     and East/West of the T-cell center point is used
!
implicit none

! input/output variables
integer                                  :: i,j,k
integer, parameter                       :: imt=3600,jmt=2400,km=42
real,    dimension(imt,jmt), intent(in)  :: DZU_k,DZT_k
real,    dimension(imt,jmt), intent(in)  :: DXT,DXU,DYT,DYU
real,    dimension(imt,jmt), intent(in)  :: UE_FIELD ! (T)UU
real,    dimension(imt,jmt), intent(in)  :: VN_FIELD ! (T)UU
real,    dimension(imt,jmt), intent(out) :: UE_NEW_FIELD ! (T)TT
real,    dimension(imt,jmt), intent(out) :: VN_NEW_FIELD ! (T)TT
real,    parameter                       :: c0=0.0
real :: A, B ! surface integrals, used for comparing old and new fields 

! southernmost gridpoints (j=1) on TTT-grid must be 0, 
! as there would not be UU gridpoints south of it otherwise

! east flux
do j = 2, jmt-1
  ! westernmost gridpoints
  if ( DZT_k(1,j).ne.c0 ) then
    if ( DZT_k(1,j-1).ne.c0 .or. DZT_k(1,j+1).ne.c0  ) then
      UE_NEW_FIELD(1,j) = 0.5 *&
       ( UE_FIELD(imt,j  ) * DYU(imt,j  )&
       + UE_FIELD(imt,j-1) * DYU(imt,j-1) )&
       / DYT(1,j)
    else
      UE_NEW_FIELD(1,j) = 0.0
    endif
  endif

  ! most of the ocean
  do i = 2,imt
    if ( DZT_k(i,j).ne.c0 ) then
      if ( DZT_k(i,j-1).ne.c0 .or. DZT_k(i,j+1).ne.c0  ) then
        UE_NEW_FIELD(i,j) = 0.5 *&
         ( UE_FIELD(i-1,j  ) * DYU(i-1,j  )&
         + UE_FIELD(i-1,j-1) * DYU(i-1,j-1) )&
         / DYT(i,j)
      else
        UE_NEW_FIELD(i,j) = 0.0
      endif
    endif
  enddo !i
enddo !j

! north flux
do j = 2, jmt
  ! western edge
  if ( DZT_k(1,j).ne.c0 ) then
    if ( DZT_k(imt,j).ne.c0 .or. DZT_k(2,j).ne.c0  ) then
      VN_NEW_FIELD(1,j) = 0.5 *&
       ( VN_FIELD(imt,j-1) * DXU(imt,j-1)&
       + VN_FIELD(1  ,j-1) * DXU(1  ,j-1) )&
       / DXT(1,j)
    else
      VN_NEW_FIELD(1,j) = 0.0
    endif
  endif
  ! rest of array
  do i = 2, imt
    if ( DZT_k(i,j).ne.c0 ) then
      if ( DZT_k(i-1,j).ne.c0 .or. DZT_k(i+1,j).ne.c0  ) then
        VN_NEW_FIELD(i,j) = 0.5 *&
         ( VN_FIELD(i-1,j-1) * DXU(i-1,j-1)&
         + VN_FIELD(i  ,j-1) * DXU(i  ,j-1) )&
         / DXT(i,j)
      else
        VN_NEW_FIELD(i,j) = 0.0
      endif
    endif
  enddo !i
enddo !j

!call surf_int(1,1,imt,jmt,UE_FIELD,    UAREA,DZT_k,A)
!call surf_int(1,1,imt,jmt,UE_NEW_FIELD,TAREA,DZT_k,B)
!      write (*,*) 'A: ', A, 'B: ', B
!write (*,*) k, 'rel. diff. of UE flux interp.: ', (A-B)/B
!call surf_int(1,1,imt,jmt,VN_FIELD,    UAREA,DZT_k,A)
!call surf_int(1,1,imt,jmt,VN_NEW_FIELD,TAREA,DZT_k,B)
!      write (*,*) 'A: ', A, 'B: ', B
!write (*,*) k, 'rel. diff. of VN flux interp.: ', (A-B)/B

end subroutine interp_mom_fluxes
!===============================================================================
