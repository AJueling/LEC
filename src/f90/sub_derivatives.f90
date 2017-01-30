!===============================================================================
subroutine nabla_hvel(imt,jmt,km,VEL,TTT_VEL,dz,DXT,DYT,DXU,DYU,DZT,DZU,&
                            TAREA,DVELDX,DVELDY,DVELDZ)
!
!     calculates 3D gradients of a horizontal velocity field
!
implicit none

! input/output variables
integer,                     intent(in)  :: imt,jmt,km
real, dimension(        km), intent(in)  :: dz
real, dimension(imt,jmt   ), intent(in)  :: DXT, DYT, DXU, DYU, TAREA
real, dimension(imt,jmt,km), intent(in)  :: DZT, DZU
real, dimension(imt,jmt,km), intent(in)  :: VEL, TTT_VEL
real, dimension(imt,jmt,km), intent(out) :: DVELDX, DVELDY, DVELDZ ! TTT
gradients
! internal variables
real,             parameter                         :: p5=0.5, c0=0.0

write (*,*) 'hdiv vel started'

! horizontal gradients, central difference,
! eq. 3.6 page 16 in 2010 POP Reference Manual
do k = 1,km
! southernmost gridpoints (j=1) on TTT-grid must be 0,
! as there would not be UU gridpoints south of it otherwise
  do j = 2,jmt
    ! westernmost gridpoints (i=1)
    if ( DZT(1,j,k).ne.c0 ) then
      DVELDX(1,j,k) = p5*( ( DYU(1  ,j  )*DZU(1  ,j  ,k)*VEL(1  ,j  ,k)     &
                           + DYU(1  ,j-1)*DZU(1  ,j-1,k)*VEL(1  ,j-1,k) )   &
                         - ( DYU(imt,j  )*DZU(imt,j  ,k)*VEL(imt,j  ,k)     &
                           + DYU(imt,j-1)*DZU(imt,j-1,k)*VEL(imt,j-1,k) ) ) &
                         / TAREA(1,j) / DZT(1,j,k)
      DVELDY(1,j,k) = p5*( ( DXU(imt,j  )*DZU(imt,j  ,k)*VEL(imt,j  ,k)     &
                           + DXU(1  ,j  )*DZU(1  ,j  ,k)*VEL(1  ,j  ,k) )   &
                         - ( DXU(imt,j-1)*DZU(imt,j-1,k)*VEL(imt,j-1,k)     &
                           + DXU(1  ,j-1)*DZU(1  ,j-1,k)*VEL(1  ,j-1,k) ) ) &
                         / TAREA(1,j) / DZT(1,j,k)
    endif
    ! all gridpoints except western- and southernmost (i=1, j=1)
    do i = 2,imt
      if ( DZT(i,j,k).ne.c0 ) then
        DVELDX(i,j,k) = p5*( ( DYU(i  ,j  )*DZU(i  ,j  ,k)*VEL(i  ,j,k)    &
                             + DYU(i  ,j-1)*DZU(i  ,j-1,k)*VEL(i,j-1,k) )  &
                           - ( DYU(i-1,j  )*DZU(i-1,j  ,k)*VEL(i-1,j,k)    &
                             +
DYU(i-1,j-1)*DZU(i-1,j-1,k)*VEL(i-1,j-1,k) )) &
                           / TAREA(i,j) / DZT(i,j,k)
        DVELDY(i,j,k) = p5*( ( DXU(i-1,j  )*DZU(i-1,j  ,k)*VEL(i-1,j,k)    &
                             + DXU(i  ,j  )*DZU(i  ,j  ,k)*VEL(i  ,j,k) )  &
                           - ( DXU(i-1,j-1)*DZU(i-1,j-1,k)*VEL(i-1,j-1,k)    &
                             + DXU(i  ,j-1)*DZU(i  ,j-1,k)*VEL(i,j-1,k) )) &
                           / TAREA(i,j) / DZT(i,j,k)
      endif
    enddo !i
  enddo !j
enddo !k

! vertical gradient
! top/bottom layer
do j = 1,jmt
  do i = 1,imt
    if ( DZT(i,j, 1).ne.c0 .and. DZT(i,j,2).ne.c0 ) then
      DVELDZ(i,j, 1) = 2*( TTT_VEL(i,j,1) - TTT_VEL(i,j,2) )&
                        / ( DZT(i,j,1) + DZT(i,j,2) )
    endif
    if ( DZT(i,j,km).ne.c0 ) then
      DVELDZ(i,j,km) = 2*( TTT_VEL(i,j,km-1) - TTT_VEL(i,j,km) )&
                        / ( DZT(i,j,km-1) + DZT(i,j,km) )
    endif
  enddo !i 
enddo !j

do k = 2,km-1
  do j = 1,jmt
    do i = 1,imt
      if ( DZT(i,j,k).ne.c0 .and. DZT(i,j,k+1).ne.c0 ) then
        DVELDZ(i,j,k) = ( ( TTT_VEL(i,j,k-1) - TTT_VEL(i,j,k) )&
                          / ( DZT(i,j,k-1) + DZT(i,j,k) )&
                        + ( TTT_VEL(i,j,k) - TTT_VEL(i,j,k+1) )&
                          / (DZT(i,j,k) + DZT(i,j,k+1) ) )
      elseif ( DZT(i,j,k).ne.c0 .and. DZT(i,j,k+1).eq.c0 ) then
        DVELDZ(i,j,k) = 2*( TTT_VEL(i,j,k-1) - TTT_VEL(i,j,k) )&
                          / ( DZT(i,j,k-1) + DZT(i,j,k) )
      endif
    enddo !i
  enddo !j
enddo !k

end subroutine nabla_hvel
!===============================================================================

!===============================================================================
subroutine grad_rho(imt,jmt,km,kmT,DZT,DZU,DXU,DYU,TAREA,UAREA,RHO,      &
                          TT_DRHODX,TT_DRHODY)
!   call grad_rho(imt,jmt,km,DZT,DZU,DXU,DYU,TAREA,UAREA,PD,DPDDX,DPDDY)
!
!     calculates gradient components of RHO/PD
!     as described in POP reference manual
!     RHO/PD on 3111 grid, output on 3221 grid
!     units [g/cm^4]
!
implicit none

! input/output variables
integer,                        intent(in)  :: imt,jmt,km
integer, dimension(imt,jmt   ), intent(in)  :: kmT
real,    dimension(imt,jmt   ), intent(in)  :: DXU, DYU ! [m]
real,    dimension(imt,jmt   ), intent(in)  :: TAREA, UAREA ! [m^2]
real,    dimension(imt,jmt,km), intent(in)  :: DZT, DZU, RHO ! [m], [g/cm^3]
real,    dimension(imt,jmt,km)              :: DRHODX, DRHODY
real,    dimension(imt,jmt,km), intent(out) :: TT_DRHODX, TT_DRHODY
! internal variables
real,    parameter                          :: p5=0.5, c0=0.0


write (*,*) 'grad rho started'

DRHODX = c0
DRHODY = c0

do k = 1,km

  ! northermost line between poles of tripolar grid
  do i = 1,imt-1
    if ( DZU(i,jmt,k).ne.c0 ) then
      if ( i.lt.imt/2 ) then      ! first half of array, before wrapping around
        DRHODX(i,jmt,k) = ( ( RHO(i+1    ,jmt,k) + RHO(imt-i  ,jmt,k) )       &
                          - ( RHO(i      ,jmt,k) + RHO(imt-i+1,jmt,k) ) )     &
                          / 2 / DXU(i,jmt)
        DRHODY(i,jmt,k) = ( ( RHO(imt-i+1,jmt,k) + RHO(imt-i  ,jmt,k) )       &
                          - ( RHO(i      ,jmt,k) + RHO(i+1    ,jmt,k) ) )     &
                          / 2 / DYU(i,jmt)
      elseif ( i.gt.imt/2 ) then  ! second half, after wrap
        DRHODX(i,jmt,k) = DRHODX(imt-i,jmt,k)
        DRHODY(i,jmt,k) = DRHODY(imt-i,jmt,k)
      endif
    endif
  enddo !i

  do j = 1,jmt-1

  ! eastern boundary
    if ( DZU(imt,j,k).ne.c0 ) then
      DRHODX(imt,j,k) = ( ( RHO(1  ,j  ,k) + RHO(1  ,j+1,k) )               &
                        - ( RHO(imt,j  ,k) + RHO(imt,j+1,k) ) )             &
                        / 2 / DXU(imt,j)
      DRHODY(imt,j,k) = ( ( RHO(imt,j+1,k) + RHO(1  ,j+1,k) )               &
                        - ( RHO(imt,j  ,k) + RHO(1  ,j  ,k) ) )             &
                        / 2 / DYU(imt,j)
    endif
  ! rest of arrays
    do i = 1,imt-1
      if ( DZU(i,j,k).ne.c0 ) then
        DRHODX(i,j,k) = ( ( RHO(i+1,j  ,k) + RHO(i+1,j+1,k) )                &
                        - ( RHO(i  ,j  ,k) + RHO(i  ,j+1,k) ) )              &
                        / 2 / DXU(i,j)
        DRHODY(i,j,k) = ( ( RHO(i  ,j+1,k) + RHO(i+1,j+1,k) )                &
                        - ( RHO(i  ,j  ,k) + RHO(i+1,j  ,k) ) )              &
                        / 2 / DYU(i,j)
      endif
    enddo !i
  enddo !j

  ! Testing
  do i = 1,imt
    do j = 1,jmt
      if ( ISNAN(DRHODX(i,j,k)) ) then
        write (*,*) i,j,k, 'DRHODX is NaN'
      endif
      if ( ISNAN(DRHODY(i,j,k)) ) then
        write (*,*) i,j,k, 'DRHODY is NaN'
        if ( i.gt.3595 ) then
          write (*,*) DRHODX(i,j,k), DRHODY(i,j,k), kmT(i,j)
          write (*,*) DZT(i,j,k), DZU(i,j,k), DYU(i,j)
          write (*,*) RHO(i,j+1,k), RHO(i+1,j+1,k)
          write (*,*) RHO(i,j,k), RHO(i+1,j,k)
          write (*,*) ' '
          !DRHODY(i,j,k) = c0
        endif
      endif
      if ( DRHODY(i,j,k) == 1.0/0 .or. DRHODY(i,j,k) == -1/0.0 ) then
        write (*,*) i,j,k, 'DRHODY is +/- Infinity'
        DRHODY(i,j,k) = c0
      endif
    enddo !j
  enddo !i

  call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,&
                     DRHODX(:,:,k),TT_DRHODX(:,:,k))
  call uu2tt_scalar(imt,jmt,k,DZT(:,:,k),DZU(:,:,k),TAREA,UAREA,&
                     DRHODY(:,:,k),TT_DRHODY(:,:,k))

  ! Testing
  do j = 1,jmt
    do i = 1,imt
      if ( ISNAN(TT_DRHODX(i,j,k)) ) then
        write (*,*) i,j,k, 'TTT_DRHODX is NaN'
        write (*,*)  DZT(i,j,k), DZU(i,j,k), kmT(i,j), DXU(i,j)
        write (*,*)  DRHODX(i-1,j,k), DRHODX(i,j,k)
        write (*,*)  DRHODX(i-1,j-1,k), DRHODX(i,j-1,k)
        TT_DRHODX(i,j,k) = c0
      endif
      if ( ISNAN(TT_DRHODY(i,j,k)) ) then
        write (*,*) i,j,k, 'TTT_DRHODY is NaN'
        write (*,*)  DZT(i,j,k), DZU(i,j,k), kmT(i,j), DYU(i,j)
        write (*,*)  DRHODY(i-1,j,k), DRHODY(i,j,k)
        write (*,*)  DRHODY(i-1,j-1,k), DRHODY(i,j-1,k)
        TT_DRHODY(i,j,k) = c0
      endif
    enddo !i
  enddo !j

enddo !k

end subroutine grad_rho
!===============================================================================

!===============================================================================
subroutine vert_der(avg,dz,ddz)
! vertical derivative
integer,              parameter   :: km=42
integer                           :: k
real,   dimension(:), intent(in)  :: avg, dz
real,   dimension(:), intent(out) :: ddz

do k = 1, km
  if ( k==1 ) then
    ddz(k) = ( avg(k  ) - avg(k+1) ) / dz(k)
  else if ( k==km ) then
    ddz(k) = ( avg(k-1) - avg(k  ) ) / dz(k)
  else
    ddz(k) = ( avg(k-1) - avg(k+1) ) / 2.0 / dz(k)
  endif
enddo !k

end subroutine vert_der
!===============================================================================

