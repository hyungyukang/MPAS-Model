! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_hadv
!
!> \brief MPAS ocean horizontal advection for thickness
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies for thickness from horizontal advection
!
!-----------------------------------------------------------------------

module ocn_thick_hadv

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_mesh
   use ocn_config

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_thick_hadv_tend, &
             ocn_thick_hadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: thickHadvOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_hadv_tend
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details
!>  This routine computes the horizontal advection tendency for
!>  thicknes based on current state and user choices of forcings.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_hadv_tend(meshPool, normalVelocity, &
      layerThicknessEdge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity    !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: thickness at edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2, k, i, iCell, nCells
      real (kind=RKIND) :: flux, invAreaCell, invAreaCell1, invAreaCell2
      real (kind=RKIND) :: temp2

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.thickHadvOn) return

      call mpas_timer_start("thick hadv")

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(areaCell,nEdgesOnCell,edgesOnCell, &
      !$acc        minLevelEdgeBot,maxLevelEdgeTop,normalVelocity,dvEdge, &
      !$acc        layerThicknessEdge,tend,edgeSignOnCell) &
      !$acc    collapse(2) 
      do iCell = 1, nCellsOwned
        do k = 1,nVertLevels
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            flux = normalVelocity(k, iEdge) * dvEdge(iEdge) &
                 * layerThicknessEdge(k, iEdge)
            temp2 = edgeSignOnCell(i, iCell) * flux / areaCell(iCell)
            !!$acc atomic update
            tend(k,iCell) = tend(k,iCell) + temp2 
          end do
        end do
      end do

#else
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell, i, iEdge, k, flux)
      do iCell = 1, nCellsOwned
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
            flux = normalVelocity(k, iEdge) * dvEdge(iEdge) &
                 * layerThicknessEdge(k, iEdge)
            tend(k, iCell) = tend(k, iCell) &
                           + edgeSignOnCell(i, iCell) &
                           * flux * invAreaCell
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("thick hadv")

   !--------------------------------------------------------------------

   end subroutine ocn_thick_hadv_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_hadv_init
!
!> \brief   Initializes ocean horizontal thickness advection
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \details
!>  This routine initializes quantities related to horizontal thickness
!>  advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_hadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      thickHadvOn = .true.

      if(config_disable_thick_hadv) thickHadvOn = .false.

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_thick_hadv_init!}}}

!***********************************************************************

end module ocn_thick_hadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
