! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_ale
!
!> \brief MPAS ocean ALE thickness driver
!> \author Mark Petersen
!> \date   August 2013
!> \details
!>  This module contains the routines for computing ALE thickness.
!
!-----------------------------------------------------------------------

module ocn_thick_ale

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_mesh

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

   public :: ocn_ALE_thickness, &
             ocn_thick_ale_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: configALEthicknessProportionality

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_ALE_thickness
!
!> \brief   Computes desired ALE thickness at new time
!> \author  Mark Petersen
!> \date    August 2013
!> \details
!>  This routine computes the desired Arbitrary Lagrangian-Eulerian (ALE)
!>  thickness at the new time. It uses the ALE formulation, and includes
!>  contributions from temp1Cell variations (z-star), high-frequency divergence
!>  (z-tilde), and imposes a minimum layer thickness.
!
!-----------------------------------------------------------------------
   !subroutine ocn_ALE_thickness(meshPool, verticalMeshPool, temp1Cell, temp2VertCell, err, newHighFreqThickness)!{{{
   subroutine ocn_ALE_thickness( &
                                restingThickness, err, &
                                newHighFreqThickness)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

!     type (mpas_pool_type), intent(in) :: &
!        meshPool           !< Input: horizonal mesh information

!     type (mpas_pool_type), intent(in) :: &
!        verticalMeshPool   !< Input: vertical mesh information

!     real (kind=RKIND), dimension(:), intent(in) :: &
!        temp1Cell     !< Input: sea surface height

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         newHighFreqThickness   !< Input: high frequency thickness.  Alters ALE thickness.

      real (kind=RKIND), dimension(:,:), intent(in), pointer :: &
         restingThickness   !>  Layer thickness when the ocean is at rest, i.e. without temp1Cell or internal perturbations.

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

!     real (kind=RKIND), dimension(:,:), intent(out) :: &
!        temp2VertCell     !< Output: desired thickness at new time

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, i, kMax, kMin
      integer :: nCells

      real (kind=RKIND) :: weightSum, thicknessSum, remainder, newThickness, thicknessWithRemainder

      logical, pointer :: thicknessFilterActive

      err = 0

      call mpas_pool_get_package(ocnPackages, 'thicknessFilterActive', thicknessFilterActive)


      ! Initialization of temporary arrays

#ifdef MPAS_OPENACC
      !$acc parallel loop present(temp2Cell)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1,nCellsAll
         temp2Cell(iCell) = 1e-14_RKIND
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsHalo( 1 )
 
      !
      ! ALE thickness alteration due to temp1Cell (z-star)
      !
      if (configALEthicknessProportionality==1) then ! restingThickness_times_weights

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,temp2Cell)  &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               temp2Cell(iCell) = temp2Cell(iCell) &
                                + vertCoordMovementWeights(k) * restingThickness(k, iCell)
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif


#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,temp2VertCell,temp1Cell,temp2Cell) &
         !$acc    private(k,thicknessSum) collapse(2)
         do iCell = 1, nCells

            ! Note that restingThickness is nonzero, and remaining terms are perturbations about zero.
            ! This is equation 4 and 6 in Petersen et al 2015, but with eqn 6
            do k = 1,nVertLevels
            thicknessSum = temp1Cell(iCell) / temp2Cell(iCell)
               temp2VertCell(k, iCell) = restingThickness(k, iCell) &
                  + ( vertCoordMovementWeights(k) * restingThickness(k, iCell) * thicknessSum )
            end do
         enddo
#else
         !$omp parallel
         !$omp do schedule(runtime) private(kMax, thicknessSum, k)
         do iCell = 1, nCells
            thicknessSum = temp1Cell(iCell) / temp2Cell(iCell)

            ! Note that restingThickness is nonzero, and remaining terms are perturbations about zero.
            ! This is equation 4 and 6 in Petersen et al 2015, but with eqn 6
            !do k = kMin, kMax
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               temp2VertCell(k, iCell) = restingThickness(k, iCell) &
                  + ( vertCoordMovementWeights(k) * restingThickness(k, iCell) * thicknessSum )
            end do
         enddo
         !$omp end do
         !$omp end parallel
#endif


         elseif (configALEthicknessProportionality==2) then ! weights_only

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,temp2Cell) &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               temp2Cell(iCell) = temp2Cell(iCell) + vertCoordMovementWeights(k)
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,temp2VertCell,temp1Cell,temp2Cell) &
         !$acc    private(weightSum,k) collapse(2)
         do iCell = 1, nCells
            do k = 1,nVertLevels
               ! Using this, we must require that sum(restingThickness(k, iCell))
               ! summed over k is equal to bottomDepth.
               ! This is equation 4 and 6 in Petersen et al 2015, but with eqn 6
               ! altered so only the W_k weights are used. The resting
               ! thickness shown in eqn 6 is not included here.
               temp2VertCell(k, iCell) = restingThickness(k, iCell) + temp1Cell(iCell) &
                                       * vertCoordMovementWeights(k) / temp2Cell(iCell)
            end do
         enddo
#else
         !$omp parallel
         !$omp do schedule(runtime) private(kMax, weightSum, k)
         do iCell = 1, nCells
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! Using this, we must require that sum(restingThickness(k, iCell))
               ! summed over k is equal to bottomDepth.
               ! This is equation 4 and 6 in Petersen et al 2015, but with eqn 6
               ! altered so only the W_k weights are used. The resting
               ! thickness shown in eqn 6 is not included here.
               temp2VertCell(k, iCell) = restingThickness(k, iCell) + temp1Cell(iCell) &
                                       * vertCoordMovementWeights(k) / temp2Cell(iCell)
            end do
         enddo
         !$omp end do
         !$omp end parallel
#endif
     endif

      if (thicknessFilterActive) then
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelCell,maxLevelCell,temp2VertCell,newHighFreqThickness) &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(kMax)
#endif
         do iCell = 1, nCells
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               temp2VertCell(k, iCell) = &
                   temp2VertCell(k, iCell) + newHighFreqThickness(k,iCell)
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      !
      ! ALE thickness alteration due to minimum and maximum thickness
      !
      if (config_use_min_max_thickness) then

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelCell,maxLevelCell,temp1Vert, &
         !$acc            temp2VertCell,restingThickness, &
         !$acc            config_min_thickness,config_max_thickness_factor) &
         !$acc private(kMin,kMax, temp1Vert, remainder, k, newThickness, &
         !$acc         temp2Vert,temp3Vert)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(kMax, temp1Vert, remainder, k, newThickness, &
         !$omp         temp2Vert,temp3Vert)
#endif
         do iCell = 1, nCells
            kMax = maxLevelCell(iCell)
            kMin = minLevelCell(iCell)

            ! go down the column:
 
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               temp1Vert(k) = temp2VertCell(k, iCell)
            end do

            remainder = 0.0_RKIND
            !$acc loop reduction(-:remainder)
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               newThickness = max( min(temp1Vert(k) + remainder, &
                                      config_max_thickness_factor * restingThickness(k,iCell) ), &
                                  config_min_thickness)
               temp2Vert(k) = newThickness - temp1Vert(k)
               remainder = remainder - temp2Vert(k)
            end do

            ! go back up the column:
            temp3Vert(kMax) = 0.0_RKIND
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               temp1Vert(k) = temp1Vert(k) + temp2Vert(k)
            end do

            !do k = kMax-1, kMin, -1
            do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
               newThickness = max( min(temp1Vert(k) + remainder, &
                                      config_max_thickness_factor * restingThickness(k,iCell) ), &
                                  config_min_thickness)
               temp3Vert(k) = newThickness - temp1Vert(k)
               remainder = remainder - temp3Vert(k)
            end do
            temp3Vert(kMin) = temp3Vert(kMin) + remainder

            do k = minLevelCell(iCell),maxLevelCell(iCell)
               temp2VertCell(k, iCell) = temp2VertCell(k, iCell) &
                                            + temp2Vert(k) &
                                            + temp3Vert(k)
            end do

         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif ! config_use_min_max_thickness

   end subroutine ocn_ALE_thickness!}}}

!***********************************************************************
!
!  routine ocn_thick_ale_init
!
!> \brief   Initializes flags used within diagnostics routines.
!> \author  Mark Petersen
!> \date    August 2013
!> \details
!>  This routine initializes flags related to quantities computed within
!>  other diagnostics routines.
!
!-----------------------------------------------------------------------
   subroutine ocn_thick_ale_init(err)!{{{
      integer, intent(out) :: err !< Output: Error flag

      if (config_ALE_thickness_proportionality=='restingThickness_times_weights') then
          configALEthicknessProportionality = 1
      elseif (config_ALE_thickness_proportionality=='weights_only') then
          configALEthicknessProportionality = 2
      else
          call mpas_log_write( &
             ' Warning: config_ALE_thickness_proportionality is not valid', &
             MPAS_LOG_CRIT)
      endif

      err = 0

    end subroutine ocn_thick_ale_init!}}}

end module ocn_thick_ale

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
