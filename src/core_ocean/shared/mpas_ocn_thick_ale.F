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
!>  contributions from SSH variations (z-star), high-frequency divergence
!>  (z-tilde), and imposes a minimum layer thickness.
!
!-----------------------------------------------------------------------
   subroutine ocn_ALE_thickness(meshPool, verticalMeshPool, SSH, ALE_thickness, err, newHighFreqThickness)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool           !< Input: horizonal mesh information

      type (mpas_pool_type), intent(in) :: &
         verticalMeshPool   !< Input: vertical mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
         SSH     !< Input: sea surface height

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         newHighFreqThickness   !< Input: high frequency thickness.  Alters ALE thickness.

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         ALE_thickness     !< Output: desired thickness at new time

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, i, kMax, kMin
      !integer, pointer :: nVertLevels
      integer :: nCells
      !integer, dimension(:), pointer :: maxLevelCell
      !integer, dimension(:), pointer :: minLevelCell
      !integer, dimension(:), pointer :: nCellsArray

      real (kind=RKIND) :: weightSum, thicknessSum, remainder, newThickness, thicknessWithRemainder
      !real (kind=RKIND), dimension(:), pointer :: vertCoordMovementWeights
      real (kind=RKIND), dimension(:), allocatable :: &
         prelim_ALE_thickness,   & !> ALE thickness at new time
         min_ALE_thickness_down, & !> ALE thickness alteration due to min/max thickness
         min_ALE_thickness_up, &   !> ALE thickness alteration due to min/max thickness
         thick_weight_sum
      real (kind=RKIND), dimension(:,:), pointer :: &
         restingThickness   !>  Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations.

      logical, pointer :: thicknessFilterActive

      err = 0

      call mpas_pool_get_package(ocnPackages, 'thicknessFilterActive', thicknessFilterActive)

      !call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      !call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      !call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)

      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

      !call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      !call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      allocate(prelim_ALE_thickness(nVertLevels),thick_weight_sum(nCellsAll), &
         min_ALE_thickness_down(nVertLevels), min_ALE_thickness_up(nVertLevels))

      !$acc enter data &
      !$acc       copyin(restingThickness,config_max_thickness_factor, &
      !$acc              config_min_thickness,config_use_min_max_thickness, &
      !$acc              configALEthicknessProportionality) &
      !$acc       create(prelim_ALE_thickness,thick_weight_sum, &
      !$acc              min_ALE_thickness_down,min_ALE_thickness_up)
                       

#ifdef MPAS_OPENACC
      !$acc parallel loop present(thick_weight_sum)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1,nCellsAll
         thick_weight_sum(iCell) = 1e-14_RKIND
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsHalo( 1 )
!
!      !
!      ! ALE thickness alteration due to SSH (z-star)
!      !
      if (configALEthicknessProportionality==1) then ! restingThickness_times_weights

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,thick_weight_sum)  &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               thick_weight_sum(iCell) = thick_weight_sum(iCell) &
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
         !$acc            restingThickness,ALE_thickness,SSH,thick_weight_sum) &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(kMax, thicknessSum, k)
#endif
         do iCell = 1, nCells
            !kMax = maxLevelCell(iCell)
            !kMin = minLevelCell(iCell)
   
            !thicknessSum = 1e-14_RKIND
            !do k = kMin, kMax
            !!$acc loop reduction(+:thicknessSum)
!           do k = minLevelCell(iCell), maxLevelCell(iCell)
!              !thicknessSum = thicknessSum + vertCoordMovementWeights(k) * restingThickness(k, iCell)
!              thick_weight_sum(iCell) = thick_weight_sum(iCell) &
!                                      + vertCoordMovementWeights(k) * restingThickness(k, iCell)
!           end do
   
            ! Note that restingThickness is nonzero, and remaining terms are perturbations about zero.
            ! This is equation 4 and 6 in Petersen et al 2015, but with eqn 6
            !do k = kMin, kMax
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ALE_thickness(k, iCell) = restingThickness(k, iCell) &
                  + ( SSH(iCell) * vertCoordMovementWeights(k) * restingThickness(k, iCell) ) &
                    / thick_weight_sum(iCell)
!                   / thicknessSum
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif


         elseif (configALEthicknessProportionality==2) then ! weights_only

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,thick_weight_sum) &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               thick_weight_sum(iCell) = thick_weight_sum(iCell) + vertCoordMovementWeights(k)
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(maxLevelCell,minLevelCell,vertCoordMovementWeights, &
         !$acc            restingThickness,ALE_thickness,SSH,thick_weight_sum) &
         !$acc    private(weightSum,k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(kMax, weightSum, k)
#endif
         do iCell = 1, nCells
            !kMax = maxLevelCell(iCell)
            !kMin = minLevelCell(iCell)
   
!           weightSum = 1e-14_RKIND
!           !do k = kMin, kMax
!           do k = minLevelCell(iCell), maxLevelCell(iCell)
!              weightSum = weightSum + vertCoordMovementWeights(k) 
!           end do
   
            !do k = kMin, kMax
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! Using this, we must require that sum(restingThickness(k, iCell))
               ! summed over k is equal to bottomDepth.
               ! This is equation 4 and 6 in Petersen et al 2015, but with eqn 6
               ! altered so only the W_k weights are used. The resting
               ! thickness shown in eqn 6 is not included here.
               ALE_thickness(k, iCell) = restingThickness(k, iCell) + ssh(iCell) * vertCoordMovementWeights(k) &
                                       / thick_weight_sum(iCell)
!                                      / weightSum
   
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
     endif

      if (thicknessFilterActive) then
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelCell,maxLevelCell,ALE_thickness,newHighFreqThickness) &
         !$acc    private(k)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(kMax)
#endif
         do iCell = 1, nCells
            !kMax = maxLevelCell(iCell)
            !kMin = minLevelCell(iCell)
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               ALE_thickness(k, iCell) = &
                   ALE_thickness(k, iCell) + newHighFreqThickness(k,iCell)
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
         !$acc    present(minLevelCell,maxLevelCell,prelim_ALE_thickness, &
         !$acc            ALE_thickness,restingThickness, &
         !$acc            config_min_thickness,config_max_thickness_factor) &
         !$acc private(kMin,kMax, prelim_ALE_Thickness, remainder, k, newThickness, &
         !$acc         min_ALE_thickness_down, min_ALE_thickness_up)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(kMax, prelim_ALE_Thickness, remainder, k, newThickness, &
         !$omp         min_ALE_thickness_down, min_ALE_thickness_up)
#endif
         do iCell = 1, nCells
            kMax = maxLevelCell(iCell)
            kMin = minLevelCell(iCell)

            ! go down the column:
 
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               prelim_ALE_thickness(k) = ALE_thickness(k, iCell)
            end do

            remainder = 0.0_RKIND
            !$acc loop reduction(-:remainder)
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               newThickness = max( min(prelim_ALE_thickness(k) + remainder, &
                                      config_max_thickness_factor * restingThickness(k,iCell) ), &
                                  config_min_thickness)
               min_ALE_thickness_down(k) = newThickness - prelim_ALE_thickness(k)
               remainder = remainder - min_ALE_thickness_down(k)
            end do

            ! go back up the column:
            min_ALE_thickness_up(kMax) = 0.0_RKIND
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               prelim_ALE_thickness(k) = prelim_ALE_thickness(k) + min_ALE_thickness_down(k)
            end do

            !do k = kMax-1, kMin, -1
            do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
               newThickness = max( min(prelim_ALE_thickness(k) + remainder, &
                                      config_max_thickness_factor * restingThickness(k,iCell) ), &
                                  config_min_thickness)
               min_ALE_thickness_up(k) = newThickness - prelim_ALE_thickness(k)
               remainder = remainder - min_ALE_thickness_up(k)
            end do
            min_ALE_thickness_up(kMin) = min_ALE_thickness_up(kMin) + remainder

            do k = minLevelCell(iCell),maxLevelCell(iCell)
               ALE_thickness(k, iCell) = ALE_thickness(k, iCell) &
                                            + min_ALE_thickness_down(k) &
                                            + min_ALE_thickness_up(k)
            end do

         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif ! config_use_min_max_thickness

      !$acc exit data &
      !$acc       delete(restingThickness,prelim_ALE_thickness,thick_weight_sum, &
      !$acc              min_ALE_thickness_down,min_ALE_thickness_up,&
      !$acc              config_max_thickness_factor,config_min_thickness, &
      !$acc              config_use_min_max_thickness, &
      !$acc              configALEthicknessProportionality)
      deallocate(prelim_ALE_thickness, min_ALE_thickness_down, min_ALE_thickness_up, &
                 thick_weight_sum)

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
