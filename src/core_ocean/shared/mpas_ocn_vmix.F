! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix
!
!> \brief MPAS ocean vertical mixing driver
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module is the main driver for
!>  vertical mixing in the ocean.
!>
!
!-----------------------------------------------------------------------

module ocn_vmix

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer

   use mpas_constants
   use ocn_constants
   use ocn_config
   use ocn_vmix_cvmix
   use ocn_vmix_coefs_redi
   use ocn_diagnostics_variables
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

   private :: tridiagonal_solve, &
              tridiagonal_solve_mult

   public :: ocn_vmix_coefs, &
             ocn_vel_vmix_tend_implicit, &
             ocn_tracer_vmix_tend_implicit, &
             ocn_vmix_init, &
             ocn_vmix_implicit, &
             ocn_compute_kpp_rhs

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: velVmixOn, tracerVmixOn
   real (kind=RKIND) :: implicitBottomDragCoef
   real (kind=RKIND) :: rayleighDampingCoef, rayleighBottomDampingCoef, &
                        rayleighDepthVariable

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs(meshPool, statePool, forcingPool, scratchPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: Scratch structure

      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool             !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool             !< Input/Output: forcing information

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

      integer :: err1, err2
      integer :: timeLevel

      integer :: k, iEdge, iCell, nEdges, nCells
!     integer, dimension(:), pointer :: nEdgesArray, nCellsArray

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing coefficients
      !
      !-----------------------------------------------------------------

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

!     call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
!     call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

!     nEdges = nEdgesArray( 2 )
      nEdges = nEdgesHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(vertViscTopOfEdge) collapse(2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1, nEdges
         do k = 1,nVertLevels
            vertViscTopOfEdge(k, iEdge) = 0.0_RKIND
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

!     nCells = nCellsArray( 2 )
      nCells = nCellsHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(vertDiffTopOfCell) collapse(2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells
         do k = 1,nVertLevels
            vertDiffTopOfCell(k, iCell) = 0.0_RKIND
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call ocn_vmix_coefs_cvmix_build(meshPool, statePool, forcingPool, err1, timeLevel)
      call ocn_vmix_coefs_redi_build(meshPool, statePool, err2, timeLevel)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_rayleigh
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!> \author  Mark Petersen, Phillip J. Wolfram
!> \date    September 2011, July 2019
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients and includes implicit rayleigh drag.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_rayleigh(meshPool, dt, kineticEnergyCell, vertViscTopOfEdge, layerThickness, & !{{{
                                         layerThicknessEdge, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp
      real (kind=RKIND) :: edgeThicknessTotal

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))
      A(1)=0.0_RKIND

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(N, cell1, cell2, edgeThicknessTotal, k, A, B, C, velTemp)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         edgeThicknessTotal = 0.0_RKIND
         do k = 1, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
            edgeThicknessTotal = edgeThicknessTotal + layerThicknessEdge(k,iEdge)
         end do

         ! A is lower diagonal term
         do k = 2, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1) = 1.0_RKIND - C(1) &
           ! added Rayleigh terms
           + dt*rayleighDampingCoef/((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal)
         do k = 2, N-1
            B(k) = 1.0_RKIND - A(k) - C(k) &
              ! added Rayleigh terms
              + dt*(rayleighDampingCoef/((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal))
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         B(N) = 1.0_RKIND - A(N) + dt*implicitBottomDragCoef &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge) &
              ! added Rayleigh terms
              + dt*(rayleighBottomDampingCoef + &
                    rayleighDampingCoef /((1.0_RKIND - rayleighDepthVariable) + rayleighDepthVariable*edgeThicknessTotal))

!        call tridiagonal_solve(A(2:N),B,C(1:N-1),normalVelocity(:,iEdge),velTemp,N)

         normalVelocity(1:N,iEdge) = velTemp(1:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do
      !$omp end parallel

      deallocate(A,B,C,velTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_rayleigh!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit(dt, layerThickness, & !{{{
                                         normalVelocity, err)
!                                        layerThicknessEdge, normalVelocity, err)

!     !$acc routine(tridiagonal_solve)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

!     type (mpas_pool_type), intent(in) :: &
!        meshPool          !< Input: mesh information

!     real (kind=RKIND), dimension(:,:), intent(in) :: &
!        kineticEnergyCell        !< Input: kinetic energy at cell

!     real (kind=RKIND), dimension(:,:), intent(in) :: &
!        vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

!     real (kind=RKIND), dimension(:,:), intent(inout) :: &
!        layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, Nsurf, nEdges,nCells,iCell
!     integer, pointer :: nVertLevels
!     integer, dimension(:), pointer :: nEdgesArray

!     integer, dimension(:), pointer :: maxLevelEdgeTop, minLevelEdgeBot

!     integer, dimension(:,:), pointer :: cellsOnEdge

!     real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp
      real (kind=RKIND) :: m

      err = 0

      !if(.not.velVmixOn) return
      return

!     call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
!     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

!     call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
!     call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
!     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      !nEdges = nEdgesArray( 1 )
      nEdges = nEdgesOwned
      nCells = nCellsOwned

!     allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))

#ifdef MPAS_OPENACC
      !$acc parallel loop &

!     !$acc    present(maxLevelEdgeTop,minLevelEdgeBot,cellsOnEdge,layerThickEdge, &
!     !$acc       layerThickness,temp1Vert,vertViscTopOfEdge,temp3Vert,temp2Vert, &
!     !$acc       kineticEnergyCell,temp4Vert,normalVelocity,implicitBottomDragCoef, &
!     !$acc       temp6Vert,temp7Vert) &
!     !$acc    private(N, cell1, cell2, temp1Vert, temp3Vert,temp4Vert,k, &
!     !$acc       temp6Vert,temp7Vert,temp2Vert)

      !$acc    present(maxLevelEdgeTop,minLevelEdgeBot, &
      !$acc       temp2Vert, &
      !$acc       normalVelocity)&
      !$acc    private(N, k, &
      !$acc       temp2Vert)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(N, cell1, cell2, temp1Vert, temp2Vert, temp3Vert,temp4Vert, temp6Vert,temp7Vert,k)
#endif
      do iEdge = 1, nEdges
!       N = maxLevelEdgeTop(iEdge)
!       Nsurf = minLevelEdgeBot(iEdge)
!       if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThickEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThickEdge here.
!        cell1 = cellsOnEdge(1,iEdge)
!        cell2 = cellsOnEdge(2,iEdge)
!        do k = Nsurf, N
!           layerThickEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
!        end do

         ! A is lower diagonal term
!        temp1Vert(1:Nsurf)=0.0_RKIND
!        do k = 1,Nsurf
!        do k = Nsurf,N
!           temp1Vert(k)=0.0_RKIND
!           temp2Vert(k)=0.0_RKIND
!           temp3Vert(k)=0.0_RKIND
!           temp4Vert(k)=0.0_RKIND
!        end do


!        do k = Nsurf+1, N
!           temp1Vert(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
!              / (layerThickEdge(k-1,iEdge) + layerThickEdge(k,iEdge)) &
!              / layerThickEdge(k,iEdge)
!        enddo



         ! C is upper diagonal term
!        temp3Vert(1:Nsurf-1)=0.0_RKIND
!        do k = 1,Nsurf-1
!           temp3Vert(k)=0.0_RKIND
!        end do

     

!        do k = Nsurf, N-1
!           temp3Vert(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
!              / (layerThickEdge(k,iEdge) + layerThickEdge(k+1,iEdge)) &
!              / layerThickEdge(k,iEdge)
!        enddo


         ! B is diagonal term
!        temp2Vert(1:Nsurf-1)=0.0_RKIND
!        do k = 1,Nsurf-1
!           temp2Vert(k)=0.0_RKIND
!        end do

         !temp2Vert(Nsurf) = 1.0_RKIND  - temp3Vert(Nsurf)
!        temp2Vert(Nsurf) = Nsurf


!        do k = Nsurf+1, N-1
!        do k = minLevelEdgeBot(iEdge)+1,maxLevelEdgeTop(iEdge)-1
!           !temp2Vert(k) = 1.0_RKIND - temp1Vert(k) - temp3Vert(k)
!           !temp2Vert(k) = implicitBottomDragCoef
!           temp2Vert(k) = real(maxLevelEdgeTop(iEdge))
!           normalVelocity(k,iEdge) = maxLevelEdgeTop(iEdge)
!        enddo
!        temp2Vert(N) = 0.0_RKIND

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
!        temp2Vert(N) = 1.0_RKIND - temp1Vert(N) + dt*implicitBottomDragCoef &
!             * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdge(N,iEdge)

!        !do k = Nsurf,N
!        do k = minLevelEdgeBot(iEdge)+1,maxLevelEdgeTop(iEdge)-1
!           !normalVelocity(k,iEdge) = temp2Vert(k)
!           normalVelocity(k,iEdge) = maxLevelEdgeTop(iEdge)
!        end do

         !-------------------------------------------------------

!        call tridiagonal_solve(temp1Vert(Nsurf+1:N),temp2Vert(Nsurf:N),temp3Vert(Nsurf:N-1), &
!             normalVelocity(Nsurf:N,iEdge),temp4Vert(Nsurf:N),N-Nsurf+1, &
!             temp6Vert(Nsurf:N),temp7Vert(Nsurf:N))

!        ! Use work variables for b and r
!        temp6Vert(Nsurf) = temp2Vert(Nsurf)
!        temp7Vert(Nsurf) = normalVelocity(Nsurf,iEdge)
!  
!        ! First pass: set the coefficients
!        do k = Nsurf+1,N
!           m = temp1Vert(k-1)/temp6Vert(k-1)
!           temp6Vert(k) = temp2Vert(k) - m*temp3Vert(k-1)
!           temp7Vert(k) = normalVelocity(k,iEdge) - m*temp7Vert(k-1)
!        end do
!  
!        temp4Vert(N) = temp7Vert(N)/temp6Vert(N)
!        ! Second pass: back-substition
!        do k = n-1, 1, -1
!           temp4Vert(k) = (temp7Vert(k) - temp3Vert(k)*temp4Vert(k+1))/temp6Vert(k)
!        end do

         !-------------------------------------------------------

!        normalVelocity(1:Nsurf-1,iEdge) = 0.0_RKIND
!        normalVelocity(Nsurf:N,iEdge) = temp4Vert(Nsurf:N)
!        normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

!        do k = 1,Nsurf-1
!           normalVelocity(k,iEdge) = 0.0_RKIND
!        end do

!        do k = Nsurf,N
!           normalVelocity(k,iEdge) = temp4Vert(k)
!        end do

!        do k = N+1,nVertLevels
!           normalVelocity(k,iEdge) = 0.0_RKIND
!        end do

!       end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc update host(normalVelocity)
      do iEdge = 1,nEdges
        N = maxLevelEdgeTop(iEdge)
        Nsurf = minLevelEdgeBot(iEdge)
        do k = Nsurf,N
           print*,k,iEdge,normalVelocity(k,iEdge)
        end do
      end do

!     do iCell = 1,nCells
!       N = maxLevelCell(iCell)
!       Nsurf = minLevelCell(iCell)
!       do k = Nsurf,N
!          print*, kineticEnergyCell(k,iCell)
!       end do
!     end do

      stop

!     deallocate(A,B,C,velTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_variable
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!>          with spatially-variable bottom drag
!> \author  Mark Petersen, Phillip J. Wolfram
!> \date    September 2011, September 2019
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients from spatially-variable bottom drag.
!>  Except for bottom drag coefficient, routine should be identifcal to
!>  ocn_vel_vmix_tend_implicit above.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_spatially_variable(meshPool, bottomDrag, dt, kineticEnergyCell, & !{{{
                                         vertViscTopOfEdge, layerThickness, &
                                         layerThicknessEdge, normalVelocity, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

       real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDrag !< Input: bottomDrag at cell centeres

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, nEdges
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: implicitCd
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      err = 0

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))
      A(1)=0.0_RKIND

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! average cell-based implicit bottom drag to edges
         implicitCd = 0.5_RKIND*(bottomDrag(cell1) + bottomDrag(cell2))

         ! A is lower diagonal term
         do k = 2, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1) = 1.0_RKIND - C(1)
         do k = 2, N-1
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         ! use implicitCd from spatially variable bottom drag
         B(N) = 1.0_RKIND - A(N) + dt*implicitCd &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge)

!        call tridiagonal_solve(A(2:N),B,C(1:N-1),normalVelocity(:,iEdge),velTemp,N)

         normalVelocity(1:N,iEdge) = velTemp(1:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_spatially_variable!}}}


!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit_variable_mannings
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!>          with spatially-variable bottom drag using Mannings friction
!> \author  Mark Petersen, Phillip J. Wolfram
!> \date    September 2011, September 2019
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients from spatially-variable bottom drag.
!>  Except for bottom drag coefficient, routine should be identifcal to
!>  ocn_vel_vmix_tend_implicit above.  Cd uses Mannings' n values for the
!>  Cd=g*n^2*h^{-1/3}.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit_spatially_variable_mannings(meshPool, forcingPool, bottomDrag, dt, & !{{{
                                         kineticEnergyCell, vertViscTopOfEdge, layerThickness, &
                                         layerThicknessEdge, normalVelocity, ssh, bottomDepth, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool          !< Input: forcing information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyCell        !< Input: kinetic energy at cell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: thickness at cell center

       real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDrag !< Input: bottomDrag at cell centeres

       real (kind=RKIND), dimension(:), pointer :: ssh

       real (kind=RKIND), dimension(:), pointer :: bottomDepth

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThicknessEdge        !< Input: thickness at edge

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

      integer :: iEdge, k, cell1, cell2, N, nEdges
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: implicitCd
      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, velTemp

      ! vegetation_drag
      real (kind=RKIND), dimension(:), pointer :: vegetationHeight
      real (kind=RKIND), dimension(:), pointer :: vegetationDiameter
      real (kind=RKIND), dimension(:), pointer ::vegetationDensity
      real (kind=RKIND), dimension(:), pointer ::vegetationManning
      integer, dimension(:), pointer ::vegetationMask
      real (kind=RKIND) :: old_bottom_Cd, lambda, beta, alpha, total_h
      real (kind=RKIND) :: inundation_depth, von_karman, cff1, cff2, cff3, cff4
      integer :: iCell, nCells
      integer, pointer :: nCellsSolve

      err = 0
      von_karman = 0.4_RKIND

      if(.not.velVmixOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(forcingPool, 'vegetationMask', vegetationMask)
      call mpas_pool_get_array(forcingPool, 'vegetationHeight', vegetationHeight)
      call mpas_pool_get_array(forcingPool, 'vegetationDensity', vegetationDensity)
      call mpas_pool_get_array(forcingPool, 'vegetationDiameter', vegetationDiameter)
      call mpas_pool_get_array(forcingPool, 'vegetationManning', vegetationManning)

      nEdges = nEdgesArray( 1 )
      nCells = nCellsSolve

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),velTemp(nVertLevels))
      A(1)=0.0_RKIND

      ! Compute bottomDrag (Manning roughness) induced by vegetation
      if (config_use_vegetation_drag .AND. config_use_vegetation_manning_equation) then
        do iCell = 1, nCells
          if (vegetationDensity(iCell) * vegetationHeight(iCell) * vegetationDiameter(iCell) .eq. 0.0_RKIND) then
            vegetationMask(iCell) = 0
          endif
          if (vegetationMask(iCell) .eq. 1) then
            total_h = bottomDepth(iCell) + ssh(iCell)
            old_bottom_Cd = gravity * bottomDrag(iCell)**2.0_RKIND * total_h**(1.0_RKIND/3.0_RKIND)
            inundation_depth = MIN(vegetationHeight(iCell), total_h)
            inundation_depth = MAX(inundation_depth, 1e-6)
            lambda = vegetationDiameter(iCell) * vegetationDensity(iCell)
            alpha = (config_vegetation_drag_coefficient*lambda/ &
                   (4.0_RKIND*von_karman**2.0_RKIND*inundation_depth**2.0_RKIND))**(1.0_RKIND/3.0_RKIND)
            beta = 0.5_RKIND*alpha*inundation_depth*(1.0_RKIND - EXP(-2.0_RKIND*alpha*inundation_depth)) &
                  / (1.0_RKIND - EXP(-alpha*inundation_depth))**2.0_RKIND
            cff1 = total_h**(2.0_RKIND/3.0_RKIND) &
                  * SQRT((0.5_RKIND*beta*lambda*config_vegetation_drag_coefficient*inundation_depth &
                  + old_bottom_Cd)/(gravity*total_h))
            cff2 = (alpha*inundation_depth)**2.0_RKIND/(1.0_RKIND - EXP(-alpha*inundation_depth))
            cff3 = (1.0_RKIND - EXP(-alpha*inundation_depth))/(alpha**2.0_RKIND*inundation_depth*total_h)
            cff4 = LOG(total_h/inundation_depth) - (1.0_RKIND-inundation_depth/total_h) &
                  * (1.0_RKIND - 1.0_RKIND/(alpha*inundation_depth))
            vegetationManning(iCell) = cff1/(cff2*(cff3+cff4))
            vegetationManning(iCell) = MAX(bottomDrag(iCell), vegetationManning(iCell))
          else
            vegetationManning(iCell) = bottomDrag(iCell)
          endif
        enddo
      endif

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        N = maxLevelEdgeTop(iEdge)
        if (N .gt. 0) then

         ! Compute A(k), B(k), C(k)
         ! layerThicknessEdge is computed in compute_solve_diag, and is not available yet,
         ! so recompute layerThicknessEdge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, N
            layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do

         ! average cell-based implicit bottom drag to edges and convert Mannings n to Cd
         if (config_use_vegetation_drag .AND. config_use_vegetation_manning_equation) then
           implicitCd = gravity*(0.5_RKIND*(vegetationManning(cell1) + vegetationManning(cell2)))**2.0 * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
         else
           implicitCd = gravity*(0.5_RKIND*(bottomDrag(cell1) + bottomDrag(cell2)))**2.0 * &
            (0.5_RKIND * (ssh(cell1) + ssh(cell2) + bottomDepth(cell1) + bottomDepth(cell2)))**(-1.0_RKIND/3.0_RKIND)
         endif

         ! A is lower diagonal term
         do k = 2, N
            A(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k,iEdge) &
               / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! C is upper diagonal term
         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (layerThicknessEdge(k,iEdge) + layerThicknessEdge(k+1,iEdge)) &
               / layerThicknessEdge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1) = 1.0_RKIND - C(1)
         do k = 2, N-1
            B(k) = 1.0_RKIND - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         ! second line uses sqrt(2.0*kineticEnergyEdge(k,iEdge))
         ! use implicitCd from spatially variable bottom drag
         B(N) = 1.0_RKIND - A(N) + dt*implicitCd &
              * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThicknessEdge(N,iEdge)

!        call tridiagonal_solve(A(2:N),B,C(1:N-1),normalVelocity(:,iEdge),velTemp,N)

         normalVelocity(1:N,iEdge) = velTemp(1:N)
         normalVelocity(N+1:nVertLevels,iEdge) = 0.0_RKIND

        end if
      end do
      !$omp end do

      deallocate(A,B,C,velTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit_spatially_variable_mannings!}}}


!***********************************************************************
!
!  routine ocn_tracer_vmix_tend_implicit
!
!> \brief   Computes tendencies for implicit tracer vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the tendencies for implicit vertical mixing for
!>  tracers using computed coefficients.
!
!-----------------------------------------------------------------------

!  subroutine ocn_tracer_vmix_tend_implicit(meshPool, dt, vertDiffTopOfCell, layerThickness, tracers, &
!                 vertNonLocalFlux, tracerGroupSurfaceFlux, config_cvmix_kpp_nonlocal_with_implicit_mix, &
!                 err)!{{{
   subroutine ocn_tracer_vmix_tend_implicit(dt, layerThickness, tracers, &
                  tracerGroupSurfaceFlux, &
                  err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

!     type (mpas_pool_type), intent(in) :: &
!        meshPool          !< Input: mesh information

!     real (kind=RKIND), dimension(:,:), intent(in) :: &
!        vertDiffTopOfCell !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, &             !< Input: thickness at cell center
         tracerGroupSurfaceFlux        !< Input: surface flux for tracers nonlocal computation

!     real (kind=RKIND), dimension(:,:,:), intent(in) :: &
!        vertNonLocalFlux             !non local flux at interfaces

!     logical, intent(in) :: config_cvmix_kpp_nonlocal_with_implicit_mix
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers        !< Input: tracers

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

      integer :: iCell, k, num_tracers, N, Nsurf, nCells,iTracer
!     integer, pointer :: nVertLevels
!     integer, dimension(:), pointer :: nCellsArray

!     integer, dimension(:), pointer :: maxLevelCell, minLevelCell

      real (kind=RKIND), dimension(:), allocatable :: A,B,C
      real (kind=RKIND), dimension(:,:), allocatable :: tracersTemp, rhs
      real (kind=RKIND) :: m

      err = 0

      !if(.not.tracerVmixOn) return
      return

!     call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
!     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      num_tracers = size(tracers, dim=1)

!     call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
!     call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

!     allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),tracersTemp(num_tracers,nVertLevels))
!     allocate(rhs(num_tracers,nVertLevels))
!     allocate(tracersTemp(num_tracers,nVertLevels))
!     allocate(rhs(num_tracers,nVertLevels))

!     !$acc enter data copyin(tracersTemp,rhs)

      !nCells = nCellsArray( 1 )
      nCells = nCellsOwned

      call mpas_timer_start('vmix tracers tend imp loop', .false.)
#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector&
      !$acc    present(minLevelCell,maxLevelCell,temp1Vert,vertDiffTopOfCell, &
      !$acc       layerThickness,temp3Vert,temp2Vert,temp4Vert,temp5Vert,tracers, &
      !$acc       temp6Vert,temp7Vert) &
      !$acc    private(N, temp1Vert, temp2Vert, temp3Vert,temp4Vert, &
      !$acc       temp5Vert, k, iTracer,temp6Vert,temp7Vert)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(N, temp1Vert, temp2Vert, temp3Vert,temp4Vert,temp5Vert,k, temp6Vert,temp7Vert,iTracer)
#endif
      do iCell = 1, nCells
         ! Compute A(k), B(k), C(k) for tracers
         Nsurf = minLevelCell(iCell)
         N = maxLevelCell(iCell)

         ! A is lower diagonal term
         temp1Vert(Nsurf)=0.0_RKIND
         do k = Nsurf+1, N
            temp1Vert(k) = -2.0_RKIND*dt*vertDiffTopOfCell(k,iCell) &
                 / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) / layerThickness(k,iCell)
         enddo

         ! C is upper diagonal term
         do k = Nsurf, N-1
            temp3Vert(k) = -2.0_RKIND*dt*vertDiffTopOfCell(k+1,iCell) &
                 / (layerThickness(k,iCell) + layerThickness(k+1,iCell)) / layerThickness(k,iCell)
         enddo
         temp3Vert(N) = 0.0_RKIND

         ! B is diagonal term
         do k = Nsurf, N
            temp2Vert(k) = 1.0_RKIND - temp1Vert(k) - temp3Vert(k)
         enddo
 

!        !$acc loop seq
         do iTracer = 1,num_tracers
!        if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) then
!           call ocn_compute_kpp_rhs(tracers(:,:,iCell), rhs(:,:), dt, N, num_tracers, &
!                            layerThickness(:,iCell), vertNonLocalFlux(:,:,iCell), &
!                            tracerGroupSurfaceFlux(:,iCell)) !mlc-arguments may need to be changed
!        else
!           rhs(:,:) = tracers(:,:,iCell)
!        endif

!           do k = 1,nVertLevels
!              temp4Vert(k) = tracers(iTracer,k,iCell)
!           end do

!        call tridiagonal_solve_mult(temp1Vert(Nsurf+1:N), temp2Vert(Nsurf:N), temp3Vert(Nsurf:N-1), rhs(:,Nsurf:N), &
!             tracersTemp(:,Nsurf:N), N-Nsurf+1, num_tracers)

!           call tridiagonal_solve(temp1Vert(Nsurf+1:N), temp2Vert(Nsurf:N), temp3Vert(Nsurf:N-1),temp4Vert(Nsurf:N), &
!                temp5Vert(Nsurf:N), N-Nsurf+1,temp6Vert(Nsurf:N),temp7Vert(Nsurf:N))
         !--------------------------------------------------------------------
         ! Use work variables for b and r
         temp6Vert(Nsurf) = temp2Vert(Nsurf)
         temp7Vert(Nsurf) = tracers(iTracer,Nsurf,iCell)
   
         ! First pass: set the coefficients
         do k = Nsurf+1,N
            m = temp1Vert(k-1)/temp6Vert(k-1)
            temp6Vert(k) = temp2Vert(k) - m*temp3Vert(k-1)
            temp7Vert(k) = tracers(iTracer,k,iCell) - m*temp7Vert(k-1)
         end do
   
         temp4Vert(N) = temp7Vert(N)/temp6Vert(N)
         ! Second pass: back-substition
         do k = n-1, 1, -1
            temp4Vert(k) = (temp7Vert(k) - temp3Vert(k)*temp4Vert(k+1))/temp6Vert(k)
         end do
         !--------------------------------------------------------------------

            do k = Nsurf,N
               tracers(iTracer,k,iCell) = temp4Vert(k)
            end do

         end do


!        do k = Nsurf,N
!           do iTracer = 1,num_tracers
!              tracers(iTracer,k,iCell) = tracersTemp(iTracer,k)
!           end do
!        end do

         do k = N+1,nVertLevels
            do iTracer = 1,num_tracers
               tracers(iTracer,k,iCell) = -1e34
            end do
         end do

         do k = 1,Nsurf-1
            do iTracer = 1,num_tracers
               tracers(iTracer,k,iCell) = -1e34
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif
      call mpas_timer_stop('vmix tracers tend imp loop')

!     !$acc exit data delete(tracersTemp,rhs)
!     deallocate(tracersTemp, rhs)

      !deallocate(A, B, C, tracersTemp, rhs)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_tend_implicit!}}}

!***********************************************************************
!
!  routine ocn_vmix_implicit
!
!> \brief   Driver for implicit vertical mixing
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine is a driver for handling implicit vertical mixing
!>  of both momentum and tracers for a block. It's intended to reduce
!>  redundant code.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, scratchPool, &
              normalVelocity,layerThickness,ssh,nonLocalSurfaceTracerFlux, &
              tracersGroup, activeTracersSurfaceFlux, err, timeLevelIn)!{{{
!  subroutine ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, scratchPool, &
!             err, timeLevelIn)!{{{
      real (kind=RKIND), intent(in) :: dt
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool
      type (mpas_pool_type), intent(inout) :: forcingPool
      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: Scratch structure
!     real (kind=RKIND), dimension(:), pointer :: ssh ! needed for depth-variable computation
!     real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness
!     real (kind=RKIND), dimension(:,:), pointer :: &
!          nonLocalSurfaceTracerFlux,activeTracersSurfaceFlux

      real (kind=RKIND), dimension(:), pointer, intent(in) :: ssh ! needed for depth-variable computation
      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: normalVelocity, layerThickness
      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: &
           nonLocalSurfaceTracerFlux,activeTracersSurfaceFlux
      real (kind=RKIND), dimension(:,:,:), pointer, intent(in) :: tracersGroup
      integer, intent(out) :: err
      integer, intent(in), optional :: timeLevelIn

      type (mpas_pool_type), pointer :: tracersPool, tracersSurfaceFluxPool

      integer :: iCell, timeLevel, k, cell1, cell2, iEdge, nCells, nEdges
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray
      real (kind=RKIND), dimension(:), pointer :: bottomDrag, bottomDepth   ! needed for depth-variable computation
      real (kind=RKIND), dimension(:,:), pointer :: tracerGroupSurfaceFlux
!     real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      real (kind=RKIND), dimension(:,:,:), allocatable :: nonLocalFluxTend
!     integer, dimension(:), pointer :: maxLevelCell, minLevelCell, maxLevelEdgeTop, minLevelEdgeTop
!     integer, dimension(:,:), pointer :: cellsOnEdge

      type (mpas_pool_iterator_type) :: groupItr

      character (len=StrKIND) :: modifiedGroupName
!     integer, pointer :: indexTempFlux, indexSaltFlux, nVertLevels
      integer ::  num_tracers, iTracer
      err = 0

      call mpas_timer_start('vmix imp')

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
!     call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
!     call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
!     call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

!     call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
!     call mpas_pool_get_array(meshPool, 'minLevelCell', maxLevelCell)
!     call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
!     call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
!     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
!     call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

!     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
!     call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
!     call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

!     call mpas_pool_get_array(forcingPool, 'bottomDrag', bottomDrag)

      call mpas_timer_start('vmix coefs', .false.)
      call ocn_vmix_coefs(meshPool, statePool, forcingPool, scratchPool, err, timeLevel)
      call mpas_timer_stop('vmix coefs')

      !nCells = nCellsArray(1)
      nCells = nCellsOwned
      num_tracers = size(tracersGroup,dim=1)

      ! if using CVMix, then viscosity has to be averaged from cell centers to cell edges
!---------------------
      if ( config_use_cvmix ) then
         nEdges = nEdgesArray( 1 )
         call mpas_timer_start('CVMix avg', .false.)
         !$omp parallel
         !$omp do schedule(runtime) private(cell1, cell2, k)
         do iEdge=1,nEdges
            vertViscTopOfEdge(:, iEdge) = 0.0_RKIND
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
               vertViscTopOfEdge(k,iEdge) = 0.5_RKIND*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1))
            end do
         end do
         !$omp end do
         !$omp end parallel
         call mpas_timer_stop('CVMix avg')
      endif
!---------------------
      !
      !  Implicit vertical solve for momentum
      !
      call mpas_timer_start('vmix solve momentum', .false.)
      if (config_use_implicit_bottom_drag_variable) then
!---------------------
        call ocn_vel_vmix_tend_implicit_spatially_variable(meshPool, bottomDrag, dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, err)
      else if (config_use_implicit_bottom_drag_variable_mannings) then
        ! update bottomDrag via Cd=g*n^2*h^(-1/3)
        call ocn_vel_vmix_tend_implicit_spatially_variable_mannings(meshPool, forcingPool, bottomDrag, &
          dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, &
          ssh, bottomDepth, err)
      else if (config_Rayleigh_friction.or. &
               config_Rayleigh_bottom_friction.or. &
               config_Rayleigh_damping_depth_variable) then
        call ocn_vel_vmix_tend_implicit_rayleigh(meshPool, dt, kineticEnergyCell, &
          vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, err)
!---------------------
      else
!       call ocn_vel_vmix_tend_implicit(meshPool, dt, kineticEnergyCell, &
!         vertViscTopOfEdge, layerThickness, layerThickEdge, normalVelocity, err)
        call ocn_vel_vmix_tend_implicit(dt, &
          layerThickness, normalVelocity, err)
      end if
      call mpas_timer_stop('vmix solve momentum')

      !
      !  Implicit vertical solve for all tracers
      !

      call mpas_timer_start('vmix solve tracers', .false.)
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
!           call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, timeLevel)
            ! store tracers
            if (trim(groupItr % memberName) == 'activeTracers') then
               if (config_compute_active_tracer_budgets) then
#ifdef MPAS_OPENACC
                  !$acc parallel loop &
                  !$acc    present(activeTracerVertMixTendency,tracersGroup) &
                  !$acc    collapse(3)
#else
                  !$omp parallel
                  !$omp do schedule(runtime) private(k,iTracer)
#endif
                  do iCell = 1, nCells
                     do k = 1, nVertLevels
                        do iTracer = 1,num_tracers
                           activeTracerVertMixTendency(iTracer,k,iCell)=tracersGroup(iTracer,k,iCell)
                        end do
                     end do
                  end do
#ifndef MPAS_OPENACC
                  !$omp end do
                  !$omp end parallel
#endif
               endif
            endif


            if ( associated(tracersGroup) ) then
!              if (trim(groupItr % memberName) == 'activeTracers') then
!                 call mpas_pool_get_array(tracersSurfaceFluxPool, 'nonLocalSurfaceTracerFlux', &
!                          tracerGroupSurfaceFlux)
!              else
!                 modifiedGroupName = trim(groupItr % memberName) // "SurfaceFlux"
!                 call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), &
!                         tracerGroupSurfaceFlux)
!              endif

!              call ocn_tracer_vmix_tend_implicit(meshPool, dt, vertDiffTopOfCell, layerThickness, tracersGroup, &
!                       vertNonLocalFlux, tracerGroupSurfaceFlux,  &
!                       config_cvmix_kpp_nonlocal_with_implicit_mix, err)

               if (trim(groupItr % memberName) == 'activeTracers') then
                  !call ocn_tracer_vmix_tend_implicit(meshPool, dt, vertDiffTopOfCell, layerThickness, tracersGroup, &
                  !         vertNonLocalFlux, nonLocalSurfaceTracerFlux,  &
                  !         config_cvmix_kpp_nonlocal_with_implicit_mix, err)
                  call ocn_tracer_vmix_tend_implicit(dt, layerThickness, tracersGroup, &
                           nonLocalSurfaceTracerFlux,  &
                           err)
               else
                  !modifiedGroupName = trim(groupItr % memberName) // "SurfaceFlux"
                  !call ocn_tracer_vmix_tend_implicit(meshPool, dt, vertDiffTopOfCell, layerThickness, tracersGroup, &
                  !         vertNonLocalFlux, activeTracersSurfaceFlux,  &
                  !         config_cvmix_kpp_nonlocal_with_implicit_mix, err)
                  call ocn_tracer_vmix_tend_implicit(dt, layerThickness, tracersGroup, &
                           activeTracersSurfaceFlux,  &
                           err)
               endif
            end if

            ! difference tracers to compute influence of vertical mixing and divide by dt
            if (trim(groupItr % memberName) == 'activeTracers') then
               if (config_compute_active_tracer_budgets) then
#ifdef MPAS_OPENACC
                  !$acc parallel loop &
                  !$acc    present(activeTracerVertMixTendency,tracersGroup) &
                  !$acc    collapse(2)
#else
                  !$omp parallel
                  !$omp do schedule(runtime)
#endif
                  do iCell = 1, nCells
                     do k = 1,nVertLevels
                        do iTracer = 1,num_tracers
                           activeTracerVertMixTendency(iTracer,k,iCell) = &
                              (tracersGroup(iTracer,k,iCell) - activeTracerVertMixTendency(iTracer,k,iCell)) / dt
                        end do
                     end do
                  end do
#ifndef MPAS_OPENACC
                  !$omp end do
                  !$omp end parallel
#endif
               endif
            endif

         end if
      end do
      call mpas_timer_stop('vmix solve tracers')

      call mpas_timer_stop('vmix imp')

   end subroutine ocn_vmix_implicit!}}}

!***********************************************************************
!
!  routine ocn_vmix_init
!
!> \brief   Initializes ocean vertical mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_init(domain, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      integer :: err_tmp

      err = 0

      rayleighDampingCoef = 0.0_RKIND

      if (config_Rayleigh_friction) then
          rayleighDampingCoef = config_Rayleigh_damping_coeff
      endif

      rayleighBottomDampingCoef = 0.0_RKIND

      if (config_Rayleigh_bottom_friction) then
          rayleighBottomDampingCoef = config_Rayleigh_bottom_damping_coeff
      endif

      rayleighDepthVariable = 0.0_RKIND
      if (config_Rayleigh_damping_depth_variable) then
        rayleighDepthVariable = 1.0_RKIND
      end if

      velVmixOn = .true.
      tracerVmixOn = .true.

      if(config_disable_vel_vmix.or.config_disable_vel_all_tend) velVmixOn = .false.
      if(config_disable_tr_vmix.or.config_disable_tr_all_tend) tracerVmixOn = .false.

      implicitBottomDragCoef = 0.0_RKIND

      if (config_use_implicit_bottom_drag) then
          implicitBottomDragCoef = config_implicit_bottom_drag_coeff
      endif

      !$acc enter data copyin(implicitBottomDragCoef)

      call ocn_vmix_cvmix_init(domain,err_tmp)
      err = ior(err, err_tmp)
      call ocn_vmix_coefs_redi_init(err_tmp)
      err = ior(err, err_tmp)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_init!}}}

!***********************************************************************
!
!  routine tridiagonal_solve
!
!> \brief   Solve the matrix equation Ax=r for x, where A is tridiagonal.
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  Solve the matrix equation Ax=r for x, where A is tridiagonal.
!>  A is an nxn matrix, with:
!>  a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!>  b diagonal, filled from 1:n
!>  c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
!-----------------------------------------------------------------------
   subroutine tridiagonal_solve(a,b,c,r,x,n) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer,intent(in) :: n
      real (KIND=RKIND), dimension(n), intent(in) :: a,b,c,r

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (KIND=RKIND), dimension(n), intent(out) :: x

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (KIND=RKIND), dimension(n) :: bTemp,rTemp
      real (KIND=RKIND) :: m
      integer i

      ! Use work variables for b and r
      bTemp(1) = b(1)
      rTemp(1) = r(1)

      ! First pass: set the coefficients
      do i = 2,n
         m = a(i-1)/bTemp(i-1)
         bTemp(i) = b(i) - m*c(i-1)
         rTemp(i) = r(i) - m*rTemp(i-1)
      end do

      x(n) = rTemp(n)/bTemp(n)
       ! Second pass: back-substition
      do i = n-1, 1, -1
         x(i) = (rTemp(i) - c(i)*x(i+1))/bTemp(i)
      end do

   end subroutine tridiagonal_solve !}}}

!***********************************************************************
!
!  routine tridiagonal_solve_mult
!
!> \brief   Solve multiple matrix equations Ax=r for x, where A is tridiagonal.
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  Solve the matrix equation Ax=r for x, where A is tridiagonal.
!>  A is an nxn matrix, with:
!>  a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!>  b diagonal, filled from 1:n
!>  c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
!-----------------------------------------------------------------------
subroutine tridiagonal_solve_mult(a,b,c,r,x,n,nSystems)!{{{
   !$acc routine

   integer,intent(in) :: n, nSystems
   real (KIND=RKIND), dimension(n), intent(in) :: a,b,c
   real (KIND=RKIND), dimension(nSystems,n), intent(in) :: r
   real (KIND=RKIND), dimension(nSystems,n), intent(out) :: x
   real (KIND=RKIND), dimension(n) :: bTemp
   real (KIND=RKIND), dimension(nSystems,n) :: rTemp
   real (KIND=RKIND) :: m
   integer i,j

   ! Use work variables for b and r
   bTemp(1) = b(1)

   !$acc loop
   do j = 1,nSystems
      rTemp(j,1) = r(j,1)
   end do

   ! First pass: set the coefficients
   !$acc loop
   do i = 2,n
      m = a(i-1)/bTemp(i-1)
      bTemp(i) = b(i) - m*c(i-1)
      do j = 1,nSystems
         rTemp(j,i) = r(j,i) - m*rTemp(j,i-1)
      end do
   end do

   !$acc loop
   do j = 1,nSystems
      x(j,n) = rTemp(j,n)/bTemp(n)
   end do
   ! Second pass: back-substition
   !$acc loop
   do i = n-1, 1, -1
      do j = 1,nSystems
         x(j,i) = (rTemp(j,i) - c(i)*x(j,i+1))/bTemp(i)
      end do
   end do

end subroutine tridiagonal_solve_mult!}}}

!***********************************************************************
!
!  subroutine ocn_compute_kpp_rhs
!
!> \brief   Computes the non local flux tendency for KPP
!> \author  Luke Van Roekel
!> \date    October 2017
!> \details
!>   Computes non local flux tendency from KPP when
!>   config_cvmix_kpp_nonlocal_with_implicit_mix = .true.
!>   otherwise this term is computed in ocn_tend_tracer
!
!-----------------------------------------------------------------------

subroutine ocn_compute_kpp_rhs(tracers, rhs, dt, maxLevelCell, nTracers, &
                         layerThickness, vertNonLocalFlux, tracerGroupSurfaceFlux)!{{{

  real (kind=RKIND), intent(in) :: dt
  real (kind=RKIND), dimension(:,:), intent(in) :: vertNonLocalFlux, tracers
  real (kind=RKIND), dimension(:), intent(in) :: layerThickness, tracerGroupSurfaceFlux
  real (kind=RKIND), dimension(:,:), intent(out) :: rhs
  integer, intent(in) :: maxLevelCell, nTracers
  integer :: iTracer, k

  do k=2,maxLevelCell-1
     do iTracer=1,nTracers
        rhs(iTracer, k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) *   &
                     (vertNonLocalFlux(1,k) - vertNonLocalFlux(1,k+1)) / layerThickness(k)
     enddo
  enddo

  k=1
  do iTracer=1,nTracers
     rhs(iTracer, k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) *  &
                              (-vertNonLocalFlux(1,k+1) )/ layerThickness(k)
  enddo

  k=maxLevelCell
  do iTracer=1,nTracers
     rhs(iTracer,k) = tracers(iTracer,k) + dt * tracerGroupSurfaceFlux(iTracer) * &
                                     vertNonLocalFlux(1,k) / layerThickness(k)
  enddo

end subroutine ocn_compute_kpp_rhs!}}}
!***********************************************************************

end module ocn_vmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
