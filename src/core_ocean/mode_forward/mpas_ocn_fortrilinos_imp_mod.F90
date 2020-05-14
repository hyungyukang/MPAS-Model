module ocn_fortrilinos_imp_mod
!---------------------------
#include "ForTrilinosInterface_config.hpp"
#include "ForTrilinos.h"
  use,intrinsic :: iso_c_binding
  use fortrilinos
  use forteuchos
  use fortpetra
  use fortest
!---------------------------
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar
  use mpas_timer
!---------------------------

  implicit none
  private  
  save

  public :: ocn_time_integration_imp_btrmode

!*********************************************************************
  contains
!*********************************************************************

  subroutine ocn_time_integration_imp_btrmode(domain,dt,n_tot_vec,stage)
   
  implicit none
  integer :: my_rank, num_procs

  integer(global_size_type) :: n_global,n_local
  integer(size_type),save :: max_entries_per_row, num_vecs = 1, lda,ione = 1,itwo=2,izero=0
  integer(int_type) :: row_nnz

  integer :: n,numvalid
  integer(global_ordinal_type) :: offset,icol,irow,icol1,icol2,irow1,irow2,gblrow

  type(TeuchosComm),save :: comm
  type(ParameterList),save:: plist_o, linear_solver_list_o, belos_list_o, solver_list_o, krylov_list_o
  type(ParameterList),save:: plist_m, linear_solver_list_m, belos_list_m, solver_list_m, krylov_list_m
  type(TrilinosSolver),save :: solver_handle_o,solver_handle_m
  type(TpetraMap),save :: map
  type(TpetraCrsMatrix),save :: A
  type(TpetraMultiVector),save :: B, X, residual

  real(scalar_type), dimension(:), allocatable :: lhs, rhs
  real(scalar_type), dimension(:), pointer :: solvec
  real(norm_type), dimension(:) :: norms(1)
  integer(global_ordinal_type) :: cols(1)
  real(scalar_type) :: vals(1)
  real(scalar_type) :: r0, sone = 1., szero = 0., tol_o,tol_m, val

  ! For MPAS-O -----------------------------------------------------------------
  type (domain_type) :: domain
  type (dm_info) :: dminfo
  real (kind=RKIND)  :: dt
  type (block_type), pointer :: block

  type (mpas_pool_type), pointer :: statePool
  type (mpas_pool_type), pointer :: meshPool
  type (mpas_pool_type), pointer :: diagnosticsPool
  type (mpas_pool_type), pointer :: tendPool

  integer :: nCells, nEdges
  integer :: i,j,iCell,iEdge,cell1,cell2,isum
  integer :: sshEdgeLag1,sshEdgeLag2
  integer :: thicknessSumLag1,thicknessSumLag2
  integer :: sshDiffNew1, sshDiffNew2,sshDiffNew

  integer,pointer :: nCellsPtr, nEdgesPtr
  integer, dimension(:),pointer :: nCellsArray
  integer, dimension(:),pointer :: nEdgesArray
  integer, dimension(:),pointer :: nEdgesOnCell
  integer, dimension(:,:),pointer :: cellsOnEdge
  integer, dimension(:,:),pointer :: edgesOnCell
  integer, dimension(:,:),pointer :: edgeSignOnCell
  real (kind=RKIND) :: sshTendb1,sshTendb2,sshTendAx,sshEdgeCur,thicknessSumCur,sshDiffCur
  real (kind=RKIND) :: sshEdgeMid,sshEdgeLag,thicknessSumMid,thicknessSumLag,sshDiffLag
  real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshCurArea,sshLagArea
  real (kind=RKIND), dimension(:),pointer :: dcEdge
  real (kind=RKIND), dimension(:),pointer :: bottomDepth
  real (kind=RKIND), dimension(:),pointer :: dvEdge
  real (kind=RKIND), dimension(:),pointer :: areaCell

  real (kind=RKIND), dimension(:),pointer :: sshCur,sshNew
  real (kind=RKIND), dimension(:),pointer :: sshSubcycleCur,sshSubcycleNew
  real (kind=RKIND), dimension(:),pointer :: normalBarotropicVelocityCur
  real (kind=RKIND), dimension(:),pointer :: barotropicForcing
  real (kind=RKIND), dimension(:),pointer :: barotropicCoriolisTerm
  real (kind=RKIND), dimension(:),pointer :: CGvec_r0,CGvec_r1
  integer(global_ordinal_type), dimension(:),allocatable :: acol,colent
  real (scalar_type), dimension(:),allocatable :: aval,valent
  character (len=*), parameter :: iterGroupName = 'iterFields'
  character (len=1) :: stage 


  ! For User defined variables -------------------------------------------------
  integer :: ierr,nvec,n_tot_vec
  integer, dimension(:), allocatable,save :: globalIdx
  integer :: sCellIdx,eCellIdx,mpi_ierr
  logical :: init_belos = .true., init_belos_o = .true., init_belos_m = .true.
  ! ----------------------------------------------------------------------------

  dminfo = domain % dminfo
   
  ! INIT belos :: Initial only ===================================================================!
  if ( init_belos ) then
    print*, 'PRINT in init'
    init_belos = .false.
    comm = TeuchosComm(dminfo % comm)
    my_rank = comm%getRank()
    num_procs = comm%getSize()

    ! Defining global index ----------------------------------------------------
    block => domain % blocklist
    do while (associated(block))
      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )

      if (num_procs == 1) then

        allocate(globalIdx(nCellsArray(4)+1))
        do iCell = 1,nCellsArray(4)
          globalIdx(iCell) = iCell
        end do
          globalIdx(nCellsArray(4)+1) = 0

      else

        call MPI_BARRIER(dminfo%comm,mpi_ierr)
        do i = 0,num_procs - 1
          if ( my_rank == i ) then
  
            if ( my_rank == 0 ) then
              call MPI_SEND(nCellsArray(1),1,MPI_INTEGER,my_rank+1,1,dminfo%comm,mpi_ierr)
              eCellIdx = nCellsArray(1)
            elseif ( my_rank == num_procs-1 ) then
              call MPI_RECV(eCellIdx,1,MPI_INTEGER,my_rank-1,1,dminfo%comm,MPI_STATUS_IGNORE,mpi_ierr)
              eCellIdx = eCellIdx + nCellsArray(1)
            else
              call MPI_RECV(eCellIdx,1,MPI_INTEGER,my_rank-1,1,dminfo%comm,MPI_STATUS_IGNORE,mpi_ierr)
              eCellIdx = eCellIdx + nCellsArray(1)
              call MPI_SEND(eCellIdx,1,MPI_INTEGER,my_rank+1,1,dminfo%comm,mpi_ierr)
            endif ! my_rank
          
        endif  ! my_rank
        call MPI_BARRIER(dminfo%comm,mpi_ierr)
        end do ! i
  
  
        do i = 0,num_procs - 1
          if ( my_rank == i ) then
            sCellIdx = eCellIdx - nCellsArray(1) + 1
          endif 
          call MPI_BARRIER(dminfo%comm,mpi_ierr)
        end do
  
  
        ! My global cell ID ------------------------------------------------------
        allocate(globalIdx(nCellsArray(4)+1))
        globalIdx(:) = 0
        CGvec_r0(:) = 0.0
        do iCell = 1,nCellsArray(1)
          globalIdx(iCell) = sCellIdx + iCell - 1
          CGvec_r0(iCell) = real(globalIdx(iCell))
        end do
  
        call mpas_dmpar_exch_group_create(domain, iterGroupName)
        call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_r0', 1)
        call mpas_threading_barrier()
        call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
        call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
  
        do iCell = 1,nCellsArray(4)+1
          globalIdx(iCell) = int(CGvec_r0(iCell)) 
        end do

        CGvec_r0(:) = 0.0_RKIND

      endif ! num_procs

      block => block % next
    end do

      ! ------------------------------------------------------------------------

      lda = nCellsArray(1)


  ! ----------------------------------------------------------------------------
  ! Step 0: Construct coefficient matrix
  n_global = -1
  map = TpetraMap(n_global, nCellsArray(1), comm) !; FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 8
  A = TpetraCrsMatrix(map, max_entries_per_row, TpetraStaticProfile)

  ! -- MPAS-O SpMV -------------------------------------------------------------
  block => domain % blocklist
  do while (associated(block))

     call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

     call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
     call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

     call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
     call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
     call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
     call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
     call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
     call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
     call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
     call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )

     call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
     call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
     call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)

     call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
     call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1' , CGvec_r1 )

     nCells = nCellsArray(1)
     nEdges = nEdgesArray(2)

     vals(1) = szero

     do iCell = 1, nCellsArray(1)
        gblrow  = globalIdx(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
          thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1),bottomDepth(cell2))
          fluxAx = edgeSignOnCell(i, iCell) * (thicknessSumLag / dcEdge(iEdge)) * dvEdge(iEdge)

          if ( globalIdx(cell2) >  0 ) then
          cols(1) = globalIdx(cell1)
          call A%insertGlobalValues(gblrow, cols, vals)

          cols(1) = globalIdx(cell2)
          call A%insertGlobalValues(gblrow, cols, vals)
          endif

        end do ! i

        cols(1) = globalIdx(iCell)
        call A%insertGlobalValues(gblrow, cols, vals)
          
     end do ! iCell

     call A%fillComplete() !; FORTRILINOS_CHECK_IERR()

     block => block % next
  end do  ! block

  ! Read in the parameterList
  plist_o = ParameterList("Stratimikos") !; FORTRILINOS_CHECK_IERR()
  plist_m = ParameterList("Stratimikos") !; FORTRILINOS_CHECK_IERR()

  call load_from_xml(plist_o, "stratimikos.xml") !; FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist_m, "stratimikos.xml") !; FORTRILINOS_CHECK_IERR()

  ! Get tolerance from the parameter list
  linear_solver_list_o = plist_o%sublist('Linear Solver Types')
  linear_solver_list_m = plist_m%sublist('Linear Solver Types')

  belos_list_o = linear_solver_list_o%sublist(plist_o%get_string('Linear Solver Type'))
  belos_list_m = linear_solver_list_m%sublist(plist_m%get_string('Linear Solver Type'))

  solver_list_o = belos_list_o%sublist('Solver Types')
  solver_list_m = belos_list_m%sublist('Solver Types')

  krylov_list_o = solver_list_o%sublist(belos_list_o%get_string('Solver Type'))
  krylov_list_m = solver_list_m%sublist(belos_list_m%get_string('Solver Type'))

  call krylov_list_o%set('Convergence Tolerance', 1e-2)
  tol_o = 1.0d-2
  call krylov_list_m%set('Convergence Tolerance', 1e-8)
  tol_m = 1.0d-8


  ! Trilinos solver handle:  'o' for outer iteration, 'm' for main iteration
  solver_handle_o = TrilinosSolver() !; FORTRILINOS_CHECK_IERR()
  solver_handle_m = TrilinosSolver() !; FORTRILINOS_CHECK_IERR()

  call solver_handle_o%init(comm) !; FORTRILINOS_CHECK_IERR()
  call solver_handle_m%init(comm) !; FORTRILINOS_CHECK_IERR()
  endif ! INIT_belos :: Initial only =============================================================!


  call mpas_timer_start("fort mat setup")

  ! -- MPAS-O SpMV -------------------------------------------------------------
  block => domain % blocklist
  do while (associated(block))

     call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

     call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
     call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

     call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
     call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
     call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
     call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
     call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
     call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
     call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
     call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )

     call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
     call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
     call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)

     call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
     call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1' , CGvec_r1 )

     nCells = nCellsArray(1)
     nEdges = nEdgesArray(2)
 

     ! A : Coefficient matrix -------------------------------------------------!

     call A%resumeFill() !; FORTRILINOS_CHECK_IERR()

     call mpas_timer_start("fort mat AllToScalar")
     call A%setAllToScalar(szero)
     call mpas_timer_stop("fort mat AllToScalar")

     do iCell = 1, nCellsArray(1)

        gblrow  = globalIdx(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Interpolation sshEdge
          sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
          thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1),bottomDepth(cell2))
          fluxAx = edgeSignOnCell(i, iCell) * (thicknessSumLag / dcEdge(iEdge)) * dvEdge(iEdge)

          if ( globalIdx(cell2) > 0) then
          cols(1) = globalIdx(cell1)
          vals(1) = -fluxAx 
          numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)

          cols(1) = globalIdx(cell2)
          vals(1) = +fluxAx 
          numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)
          endif

        end do ! i

        cols(1) = globalIdx(iCell)
        vals(1) = (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
        numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)

     end do ! iCell

  call A%fillComplete() !; FORTRILINOS_CHECK_IERR()

  call mpas_timer_stop("fort mat setup")

  call mpas_timer_start("fort resid")

     ! B : Right Hand Side ----------------------------------------------------!
     do iCell = 1, nCellsArray(1)

       sshTendb1 = 0.0_RKIND
       sshTendb2 = 0.0_RKIND

       do i = 1, nEdgesOnCell(iCell)
         iEdge = edgesOnCell(i, iCell)

         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         ! Interpolation sshEdge
         sshEdgeCur = 0.5_RKIND * (        sshCur(cell1) + sshCur(cell2))
         sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
         sshEdgeMid = 0.5_RKIND * (           sshEdgeLag + sshEdgeCur           )

         ! method 1, matches method 0 without pbcs, works with pbcs.
         thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
         thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
         thicknessSumMid = sshEdgeMid + min(bottomDepth(cell1), bottomDepth(cell2))

         ! nabla (ssh^0)
         sshDiffCur = (        sshCur(cell2)-        sshCur(cell1)) / dcEdge(iEdge)
         sshDiffLag = (sshSubcycleCur(cell2)-sshSubcycleCur(cell1)) / dcEdge(iEdge)

         !--------------------------------------------------------------!
         fluxb1 = thicknessSumMid * normalBarotropicVelocityCur(iEdge)
         fluxb2 = thicknessSumLag * (0.5*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)) )

         sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
         sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
         !--------------------------------------------------------------!

       end do ! i

         sshTendb1  = (2.0_RKIND/(gravity*dt*0.5)) * sshTendb1
         sshTendb2  = (2.0_RKIND/(gravity       )) * sshTendb2
         sshCurArea = (1.0_RKIND/(gravity*dt**2.0*0.5**2.0)) * sshCur(iCell) * areaCell(iCell)

         CGvec_r0(iCell) = -(-sshCurArea - sshTendb1 + sshTendb2) 
     end do ! iCell

  call mpas_timer_stop("fort resid")

  call mpas_timer_start("fort Vector")
  B = TpetraMultiVector(map, CGvec_r0(1:nCellsArray(1)), lda, num_vecs) !; FORTRILINOS_CHECK_IERR()
  X = TpetraMultiVector(map, sshSubcycleCur(1:nCellsArray(1)), lda, num_vecs) !; FORTRILINOS_CHECK_IERR()

  residual = TpetraMultiVector(map,num_vecs,.false.) !; FORTRILINOS_CHECK_IERR()
  call mpas_timer_stop("fort Vector")


  call mpas_timer_start("fort solver setup")

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------

  ! Step 2: setup the problem
  if ( stage == 'o' ) then
    call solver_handle_o%setup_matrix(A) !; FORTRILINOS_CHECK_IERR()
  elseif ( stage == 'm' ) then
    call solver_handle_m%setup_matrix(A) !; FORTRILINOS_CHECK_IERR()
  endif


  ! Step 3: setup the solver - Initial only
  if ( init_belos_o .and. stage == 'o') then
    call solver_handle_o%setup_solver(plist_o) !; FORTRILINOS_CHECK_IERR()
    init_belos_o = .false.
  elseif ( init_belos_m .and. stage == 'm') then
    call solver_handle_m%setup_solver(plist_m) !; FORTRILINOS_CHECK_IERR()
    init_belos_m = .false.
  endif


  call mpas_timer_stop("fort solver setup")

  call mpas_timer_start("fort init resid")
  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero) !; FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone) !; FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms) !; FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
  call mpas_timer_stop("fort init resid")

  call mpas_timer_start("fort solve")


  if ( stage == 'o' ) then
    call solver_handle_o%solve(B, X) !; FORTRILINOS_CHECK_IERR()
  elseif ( stage == 'm' ) then
    call solver_handle_m%solve(B, X) !; FORTRILINOS_CHECK_IERR()
  endif


  call mpas_timer_stop("fort solve")

  call mpas_timer_start("fort check")
  ! Check the solution
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero) !; FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone) !; FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms) !; FORTRILINOS_CHECK_IERR()
  call mpas_timer_stop("fort check")


  if ( stage == 'o' .and. norms(1)/r0 > tol_o) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  elseif ( stage == 'm' .and. norms(1)/r0 > tol_m) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  end if


  ! Get solution
  solvec => X%getData(ione)
  sshSubcycleCur(1:nCellsArray(1)) = solvec(1:nCellsArray(1))
  nullify(solvec)


  call mpas_timer_start("si halo ssh")
  call mpas_dmpar_exch_group_create(domain, iterGroupName)
  call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'sshSubcycle', 1 )
  call mpas_threading_barrier()
  call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
  call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
  call mpas_timer_stop("si halo ssh")


     block => block % next
  end do  ! block

! call mpas_timer_start("fort final")
! ! Step 5: clean up
! call solver_handle%finalize() !; FORTRILINOS_CHECK_IERR()
! call mpas_timer_stop("fort final")

! ------------------------------------------------------------------

! call krylov_list%release()
! call solver_list%release!; FORTRILINOS_CHECK_IERR()
! call belos_list%release!; FORTRILINOS_CHECK_IERR()
! call linear_solver_list%release!; FORTRILINOS_CHECK_IERR()
! call solver_handle%release()!; FORTRILINOS_CHECK_IERR()
! call plist%release()!; FORTRILINOS_CHECK_IERR()

! call X%release()!; FORTRILINOS_CHECK_IERR()
! call B%release()!; FORTRILINOS_CHECK_IERR()
! call A%release()!; FORTRILINOS_CHECK_IERR()
! call map%release()!; FORTRILINOS_CHECK_IERR()
! call comm%release()!; FORTRILINOS_CHECK_IERR()

! deallocate(globalIdx)
! deallocate(aval)
! deallocate(acol)
! deallocate(colent)
! deallocate(valent)

! deallocate(norms)
! deallocate(cols)
! deallocate(vals)
! deallocate(lhs)
! deallocate(rhs)

  end subroutine ocn_time_integration_imp_btrmode

end module

