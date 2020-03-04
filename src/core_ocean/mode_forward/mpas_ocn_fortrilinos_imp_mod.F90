module ocn_fortrilinos_imp_mod
!---------------------------
#include "ForTrilinosInterface_config.hpp"
#include "ForTrilinos.h"
  !use ISO_FORTRAN_ENV
  use,intrinsic :: iso_c_binding
  use fortrilinos
  use forteuchos
  use fortpetra
  use fortest
  use ocn_belos_mod
  use ocn_model_evaluator_mod
!---------------------------
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar
!---------------------------

  implicit none
  private  
  save

  public :: ocn_time_integration_imp_btrmode

!*********************************************************************
  contains
!*********************************************************************

  subroutine ocn_time_integration_imp_btrmode(domain,dt,xstate,nvec,n_tot_vec)
   
  implicit none
  integer :: my_rank, num_procs

  integer(global_size_type) :: n_global
  integer(size_type) :: max_entries_per_row, num_vecs = 1, lda
  integer(int_type) :: row_nnz

  integer :: n,numvalid
  integer(global_ordinal_type) :: offset,icol,irow,icol1,icol2,irow1,irow2,gblrow

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist, linear_solver_list, belos_list, solver_list, krylov_list
  type(TrilinosSolver) :: solver_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A
  type(TpetraMultiVector) :: B, X, residual
! class(ForTpetraOperator), allocatable :: op

  real(scalar_type), dimension(:), allocatable :: lhs, rhs
  real(norm_type), dimension(:), allocatable :: norms
  integer(global_ordinal_type), dimension(1) :: cols
  real(scalar_type), dimension(1) :: vals
  real(scalar_type) :: r0, sone = 1., szero = 0., tol, val

  ! For MPAS-O -----------------------------------------------------------------
  type (domain_type) :: domain
  type (dm_info) :: dminfo
  real (kind=RKIND)  :: dt
  type (block_type), pointer :: block

  type (mpas_pool_type), pointer :: statePool
  type (mpas_pool_type), pointer :: meshPool
  type (mpas_pool_type), pointer :: diagnosticsPool
  type (mpas_pool_type), pointer :: tendPool
! type (mpas_pool_type), pointer :: tracersPool
! type (mpas_pool_type), pointer :: tracersTendPool

  integer :: nCells, nEdges
  integer :: i,j,iCell,iEdge,cell1,cell2
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
  character (len=*), parameter :: iterGroupName = 'iterFields'


  ! For User defined variables -------------------------------------------------
  real (c_double), dimension(nvec) :: xstate
  integer (c_int) :: ierr,nvec,n_tot_vec
  integer, dimension(:), allocatable :: globalIdx
  integer :: sCellIdx,eCellIdx,mpi_ierr
  logical :: init_belos = .true.
  ! ----------------------------------------------------------------------------
   
  ! INIT belos -----------------------------------------------------------------
  if ( init_belos ) then
    print*, 'PRINT in init'
    init_belos = .false.
    dminfo = domain % dminfo
    comm = TeuchosComm(dminfo % comm)
    my_rank = comm%getRank()
    num_procs = comm%getSize()

    ! Defining global index ----------------------------------------------------
      block => domain % blocklist
      do while (associated(block))
        call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
        nCells = nCellsArray(1)
        n = nCells
        block => block % next
      end do

      call MPI_BARRIER(dminfo%comm,mpi_ierr)

      do i = 0,num_procs - 1
        if ( my_rank == i ) then

          if ( my_rank == 0 ) then
            call MPI_SEND(nCells,1,MPI_INTEGER,my_rank+1,1,dminfo%comm,mpi_ierr)
            eCellIdx = nCells
          elseif ( my_rank == num_procs-1 ) then
            call MPI_RECV(eCellIdx,1,MPI_INTEGER,my_rank-1,1,dminfo%comm,MPI_STATUS_IGNORE,mpi_ierr)
            eCellIdx = eCellIdx + nCells
          else
            call MPI_RECV(eCellIdx,1,MPI_INTEGER,my_rank-1,1,dminfo%comm,MPI_STATUS_IGNORE,mpi_ierr)
            eCellIdx = eCellIdx + nCells 
            call MPI_SEND(eCellIdx,1,MPI_INTEGER,my_rank+1,1,dminfo%comm,mpi_ierr)
          endif
        
        endif 
        call MPI_BARRIER(dminfo%comm,mpi_ierr)
      end do

      do i = 0,num_procs - 1
        if ( my_rank == i ) then
          sCellIdx = eCellIdx - nCells + 1
        endif 
        call MPI_BARRIER(dminfo%comm,mpi_ierr)
      end do

      ! My global cell ID ------------------------------------------------------
      allocate(globalIdx(nCellsArray(4)))
      do iCell = 1,nCells
        globalIdx(iCell) = sCellIdx + iCell - 1
        CGvec_r0(iCell) = real(globalIdx(iCell))
      end do

      call mpas_dmpar_exch_group_create(domain, iterGroupName)
      call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_r0', 1)
      call mpas_threading_barrier()
      call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
      call mpas_dmpar_exch_group_destroy(domain, iterGroupName)

      do iCell = 1,nCellsArray(4)
        globalIdx(iCell) = int(CGvec_r0(iCell)) 
      end do

      ! ------------------------------------------------------------------------

  endif ! INIT_belos

  ! ----------------------------------------------------------------------------


  ! Read in the parameterList
  plist = ParameterList("Stratimikos"); FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist, "stratimikos.xml"); FORTRILINOS_CHECK_IERR()


  ! Get tolerance from the parameter list
  linear_solver_list = plist%sublist('Linear Solver Types')
  belos_list = linear_solver_list%sublist(plist%get_string('Linear Solver Type'))
  solver_list = belos_list%sublist('Solver Types')
  krylov_list = solver_list%sublist(belos_list%get_string('Solver Type'))

  tol = krylov_list%get_real("Convergence Tolerance")


  ! ------------------------------------------------------------------
  ! Step 0: Construct tri-diagonal matrix
  n_global = -1
  map = TpetraMap(n_global, n, comm); FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 8
  A = TpetraCrsMatrix(map, max_entries_per_row, TpetraStaticProfile)

! allocate(cols(max_entries_per_row))
! allocate(vals(max_entries_per_row))


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

     call A%setAllToScalar(szero)

     do iCell = 1, nCells

        gblrow  = globalIdx(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Interpolation sshEdge
          sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))

          ! method 1, matches method 0 without pbcs, works with pbcs.
          thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))

          !--------------------------------------------------------------!
          fluxAx = edgeSignOnCell(i,iCell)*dvEdge(iEdge)*thicknessSumLag / dcEdge(iEdge)

!         CGvec_r0(cell1) = CGvec_r0(cell1) + fluxAx
!         CGvec_r0(cell2) = CGvec_r0(cell2) - fluxAx
          !--------------------------------------------------------------!
     
            cols(1) = globalIdx(cell1)
            vals(1) = fluxAx 
            numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)

            cols(1) = globalIdx(cell2)
            vals(1) = -fluxAx 
            numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)
        end do ! i

!       CGvec_r0(iCell) = CGvec_r0(iCell) - (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
!       call A%insertGlobalValues(offset + i, cols(1:row_nnz-1), vals(1:row_nnz-1))  

        cols(1) = globalIdx(iCell)
        vals(1) = - (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
        numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)

     end do ! iCell

     ! End reduction ----------------------------------------------------------!

     block => block % next
  end do  ! block

  ! ---------------------------------------------------------------------------!

! do i = 1, n
!   row_nnz = 1
!   if (i .ne. 1 .or. my_rank > 0) then
!     cols(row_nnz) = offset + i-1
!     vals(row_nnz) = -1.0
!     row_nnz = row_nnz + 1
!   end if
!   cols(row_nnz) = offset + i
!   vals(row_nnz) = 2.0
!   row_nnz = row_nnz + 1
!   if (i .ne. n .or. my_rank .ne. num_procs-1) then
!     cols(row_nnz) = offset + i+1
!     vals(row_nnz) = -1.0
!     row_nnz = row_nnz + 1
!   end if

!   call A%insertGlobalValues(offset + i, cols(1:row_nnz-1), vals(1:row_nnz-1)); FORTRILINOS_CHECK_IERR()
! end do
  call A%fillComplete(); FORTRILINOS_CHECK_IERR()


  ! The solution X(i) = i-1
  allocate(lhs(nCells))
  allocate(rhs(nCells))


! if (my_rank > 0) then
!   rhs(1) = 0.0
! else
!   rhs(1) = -1.0
! end if
! if (my_rank .ne. num_procs-1) then
!   rhs(n) = 0.0
! else
!   rhs(n) = offset+n
! end if
! do i = 2, n-1
!   rhs(i) = 0.0
! end do

  lda = nCells
  B = TpetraMultiVector(map, rhs, lda, num_vecs); FORTRILINOS_CHECK_IERR()
  X = TpetraMultiVector(map, num_vecs); FORTRILINOS_CHECK_IERR()
  residual = TpetraMultiVector(map, num_vecs, .false.); FORTRILINOS_CHECK_IERR()

  allocate(norms(1))


  ! Step 0: create a handle
  solver_handle = TrilinosSolver(); FORTRILINOS_CHECK_IERR()


  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------

  ! Step 1: initialize a handle
  call solver_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  call solver_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver
  call solver_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  ! Step 4: solve the system
! call X%randomize()
  call X%putScalar(szero)

  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()

  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  r0 = norms(1)

  print*, r0
  call MPI_BARRIER(dminfo%comm,mpi_ierr)
  stop

  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()


  ! Check the solution
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()

  if (norms(1)/r0 > tol) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  end if

  ! Step 5: clean up
  call solver_handle%finalize(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------

  call solver_handle%release(); FORTRILINOS_CHECK_IERR()
  call plist%release(); FORTRILINOS_CHECK_IERR()
  call X%release(); FORTRILINOS_CHECK_IERR()
  call B%release(); FORTRILINOS_CHECK_IERR()
  call A%release(); FORTRILINOS_CHECK_IERR()
  call map%release(); FORTRILINOS_CHECK_IERR()
  call comm%release(); FORTRILINOS_CHECK_IERR()
  deallocate(norms)
! deallocate(cols)
! deallocate(vals)
  deallocate(lhs)
  deallocate(rhs)

  stop
! n_global_vec = nvec
! n_total_vec  = n_tot_vec
! !--- Only for initial
! if ( init_belos ) then
! print*, 'PRINT in init'
! init_nox = .false.

  end subroutine ocn_time_integration_imp_btrmode

end module










!  public :: noxsolve,noxfinish
!  class(ForModelEvaluator),private,allocatable :: model_evaluator
!
!!*********************************************************************
!                                contains
!!*********************************************************************
!
!  subroutine ocn_time_integration_imp_btrmode(domain,dt,xstate,nvec,n_tot_vec)
!   
!    implicit none
!    !--------------
!    type (TeuchosComm) :: comm
!    type (ParameterList) :: params
!    !--------------
!    type (domain_type) :: domain
!    type (dm_info) :: dminfo
!    real (kind=RKIND)  :: dt
!    real (c_double), dimension(nvec) :: xstate
!    integer (c_int) :: ierr,nvec,n_tot_vec
!    integer(global_size_type) :: n_global_vec,n_global2_vec,n_total_vec
!    integer(global_size_type) :: nCells1_vec,nCells4_vec
!    logical :: init_nox = .true.
!    ! -------------------------------------------------------------------
!   
!    n_global_vec = nvec
!     n_total_vec = n_tot_vec
!
!    !--- Only for initial
!    if ( init_nox ) then
!   
!!     print*, 'PRINT in init'
!
!      init_nox = .false.
!
!      dminfo = domain % dminfo
!
!      comm = TeuchosComm(dminfo % comm)
!      params = ParameterList("ocnModelEvaluator")
!      call load_from_xml(params,'nox_params.xml')
!
!      allocate(model_evaluator, source=ocnModelEvaluator(comm,domain,xstate,dt,n_global_vec,n_total_vec))
!      call init_ForModelEvaluator(model_evaluator)
!      call model_evaluator%setup(params)
!      call params%release()
!    endif ! init_nox
!    
!
!    call noxsolve(comm,domain,dt,params,xstate,nvec,ierr)
!
!  end subroutine ocn_time_integration_imp_btrmode
!
!!------------------------------------------------------------------------
!
!  subroutine noxsolve(comm,domain,dt,params,xstate,nvec,ierr)
!    ! -------------------------------------------------------------------
!    use,intrinsic :: iso_c_binding
!    implicit none
!
!    type (TeuchosComm) :: comm
!    type (domain_type) :: domain
!    real (kind=RKIND)  :: dt
!    integer(c_int) :: vector_size,ierr,nvec
!    real(c_double),dimension(nvec) :: xstate
!
!    type(NOXSolver) :: nox_solver
!    type(ocnModelEvaluator) :: self
!    integer(kind(NOXStatusType)) :: status
!    type(ParameterList) :: params
!    integer(global_size_type) :: n_global_vec
!    ! -------------------------------------------------------------------
!
!    ierr = 0
!
!    n_global_vec = nvec
!
!!     print*, 'PRINT in noxsolve'
!    
!    ! ------------------------------------------------- !
!    call mpas_timer_start("si fortrilinos params")
!    params = ParameterList("ocnModelEvaluator")
!    call load_from_xml(params,'nox_params.xml')
!    call mpas_timer_stop ("si fortrilinos params")
!    ! ------------------------------------------------- !
!
!    ! ------------------------------------------------- !
!    call mpas_timer_start("si fortrilinos evalse")
!    call model_evaluator%setup(params)
!    call mpas_timer_stop ("si fortrilinos evalse")
!    ! ------------------------------------------------- !
!
!    nox_solver = NOXSolver(model_evaluator)
!
!    ! ------------------------------------------------- !
!    call mpas_timer_start("si fortrilinos nox setup")
!    call nox_solver%setup(params)
!    call mpas_timer_stop ("si fortrilinos nox setup")
!    ! ------------------------------------------------- !
!
!    ! ------------------------------------------------- !
!    call mpas_timer_start("si fortrilinos nox solver")
!    status = nox_solver%solve()
!    call mpas_timer_stop ("si fortrilinos nox solver")
!    ! ------------------------------------------------- !
!
!    ! ------------------------------------------------- !
!    call mpas_timer_start("si fortrilinos nox release")
!    call nox_solver%release()
!    call mpas_timer_stop ("si fortrilinos nox release")
!    ! ------------------------------------------------- !
!
!    if (status /= NOXConverged) ierr = 1
!!   print*, "Status@@",status,NOXConverged
!
!    !--- Get the final solution
!    call getFinalSolution(self,xstate,n_global_vec)
!
!    call params%release()
!
!  end subroutine noxsolve
!
!!------------------------------------------------------------------------
!
!  subroutine noxfinish()
!    call model_evaluator%release()
!    deallocate(model_evaluator)
!  end subroutine noxfinish
!
!end module
