module ocn_fortrilinos_imp_mod
!---------------------------
#include "ForTrilinosInterface_config.hpp"
#include "ForTrilinos.h"
  use ISO_FORTRAN_ENV
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

  integer :: ierr
  integer :: i, n
  integer(global_ordinal_type) :: offset

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist, linear_solver_list, belos_list, solver_list, krylov_list
  type(TrilinosSolver) :: solver_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A
  type(TpetraMultiVector) :: B, X, residual
  class(ForTpetraOperator), allocatable :: op

  real(scalar_type), dimension(:), allocatable :: lhs, rhs
  real(norm_type), dimension(:), allocatable :: norms
  integer(global_ordinal_type), dimension(:), allocatable :: cols
  real(scalar_type), dimension(:), allocatable :: vals
  real(scalar_type) :: r0, sone = 1., szero = 0., tol

  !--------------
  type (domain_type) :: domain
  type (dm_info) :: dminfo
  real (kind=RKIND)  :: dt
  real (c_double), dimension(nvec) :: xstate
  integer (c_int) :: ierr,nvec,n_tot_vec
! integer(global_size_type) :: n_global_vec,n_global2_vec,n_total_vec
! integer(global_size_type) :: nCells1_vec,nCells4_vec
  logical :: init_belos = .true.
  ! -------------------------------------------------------------------
  n = 10
   
  if ( init_belos ) then
    print*, 'PRINT in init'
    init_belos = .false.
    dminfo = domain % dminfo
    comm = TeuchosComm(dminfo % comm)
    my_rank = comm%getRank()
    num_procs = comm%getSize()
  endif


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

  max_entries_per_row = 3
  A = TpetraCrsMatrix(map, max_entries_per_row, TpetraStaticProfile)

  allocate(cols(max_entries_per_row))
  allocate(vals(max_entries_per_row))
  offset = n * my_rank
  do i = 1, n
    row_nnz = 1
    if (i .ne. 1 .or. my_rank > 0) then
      cols(row_nnz) = offset + i-1
      vals(row_nnz) = -1.0
      row_nnz = row_nnz + 1
    end if
    cols(row_nnz) = offset + i
    vals(row_nnz) = 2.0
    row_nnz = row_nnz + 1
    if (i .ne. n .or. my_rank .ne. num_procs-1) then
      cols(row_nnz) = offset + i+1
      vals(row_nnz) = -1.0
      row_nnz = row_nnz + 1
    end if

    call A%insertGlobalValues(offset + i, cols(1:row_nnz-1), vals(1:row_nnz-1)); FORTRILINOS_CHECK_IERR()
  end do
  call A%fillComplete(); FORTRILINOS_CHECK_IERR()

  ! The solution X(i) = i-1
  allocate(lhs(n))
  allocate(rhs(n))
  if (my_rank > 0) then
    rhs(1) = 0.0
  else
    rhs(1) = -1.0
  end if
  if (my_rank .ne. num_procs-1) then
    rhs(n) = 0.0
  else
    rhs(n) = offset+n
  end if
  do i = 2, n-1
    rhs(i) = 0.0
  end do
  lda = n

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
  call X%randomize()
  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
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
  ! Implicit (inversion-of-control) setup [ no solve ]
  ! ------------------------------------------------------------------
  ! We cannot use most preconditioners without a matrix, so we remove any from
  ! the parameter list. We also adjust the number of iterations so that it is
  ! sufficient for convergence
  call plist%set('Preconditioner Type', 'None')
  call krylov_list%set('Maximum Iterations', 333)

  allocate(op, source=TriDiagOperator(map, A%getColMap()))
  call init_ForTpetraOperator(op); FORTRILINOS_CHECK_IERR()

  ! Step 1: initialize a handle
  call solver_handle%init(comm); FORTRILINOS_CHECK_IERR()
  ! Step 2: setup the problem
  ! Implicit (inversion-of-control) setup
  call solver_handle%setup_operator(op); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver

  call solver_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  call krylov_list%release()

  ! Step 4: solve the system
  call X%randomize()
  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()

  ! Check the solution
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  if (norms(1)/r0 > tol) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 666
  end if

  ! Step 5: clean up
  call solver_handle%finalize(); FORTRILINOS_CHECK_IERR()

  call krylov_list%release; FORTRILINOS_CHECK_IERR()
  call solver_list%release; FORTRILINOS_CHECK_IERR()
  call belos_list%release; FORTRILINOS_CHECK_IERR()
  call linear_solver_list%release; FORTRILINOS_CHECK_IERR()
  call op%release(); FORTRILINOS_CHECK_IERR()
  deallocate(op)
  ! ------------------------------------------------------------------

  call solver_handle%release(); FORTRILINOS_CHECK_IERR()
  call plist%release(); FORTRILINOS_CHECK_IERR()
  call X%release(); FORTRILINOS_CHECK_IERR()
  call B%release(); FORTRILINOS_CHECK_IERR()
  call A%release(); FORTRILINOS_CHECK_IERR()
  call map%release(); FORTRILINOS_CHECK_IERR()
  call comm%release(); FORTRILINOS_CHECK_IERR()
  deallocate(norms)
  deallocate(cols)
  deallocate(vals)
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
