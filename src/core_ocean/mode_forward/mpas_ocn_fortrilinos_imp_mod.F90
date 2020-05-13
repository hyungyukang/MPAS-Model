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

  subroutine ocn_time_integration_imp_btrmode(domain,dt,n_tot_vec,area_mean,stage)
   
  implicit none
  integer :: my_rank, num_procs

  integer(global_size_type),save :: n_global,n_local
  integer(size_type),save :: max_entries_per_row, num_vecs = 1, lda,ione = 1,itwo=2,izero=0
  integer,save :: max_entries_per_row2
  integer(int_type),save :: row_nnz

  integer,save :: n,numvalid,numnnz_a,numrow_a
  integer,save :: numnnz_c,numrow_c
  integer,dimension(:),allocatable :: indchk_a,indlcl_a
  integer,dimension(:),allocatable :: indchk_c,indlcl_c

  integer(size_type), dimension(:), allocatable,save :: row_ptrs_a, rp_res_a
  integer(int_type), dimension(:), allocatable,save :: colind_a, col_res_a, col_gbl_a
  real(scalar_type), dimension(:), allocatable,save :: values_a, val_res_a

  integer(global_ordinal_type), dimension(:), allocatable,save :: row_ptrs_c, rp_res_c
  !integer(size_type), dimension(:), allocatable,save :: row_ptrs_c, rp_res_c
  integer(int_type), dimension(:), allocatable,save :: colind_c, col_res_c, col_gbl_c
  real(scalar_type), dimension(:), allocatable,save :: values_c, val_res_c

  integer :: row,col(1)
  integer,dimension(:),allocatable :: cole
  integer(global_ordinal_type) :: offset,icol,irow,icol1,icol2,irow1,irow2,gblrow

  type(TeuchosComm),save :: comm
  type(ParameterList),save:: plist_o, linear_solver_list_o, belos_list_o, solver_list_o, krylov_list_o
  type(ParameterList),save:: plist_m, linear_solver_list_m, belos_list_m, solver_list_m, krylov_list_m
  type(TrilinosSolver),save :: solver_handle_o,solver_handle_m
  type(TpetraMap),save :: map,colmap_a,colmap_c,mape
  type(TpetraCrsMatrix),save :: A,C
  type(TpetraCrsGraph),save :: graph,graph1
  type(TpetraMultiVector),save :: B, X, residual

  real(scalar_type), dimension(:), allocatable :: lhs, rhs
  real(scalar_type), dimension(:), pointer :: solvec
  real(norm_type), dimension(:) :: norms(1)
  integer(global_ordinal_type) :: cols(1)
  !integer(global_ordinal_type),dimension(:),allocatable,save :: colse
  integer,dimension(:),allocatable,save :: colse
  real(scalar_type),dimension(:),allocatable,save :: valse
  real(scalar_type) :: vals(1)
  real(scalar_type) :: r0, sone = 1., szero = 0., tol_o,tol_m, val

  ! For MPAS-O -----------------------------------------------------------------
  type (domain_type) :: domain
  type (dm_info) :: dminfo
  real (kind=RKIND)  :: dt,area_mean
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
  integer(global_ordinal_type), dimension(:), allocatable,save :: globalIdx_fort
  integer :: sCellIdx,eCellIdx,mpi_ierr
  logical,save :: init_belos = .true., init_belos_o = .true., init_belos_m = .true.

  integer :: numrow, numnnz, id
  integer(global_ordinal_type), dimension(:), allocatable :: rowptr, rowchk
  integer, dimension(:), allocatable :: ind, indchk, indgbl, indlcl,rowptrs
  real(scalar_type), dimension(:), allocatable :: vala, valchk

  type(TpetraMap),save :: colmap
  integer(int_type), dimension(:), allocatable :: colind
  real(kind=RKIND) :: sum1,sum2
  character(kind=C_CHAR, len=:), allocatable :: description
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

      allocate(globalIdx_fort(nCellsArray(4)+1))
  
      do i = 1,nCellsArray(4)+1
        globalIdx_fort(i) = globalIdx(i)
      end do

      block => block % next
    end do

      ! ------------------------------------------------------------------------

      lda = nCellsArray(1)

  ! ----------------------------------------------------------------------------
  ! Step 0: Construct coefficient matrix
  n_global = -1
  !n_global = n_tot_vec
  map = TpetraMap(n_global, nCellsArray(1), comm)
  !mape = TpetraMap(n_global,nCellsArray(2), comm) !; FORTRILINOS_CHECK_IERR()
  !map = TpetraMap(n_global, globalIdx_fort(1:nCellsArray(1)), comm)
  mape = TpetraMap(n_global,globalIdx_fort(1:nCellsArray(2)), comm) !; FORTRILINOS_CHECK_IERR()
  max_entries_per_row = 8
  C = TpetraCrsMatrix(map,mape,max_entries_per_row, TpetraStaticProfile)
  graph = TpetraCrsGraph(map,mape,max_entries_per_row, TpetraStaticProfile)
  graph1 = TpetraCrsGraph(map,mape,max_entries_per_row, TpetraStaticProfile)
  !A = TpetraCrsMatrix(map,mape,max_entries_per_row, TpetraStaticProfile)
  !C = TpetraCrsMatrix(map,max_entries_per_row, TpetraStaticProfile)

  print*, 'AAA'

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

     !vals(1) = szero
     vals(1) = 1.d0

     do iCell = 1, nCellsArray(1)
        gblrow  = globalIdx(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)
          thicknessSumLag = min(bottomDepth(cell1),bottomDepth(cell2))
          fluxAx = edgeSignOnCell(i, iCell) * (thicknessSumLag / dcEdge(iEdge)) * dvEdge(iEdge)

          if ( globalIdx(cell2) >  0 ) then
          cols(1) = globalIdx(cell1)
          vals(1) = -fluxAx
          call C%insertGlobalValues(gblrow, cols, vals)
!         call graph%insertGlobalIndices(gblrow, cols)

          cols(1) = globalIdx(cell2)
          vals(1) = +fluxAx
          call C%insertGlobalValues(gblrow, cols, vals)
!         call graph%insertGlobalIndices(gblrow, cols)
          endif

        end do ! i

        cols(1) = globalIdx(iCell)
        vals(1) = (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
        call C%insertGlobalValues(gblrow, cols, vals)
!       call graph%insertGlobalIndices(gblrow, cols)
          
     end do ! iCell

     call C%fillComplete() !; FORTRILINOS_CHECK_IERR()
!    call graph%fillComplete() !; FORTRILINOS_CHECK_IERR()
     !call A%fillComplete() !; FORTRILINOS_CHECK_IERR()

! print*, 'BBB'
     graph = C%getCrsGraph()
     C = TpetraCrsMatrix(graph)

!    call C%resumeFill() !; FORTRILINOS_CHECK_IERR()

!    do iCell = 1, nCellsArray(1)

!       gblrow = globalIdx(iCell)
!       row = iCell

!       do i = 1, nEdgesOnCell(iCell)
!         iEdge = edgesOnCell(i, iCell)
!         cell1 = cellsOnEdge(1, iEdge)
!         cell2 = cellsOnEdge(2, iEdge)

!         ! Interpolation sshEdge
!         sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
!         thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1),bottomDepth(cell2))
!         fluxAx = edgeSignOnCell(i, iCell) * (thicknessSumLag / dcEdge(iEdge)) * dvEdge(iEdge)

!         if ( cell2 > 0 ) then
!           col(1) = cell1
!           vals(1) = 1.0
!           numvalid = C%sumIntoLocalValues(row,col,vals)
! 
!           col(1) = cell2
!           vals(1) = 1.0
!           numvalid = C%sumIntoLocalValues(row,col,vals)
!         endif

!       end do ! i

!       col(1) = iCell !map%getLocalElement((map%getGlobalElement(iCell)))
!       vals(1) = (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
!       numvalid = C%sumIntoLocalValues(row,col,vals)

!    end do ! iCell

!    call C%fillComplete() !; FORTRILINOS_CHECK_IERR()

     block => block % next
  end do  ! block


! Get all mat elements in CRS format using 'getAllValues' -------------------------------
  numrow_c = C%getNodeNumRows()
  numnnz_c = C%getNodeNumEntries()
 
  allocate(row_ptrs_c(numrow_c+1),rp_res_c(numrow_c+1) )
  allocate(  colind_c(numnnz_c)  ,col_res_c(numnnz_c),  col_gbl_c(numnnz_c) )
  allocate(  values_c(numnnz_c)  ,val_res_c(numnnz_c) )

  allocate(row_ptrs_a(numrow_c+1),rp_res_a(numrow_c+1) )
  allocate(  colind_a(numnnz_c)  ,col_res_a(numnnz_c),  col_gbl_a(numnnz_c) )
  allocate(  values_a(numnnz_c)  ,val_res_a(numnnz_c) )

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

! call mpas_timer_start("fort release")
! call A%release()!; FORTRILINOS_CHECK_IERR()
! call mpas_timer_stop("fort release")
!
! call mpas_timer_start("fort map")
! A = TpetraCrsMatrix(map,map,max_entries_per_row, TpetraStaticProfile)
! call mpas_timer_stop("fort map")

! call mpas_timer_start("fort setAllValues")
! call A%setAllValues(row_ptrs_c, colind_c, values_c)
! call mpas_timer_stop("fort setAllValues")


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

!    call A%resumeFill() !; FORTRILINOS_CHECK_IERR()

!    call mpas_timer_start("fort mat AllToScalar")
!    call A%setAllToScalar(szero)
!    call mpas_timer_stop("fort mat AllToScalar")


!    do iCell = 1, nCellsArray(1)

!       gblrow  = globalIdx(iCell)
!       row = iCell

!       do i = 1, nEdgesOnCell(iCell)
!         iEdge = edgesOnCell(i, iCell)
!         cell1 = cellsOnEdge(1, iEdge)
!         cell2 = cellsOnEdge(2, iEdge)

!         ! Interpolation sshEdge
!         sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
!         thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1),bottomDepth(cell2))
!         fluxAx = edgeSignOnCell(i, iCell) * (thicknessSumLag / dcEdge(iEdge)) * dvEdge(iEdge)
!

!         if ( globalIdx(cell2) > 0 ) then
!         cols(1) = globalIdx(cell1)
!         vals(1) = -fluxAx 
!         numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)

!         cols(1) = globalIdx(cell2)
!         vals(1) = +fluxAx 
!         numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)
!         endif

!       end do ! i

!       cols(1) = globalIdx(iCell)
!       vals(1) = (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
!       numvalid = A%sumIntoGlobalValues(gblrow,cols,vals)
!    end do ! iCell

! call A%fillComplete() !; FORTRILINOS_CHECK_IERR()

     ! C : Coefficient matrix -------------------------------------------------!

     call C%resumeFill() !; FORTRILINOS_CHECK_IERR()
     !C = TpetraCrsMatrix(graph)

     call mpas_timer_start("fort mat AllToScalar")
     call C%setAllToScalar(szero)
     call mpas_timer_stop("fort mat AllToScalar")


     do iCell = 1, nCellsArray(1)

        gblrow = globalIdx(iCell)
        row = iCell

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Interpolation sshEdge
          sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
          thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1),bottomDepth(cell2))
          fluxAx = edgeSignOnCell(i, iCell) * (thicknessSumLag / dcEdge(iEdge)) * dvEdge(iEdge)

!         if ( cell2 > 0 ) then
!           cols(1) = globalIdx(cell1) !map%getLocalElement((map%getGlobalElement(cell1)))
!           vals(1) = -fluxAx 
!           numvalid = C%sumIntoGlobalValues(gblrow,cols,vals)
  
            col(1) = cell1
            vals(1) = -fluxAx 
            numvalid = C%sumIntoLocalValues(row,col,vals)
  
            !cols(1) = globalIdx(cell2) !map%getLocalElement((map%getGlobalElement(cell1)))
            !vals(1) = +fluxAx 
            !numvalid = C%sumIntoGlobalValues(gblrow,cols,vals)
  
            col(1) = cell2
            vals(1) = +fluxAx 
            numvalid = C%sumIntoLocalValues(row,col,vals)
!         endif

        end do ! i

!       cols(1) = globalIdx(iCell) !map%getLocalElement((map%getGlobalElement(iCell)))
!       vals(1) = (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
!       numvalid = C%sumIntoGlobalValues(gblrow,cols,vals)

        col(1) = iCell !map%getLocalElement((map%getGlobalElement(iCell)))
        vals(1) = (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
        numvalid = C%sumIntoLocalValues(row,col,vals)

     end do ! iCell

  call mpas_timer_start("fort mat complete")
  call C%fillComplete() !; FORTRILINOS_CHECK_IERR()
  call mpas_timer_stop("fort mat complete")

  call mpas_timer_stop("fort mat setup")

! numrow_c = C%getNodeNumRows()
! numnnz_c = C%getNodeNumEntries()

! print*, numnnz_c
! print*, numrow_c

! call C%getAllValues(row_ptrs_c, colind_c, values_c)

! row_ptrs_c(:) = 0
! colind_c(:) = 0
! values_c(:) = 0.d0

! call mpi_allreduce(numnnz_c,numrow_c,1,mpi_integer,mpi_sum,dminfo%comm,mpi_ierr)

! print*, size(colind_c)
! print*, 'CCC'

! print*, 'AAA'
  call mpas_timer_start("fort map")
  !A = TpetraCrsMatrix(map,mape,max_entries_per_row, TpetraStaticProfile)
  !A = TpetraCrsMatrix(graph)
  !call A%replaceColMap(graph%getColMap())
  !call graph%replaceColMap

  !colmap = C%getColMap()
  !call A%replaceColMap(colmap)
  call mpas_timer_stop("fort map")
! call mpi_barrier(dminfo%comm,mpi_ierr)
! stop
! print*, 'DDD'

! call mpas_timer_start("fort setAllValues")
! call A%setAllValues(row_ptrs_c, colind_c, values_c)
! call mpas_timer_stop("fort setAllValues")

! print*, 'EEE'
!   call mpi_barrier(dminfo%comm,mpi_ierr)
!   stop 

! print*, 'AAA'
! call mpas_timer_start("fort expert")
! call A%expertStaticFillComplete(map,map)
! call mpas_timer_stop("fort expert")
  
! print*, 'BBB'

! call A%fillComplete() !; FORTRILINOS_CHECK_IERR()

! print*, 'CCC'
! numrow_c = A%getNodeNumEntries()

! print*, numnnz_c,numrow_c

! call A%getAllValues(row_ptrs_a, colind_a, values_a)
! numrow_a = A%getNodeNumRows()
! numnnz_a = A%getNodeNumEntries()

! print*, numrow_a- numrow_c, numnnz_a-numnnz_c
!  
  
! print*, C%isGloballyIndexed()
! print*, C%isLocallyIndexed()
! print*, minval(values_a),minval(values_c),maxval(values_a),maxval(values_c)

! description = C%description()

! print*, description
! call mpi_barrier(dminfo%comm,mpi_ierr)
! print*
! print*
! call mpi_barrier(dminfo%comm,mpi_ierr)
! description = A%description()
! print*, description

! print*, A%get
  
! sum1 = 0.d0
! do i = 1,numnnz_c
!   sum1 = sum1 + abs(values_a(i) - values_c(i)) 
! end do
! print*, sum1

! call mpi_barrier(dminfo%comm,mpi_ierr)
! stop

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

! print*, 'DDD'

  ! Step 2: setup the problem
  if ( stage == 'o' ) then
    call solver_handle_o%setup_matrix(C) !; FORTRILINOS_CHECK_IERR()
  elseif ( stage == 'm' ) then
    call solver_handle_m%setup_matrix(C) !; FORTRILINOS_CHECK_IERR()
  endif

! print*, 'EEE'

  ! Step 3: setup the solver - Initial only
  if ( init_belos_o .and. stage == 'o') then
    call solver_handle_o%setup_solver(plist_o) !; FORTRILINOS_CHECK_IERR()
    init_belos_o = .false.
!   init_belos_o = .true.
  elseif ( init_belos_m .and. stage == 'm') then
    call solver_handle_m%setup_solver(plist_m) !; FORTRILINOS_CHECK_IERR()
    init_belos_m = .false.
!   init_belos_m = .true.
  endif

! print*, 'FFF'

  call mpas_timer_stop("fort solver setup")

  call mpas_timer_start("fort init resid")
  ! Calculate initial residual
  call C%apply(X, residual, TeuchosNO_TRANS, sone, szero) !; FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone) !; FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms) !; FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
! print*, r0
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
  call C%apply(X, residual, TeuchosNO_TRANS, sone, szero) !; FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone) !; FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms) !; FORTRILINOS_CHECK_IERR()
  call mpas_timer_stop("fort check")

! print*, norms(1)


  if ( stage == 'o' .and. norms(1)/r0 > tol_o) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    !print*,r0,'tol = ',norms(1)/r0 
    stop 1
  elseif ( stage == 'm' .and. norms(1)/r0 > tol_m) then
    !print*,r0,'tol = ',norms(1)/r0 
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  end if

! r0 = norms(1)

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


! call mpas_timer_start("si halo ssh")
! call mpas_dmpar_exch_group_create(domain, iterGroupName)
! call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'sshSubcycle', 1 )
! call mpas_threading_barrier()
! call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
! call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
! call mpas_timer_stop("si halo ssh")
       
! call mpas_timer_start("fort release")
! call A%release()
! call mpas_timer_stop("fort release")

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

