module ocn_belos_mod
!********************************************************8
!********************************************************8

! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
!module myoperators
  use forteuchos
  use fortpetra
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar

  use mpas_timer
  use mpas_threading
  use mpas_timekeeping
  use mpas_log

  implicit none

  type, extends(ForTpetraOperator) :: BtrOperator
    type(TpetraMap) :: row_map, col_map, domain_map, range_map
  contains
!   procedure :: apply => my_apply
!   procedure :: getDomainMap => my_getDomainMap
!   procedure :: getRangeMap => my_getRangeMap
    procedure :: release => delete_BtrOperator
  end type
! interface BtrOperator
!   procedure new_BtrOperator
! end interface

contains
! function new_BtrOperator(row_map, col_map) &
!     result(self)
!   use, intrinsic :: ISO_C_BINDING
!   type(TpetraMap), intent(in) :: row_map, col_map
!   type(BtrOperator) :: self

!   self = ForTpetraOperator()
!   self%row_map = row_map
!   self%col_map = col_map
!   self%domain_map = row_map
!   self%range_map = row_map
! end function

! subroutine my_apply(self, x, y, mode, alpha, beta)
!   use, intrinsic :: ISO_C_BINDING
!   implicit none
!   class(BtrOperator), intent(in) :: self
!   class(TpetraMultiVector), intent(in) :: x
!   class(TpetraMultiVector), intent(in) :: y
!   integer(kind(TeuchosETransp)), intent(in) :: mode
!   real(scalar_type), intent(in) :: alpha
!   real(scalar_type), intent(in) :: beta
!   integer :: lid
!   integer(global_ordinal_type) :: gid

!   type(TeuchosComm) :: comm
!   type(TpetraImport) :: import
!   type(TpetraMultiVector) :: x_ghosted
!   integer :: my_rank, num_procs
!   integer :: i, n
!   real(scalar_type), dimension(:), pointer :: xdata
!   real(scalar_type), dimension(:), pointer :: ydata

!   integer, save :: counter = 0

!   counter = counter + 1

!   comm = self%row_map%getComm()
!   my_rank = comm%getRank()
!   num_procs = comm%getSize()

!   import = TpetraImport(self%domain_map, self%col_map)
!   x_ghosted = TpetraMultiVector(self%col_map, 1_size_type)
!   call x_ghosted%doImport(x, import, TpetraINSERT)
!   call import%release()

!   xdata => x_ghosted%getData        (1_size_type)
!   ydata => y        %getDataNonConst(1_size_type)
!   n = y%getLocalLength()

!   ! Sometimes, ydata may be unitialized (when beta is 0), potentially containing
!   ! signaling NaNs. Therefore, for beta = 0, we explicitly zero it out.
!   if (beta .eq. 0) then
!     do i = 1, n
!       ydata(i) = 0
!     end do
!   else
!     do i = 1, n
!       ydata(i) = beta * ydata(i)
!     end do
!   end if

!   ! y = alpha * A*x + beta * y
!   do i = 1, n
!     gid = self%range_map%getGlobalElement(i)

!     ! A has [-1 2 -1] stencil
!     if (i > 1 .or. my_rank > 0) then
!       lid = self%col_map%getLocalElement(gid-1)
!       ydata(i) = ydata(i) - alpha*xdata(lid)
!     end if

!     lid = self%col_map%getLocalElement(gid)
!     ydata(i) = ydata(i) + 2*alpha*xdata(lid)

!     if (i < n .or. my_rank .ne. num_procs-1) then
!       lid = self%col_map%getLocalElement(gid+1)
!       ydata(i) = ydata(i) - alpha*xdata(lid)
!     end if
!   end do

!   ! Residual --------------------------------------------------------

!   ! -----------------------------------------------------------------

!   nullify(xdata)
!   nullify(ydata)

!   call x_ghosted%release()
!   call comm%release()

! end subroutine

! function my_getDomainMap(self) &
!     result(domain_map)
!   use, intrinsic :: ISO_C_BINDING
!   implicit none
!   class(BtrOperator), intent(in) :: self
!   type(TpetraMap) :: domain_map

!   domain_map = self%domain_map
! end function

! function my_getRangeMap(self) &
!     result(range_map)
!   use, intrinsic :: ISO_C_BINDING
!   implicit none
!   class(BtrOperator), intent(in) :: self
!   type(TpetraMap) :: range_map

!   range_map = self%range_map
! end function

  subroutine delete_BtrOperator(self)
    class(BtrOperator), intent(inout) :: self

    call self%row_map%release()
    call self%col_map%release()
    call self%domain_map%release()
    call self%range_map%release()

#ifdef __GNUC__
    ! FIXME This segfaults with Flang
    ! Call base class release()
    call self%ForTpetraOperator%release()
#endif
  end subroutine


end module

!********************************************************8
!********************************************************8


!!---------------------------
!  use forteuchos
!  use fortpetra
!  use fortrilinos
!! use ocn_residual_mod
!!---------------------------
!  use, intrinsic :: iso_c_binding
!  use mpas_derived_types
!  use mpas_pool_routines
!  use mpas_constants
!  use mpas_dmpar
!
!  use mpas_timer
!  use mpas_threading
!  use mpas_timekeeping
!  use mpas_log
!!---------------------------
!
!  implicit none
!
!  real(scalar_type), parameter :: zero=0., one=1.
!  logical,save :: justDoIt  = .true.
!
!  character (len=*), parameter :: iterGroupName = 'iterFields'
!
!  type(TpetraMultiVector),private,allocatable :: solnvec
!  type, extends(ForModelEvaluator) :: ocnModelEvaluator
!    type(TeuchosComm), private :: comm
!    type(TpetraMap), private :: x_map, f_map
!    type(TpetraImport), private :: importer
!    type(TpetraMultiVector), private :: x,f
!    ! -------------------------------------------------------
!    ! For MPASO 
!    type(domain_type) :: domain
!    integer :: nCells,nEdges,nCells1,nCells2,nCells4
!    integer, dimension(:), allocatable   :: nCellsArray, nEdgesArray
!    integer, dimension(:), allocatable   :: nEdgesOnEdge, nEdgesOnCell,indexToCellID
!    integer, dimension(:,:), allocatable :: cellsOnEdge,edgesOnCell, edgeSignOnCell
!    real (kind=RKIND) :: dt
!    real (kind=RKIND), dimension(:), allocatable :: dcEdge, bottomDepth, dvEdge, areaCell
!    real (kind=RKIND), dimension(:), allocatable :: sshCur,sshSubcycleCur,normalBarotropicVelocityCur
!    real (kind=RKIND), dimension(:), allocatable :: barotropicForcing,barotropicCoriolisTerm
!    real (kind=RKIND), dimension(:), allocatable :: CGvec_r0
!    integer :: s1_nCellsArray,s1_nEdgesArray
!    integer :: s1_nEdgesOnCell
!    integer :: s1_cellsOnEdge, s2_cellsOnEdge, s1_edgesOnCell, s2_edgesOncell
!    integer :: s1_edgeSignOnCell, s2_edgeSignOnCell
!    integer :: s1_dcEdge,s1_bottomDepth,s1_dvEdge,s1_areaCell
!    integer :: s1_sshCur,s1_sshSubcycleCur,s1_normalBarotropicVelocityCur
!    integer :: s1_barotropicForcing, s1_barotropicCoriolisTerm,s1_CGvec_r0
!    ! -------------------------------------------------------
!   
!  contains
!    procedure, private :: getFinalSolution
!    procedure :: evaluate_residual       => ocnModelEvaluator_eval_resid
!    procedure :: update_solution_vector  => ocnModelEvaluator_update_solution_vector
!    procedure :: get_x_map               => ocnModelEvaluator_get_x_map 
!    procedure :: get_f_map               => ocnModelEvaluator_get_f_map 
!    procedure :: release                 => delete_ocnModelEvaluator
!  end type
!
!  interface ocnModelEvaluator
!    procedure new_ocnModelEvaluator
!  end interface 
!
!!*********************************************************************
!                                contains
!!*********************************************************************
!
!  function new_ocnModelEvaluator(comm,domain,xstate,dt,n_global_vec,n_total_vec) result(self)
!    ! -------------------------------------------------------------- !
!    use,intrinsic :: iso_c_binding
!    type(domain_type) :: domain
!    real(kind=RKIND) :: dt
!    type(ocnMOdelEvaluator) :: self
!    type(TeuchosComm), intent(in) :: comm 
!    integer(global_size_type), intent(in) :: n_global_vec,n_total_vec
!    integer :: i,j
!    integer(size_type) :: num_vecs=1
!    integer(size_type) :: col
!    integer :: lclrow
!    real(scalar_type) :: val
!    ! -------------------------------------------------------------- !
!    type (block_type), pointer :: block
!    type (mpas_pool_type), pointer :: statePool
!    type (mpas_pool_type), pointer :: tracersPool
!    type (mpas_pool_type), pointer :: meshPool
!    type (mpas_pool_type), pointer :: diagnosticsPool
!    type (mpas_pool_type), pointer :: tendPool
!    type (mpas_pool_type), pointer :: tracersTendPool
!
!    integer :: nCells, nEdges
!    integer :: iCell,iEdge,cell1,cell2
!    integer, pointer :: nCellsPtr,nEdgesPtr
!    integer, dimension(:), pointer :: nCellsArray, nEdgesArray
!    integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell, indexToCellID
!    integer, dimension(:,:), pointer :: cellsOnEdge,edgesOnCell, edgeSignOnCell
!    real (kind=RKIND), dimension(:), pointer :: dcEdge, bottomDepth, dvEdge, areaCell
!    real (kind=RKIND), dimension(:), pointer :: sshCur,sshSubcycleCur,normalBarotropicVelocityCur
!    real (kind=RKIND), dimension(:), pointer :: barotropicForcing,barotropicCoriolisTerm
!    real (kind=RKIND), dimension(:), pointer :: CGvec_r0
!
!    real (kind=RKIND), dimension(n_global_vec) :: xstate
!
!    integer(global_ordinal_type),dimension(n_global_vec) :: elements
!    ! -------------------------------------------------------------- !
!    real(scalar_type), dimension(:), pointer :: xdata
!    integer(size_type), parameter :: ione = 1
!    integer(size_type) :: nState
!    integer :: itmp1,itmp2
!    real (kind=RKIND) :: temp1,temp2
!    integer(global_ordinal_type),dimension(:),allocatable :: iTmpArray
!    ! -------------------------------------------------------------- !
!
!!   print*, 'PRINT in new_'
!   
!    self = ForModelEvaluator()
!    self%comm = comm
!    self%domain = domain
!    self%dt = dt
!
!    block => domain % blocklist
!    do while (associated(block))
!
!      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
!      call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
!      call mpas_pool_get_subpool(tendPool       , 'tracersTend', tracersTendPool)
!      call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
!      call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
!      call mpas_pool_get_subpool(statePool      , 'tracers'    , tracersPool    )
!      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
!
!      call mpas_pool_get_array(meshPool, 'nEdgesOnCell'  , nEdgesOnCell  )
!      call mpas_pool_get_array(meshPool, 'edgesOnCell'  , edgesOnCell   )
!      call mpas_pool_get_array(meshPool, 'cellsOnEdge'  , cellsOnEdge   )
!      call mpas_pool_get_array(meshPool, 'indexToCellID'   , indexToCellID   )
!      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur,1)
!      call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
!
!      nCells = nCellsArray(1)
!      self%nCells = nCells
!      self%nCells2 = nCellsArray(2)
!      self%nCells4 = nCellsArray(4)
!
!      do iCell = 1,nCells
!        elements(iCell) = indexToCellID(iCell)
!      end do
!
!      ! owned_space
!      self%x_map = TpetraMap(n_total_vec, elements, comm)
!
!      block => block % next
!    end do  ! block
!
!    !-----------------------------------------------------------------------------!
!
!    self%importer = TpetraImport(self%x_map,self%x_map)
!
!    ! residual space
!    self%f_map = self%x_map
!
!    self%x  = TpetraMultiVector(self%x_map,num_vecs)
!
!    ! -------------------------------------------------------------- !  ! For MPASO variables input
! 
!  end function new_ocnModelEvaluator
!
!  !===================================================================
!
!  subroutine ocnModelEvaluator_eval_resid(self,x,f)
!    
!    !-----------------------------------------------------------------
!    class(ocnModelEvaluator),intent(in) :: self
!    class(TpetraMultiVector),intent(in) :: x
!    class(TpetraMultiVector),intent(in) :: f
!    integer(size_type), parameter :: ione = 1
!    real(scalar_type), dimension(:), pointer :: xdata
!    real(scalar_type), dimension(:), pointer :: udata
!    real(scalar_type), dimension(:), pointer :: fdata
!    integer(c_int)                        :: nvec
!    integer(size_type) :: col
!    integer :: lclrow
!    real(scalar_type) :: val
!    !-------
!!   type(domain_type) :: domain
!    type(MPAS_TimeInterval_type) :: timeStep
!    real(kind=RKIND) :: dt
!    integer :: ierr,err_tmp
!    type (block_type),pointer :: block
!    type (mpas_pool_type), pointer :: statePool
!    type (mpas_pool_type), pointer :: tracersPool
!    type (mpas_pool_type), pointer :: meshPool
!    type (mpas_pool_type), pointer :: diagnosticsPool
!    type (mpas_pool_type), pointer :: tendPool
!    type (mpas_pool_type), pointer :: tracersTendPool
!
!    integer :: nCells, nEdges
!    integer :: s1_nCellsArray,s1_nEdgesArray
!    integer :: s1_nEdgesOnCell,s1_dvEdge,s1_dcEdge,s1_areaCell,s1_bottomDepth
!    integer :: s1_cellsOnEdge,s2_cellsOnEdge,s1_edgesOnCell,s2_edgesOnCell
!    integer :: s1_edgeSignOnCell,s2_edgeSignOnCell
!    integer :: s1_sshCur,s1_sshSubcycleCur,s1_normalBarotropicVelocityCur
!    integer :: s1_barotropicForcing, s1_barotropicCoriolisTerm,s1_CGvec_r0
!    integer :: i,j,iCell,iEdge,cell1,cell2
!
!    integer,pointer :: nCellsPtr, nEdgesPtr
!    integer, dimension(:),pointer :: nCellsArray
!    integer, dimension(:),pointer :: nEdgesArray
!    integer, dimension(:),pointer :: nEdgesOnCell
!    integer, dimension(:,:),pointer :: cellsOnEdge
!    integer, dimension(:,:),pointer :: edgesOnCell
!    integer, dimension(:,:),pointer :: edgeSignOnCell
!    real (kind=RKIND) :: sshTendb1,sshTendb2,sshTendAx,sshEdgeCur,thicknessSumCur,sshDiffCur
!    real (kind=RKIND) :: sshEdgeMid,sshEdgeLag,thicknessSumMid,thicknessSumLag,sshDiffLag
!    real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshCurArea,sshLagArea
!    real (kind=RKIND), dimension(:),pointer :: dcEdge
!    real (kind=RKIND), dimension(:),pointer :: bottomDepth
!    real (kind=RKIND), dimension(:),pointer :: dvEdge
!    real (kind=RKIND), dimension(:),pointer :: areaCell
!
!    real (kind=RKIND), dimension(:),pointer :: sshCur
!    real (kind=RKIND), dimension(:),pointer :: sshSubcycleCur
!    real (kind=RKIND), dimension(:),pointer :: normalBarotropicVelocityCur
!    real (kind=RKIND), dimension(:),pointer :: barotropicForcing
!    real (kind=RKIND), dimension(:),pointer :: barotropicCoriolisTerm
!    real (kind=RKIND), dimension(:),pointer :: CGvec_r0,CGvec_r1
!
!    character (len=*), parameter :: iterGroupName = 'iterFields'
!    !-----------------------------------------------------------------
!
!    call mpas_timer_start("si fortrilinos residu")
!
!    call self%update_solution_vector(x)
!
!    call f%putScalar(zero)
!
!    !**************************************************************!
! 
!    block => self%domain % blocklist
!    do while (associated(block))
!      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
!      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
!      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
!      call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
!
!      call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
!      call mpas_pool_get_subpool(tendPool       , 'tracersTend', tracersTendPool)
!      call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
!      call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
!      call mpas_pool_get_subpool(statePool      , 'tracers'    , tracersPool    )
!      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
!
!      call mpas_pool_get_array(meshPool, 'nEdgesOnCell'  , nEdgesOnCell  )
!      call mpas_pool_get_array(meshPool, 'edgesOnCell'   , edgesOnCell   )
!      call mpas_pool_get_array(meshPool, 'cellsOnEdge'   , cellsOnEdge   )
!      call mpas_pool_get_array(meshPool, 'dcEdge'        , dcEdge        )
!      call mpas_pool_get_array(meshPool, 'bottomDepth'   , bottomDepth   )
!      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
!      call mpas_pool_get_array(meshPool, 'dvEdge'        , dvEdge        )
!      call mpas_pool_get_array(meshPool, 'areaCell'      , areaCell      )
!
!      call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
!      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur,1)
!      call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
!
!      call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
!      call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
!      call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
!      call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1', CGvec_r1)
!
!      nCells = nCellsArray( 1 )
!      nEdges = nEdgesArray( 2 )
!
!      ! Initial value change ------------ one time per timestep
!
!      if ( justDoIt ) then
!        justDoIt = .false.
!  
!        do iCell = 1,self%nCells
!          lclrow = iCell
!          col = 1
!          val = sshSubcycleCur(iCell)
!          call x%replaceLocalValue(lclrow,col,val)
!          call solnvec%replaceLocalValue(lclrow,col,val)
!        end do
!  
!      endif ! justDoIt
!  
!      udata => solnvec%getData(ione)
!
!      CGvec_r1(1:nCells) = udata(1:nCells)
!
!
!      ! Halo exchange 
!      call ocnHaloExchange (self)
!
!        ! SpMV : A*x --------------------------------------------------------------------
!  
!        do iCell = 1, nCells
!  
!          sshTendb1 = 0.0_RKIND
!          sshTendb2 = 0.0_RKIND
!          sshTendAx = 0.0_RKIND
!  
!          do i = 1, nEdgesOnCell(iCell)
!            iEdge = edgesOnCell(i, iCell)
!  
!            cell1 = cellsOnEdge(1, iEdge)
!            cell2 = cellsOnEdge(2, iEdge)
!           
!            ! Interpolation sshEdge
!            sshEdgeCur = 0.5_RKIND * (  sshCur(cell1) +   sshCur(cell2))
!            sshEdgeLag = 0.5_RKIND * (CGvec_r1(cell1) + CGvec_r1(cell2))
!            sshEdgeMid = 0.5_RKIND * (sshEdgeLag + sshEdgeCur)
!  
!            ! method 1, matches method 0 without pbcs, works with pbcs.
!            thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
!            thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
!            thicknessSumMid = sshEdgeMid + min(bottomDepth(cell1), bottomDepth(cell2))
!  
!            ! nabla (ssh^0)
!            sshDiffCur = (  sshCur(cell2) -    sshCur(cell1)) / dcEdge(iEdge)
!            sshDiffLag = (CGvec_r1(cell2) -  CGvec_r1(cell1)) / dcEdge(iEdge)
!
!            !--------------------------------------------------------------!
!            fluxb1 = thicknessSumMid * normalBarotropicVelocityCur(iEdge)
!            fluxb2 = thicknessSumLag * (0.5_RKIND*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)) )
!            fluxAx = thicknessSumLag * sshDiffLag
!  
!            sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
!            sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
!            sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
!
!            !--------------------------------------------------------------!
!  
!          end do ! i
!         
!          !-------------------------------------------------------------------------------! 
!          sshTendb1 = (4.0_RKIND/(gravity*self%dt)) * sshTendb1
!          sshTendb2 = (2.0_RKIND/(gravity        )) * sshTendb2
!          sshTendAx =                                 sshTendAx
!
!          sshCurArea = (4.0_RKIND/(gravity*self%dt**2.0)) * sshCur(iCell) * areaCell(iCell)
!          sshLagArea = (4.0_RKIND/(gravity*self%dt**2.0)) *  udata(iCell) * areaCell(iCell)
!  
!          CGvec_r0(iCell) =+(-sshCurArea - sshTendb1 + sshTendb2)   &
!                           -(-sshLagArea - sshTendAx) 
!          !-------------------------------------------------------------------------------! 
!  
!          lclrow = iCell
!          col = 1
!          val = CGvec_r0(iCell)
!  
!          ! Replace value ------------------------------------------------ 
!          call f%replaceLocalValue(lclrow,col,val)
!  
!        end do ! iCell
!  
!      ! -------------------------------------------------------------------------------
!
!      block => block % next
!    end do  ! block
!
!    !**************************************************************!
!
!    nullify(udata)
!
!!   print*, 'PRINT in resid end'
!
!    call mpas_timer_stop ("si fortrilinos residu")
!
!  end subroutine ocnModelEvaluator_eval_resid
!
!  !===================================================================
!
!  subroutine ocnModelEvaluator_update_solution_vector(self,xp)
!
!    !----------------------------------------------------------------
!    class(ocnModelEvaluator),intent(in) :: self
!    class(TpetraMultiVector),intent(in) :: xp
!    integer(size_type) :: num_vecs=1
!    !----------------------------------------------------------------
!
!    if (.not.allocated(solnvec)) then
!      allocate(solnvec, source=TpetraMultiVector(self%x_map,num_vecs))
!    endif
!
!    call solnvec%doImport(xp, self%importer, TpetraREPLACE)
!
!  end subroutine ocnModelEvaluator_update_solution_vector
!
!  !==================================================================
!
!  function ocnModelEvaluator_get_x_map(self) result(map)
!    class(ocnModelEvaluator),intent(in) :: self
!    type(TpetraMap) :: map
!
!    map = self%x_map
!
!  end function ocnModelEvaluator_get_x_map
!
!  !==================================================================
!
!  function ocnModelEvaluator_get_f_map(self) result(map)
!    class(ocnModelEvaluator),intent(in) :: self
!    type(TpetraMap) :: map
!
!    map = self%x_map
!
!  end function ocnModelEvaluator_get_f_map
!
!  !==================================================================
!
!  subroutine getFinalSolution (self,fvec,n_global_vec)
!
!    class(ocnModelEvaluator),intent(in) :: self
!    integer(global_size_type), intent(in) :: n_global_vec
!    real (kind=RKIND), dimension(n_global_vec) :: fvec
!    real(scalar_type), dimension(:), pointer :: udata
!    integer(size_type), parameter :: ione = 1
!   
!    ! This is for the initial value replacement in 'evaluate_resid'.
!    justDoIt = .true.
!
!    udata => solnvec%getData(ione)
!    fvec(:) = udata(:)
!
!  end subroutine getFinalSolution
!
!  !==================================================================
!
!  subroutine ocnHaloExchange (self)
!
!    ! ----------------------------------------------------
!    class(ocnModelEvaluator) :: self
!    real(scalar_type), dimension(:), pointer :: udata
!    integer(size_type), parameter :: ione = 1
!    ! ----------------------------------------------------
!    type (block_type), pointer :: block
!    type (mpas_pool_type), pointer :: statePool
!    type (mpas_pool_type), pointer :: tracersPool
!    type (mpas_pool_type), pointer :: meshPool
!    type (mpas_pool_type), pointer :: diagnosticsPool
!    type (mpas_pool_type), pointer :: tendPool
!    type (mpas_pool_type), pointer :: tracersTendPool
!    ! ----------------------------------------------------
!    integer :: nCells, nEdges
!    integer :: iCell,iEdge,cell1,cell2
!    integer, pointer :: nCellsPtr,nEdgesPtr
!    integer, dimension(:), pointer :: nCellsArray, nEdgesArray
!    real (kind=RKIND), dimension(:), pointer :: CGvec_r1
!    ! ----------------------------------------------------
!
!    block => self%domain % blocklist
!    do while (associated(block))
!      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
!      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
!      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
!      call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
!
!      call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
!      call mpas_pool_get_subpool(tendPool       , 'tracersTend', tracersTendPool)
!      call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
!      call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
!      call mpas_pool_get_subpool(statePool      , 'tracers'    , tracersPool    )
!      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
!      call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1', CGvec_r1)
!
!      call mpas_threading_barrier()
!      call mpas_timer_start("si halo iter r0")
!      call mpas_dmpar_exch_group_create(self%domain, iterGroupName)
!      call mpas_dmpar_exch_group_add_field(self%domain, iterGroupName, 'CGvec_r1', 1)
!      call mpas_threading_barrier()
!      call mpas_dmpar_exch_group_full_halo_exch(self%domain, iterGroupName)
!      call mpas_dmpar_exch_group_destroy(self%domain, iterGroupName)
!      call mpas_timer_stop("si halo iter r0")
!
!      block => block % next
!    end do  ! block
!   
!  end subroutine ocnHaloExchange
!
!  !==================================================================
!
!  subroutine delete_ocnModelEvaluator(self)
!    class(ocnModelEvaluator),intent(inout) :: self
!
!    call self%comm%release()
!    call self%x_map%release()
!    call self%f_map%release()
!    call self%importer%release()
!    call self%x%release()
!    call solnvec%release()
!
!    deallocate(solnvec)
!   
!    call self%ForModelEvaluator%release()
!
!  end subroutine delete_ocnModelEvaluator
!
!  !==================================================================
!
!end module
