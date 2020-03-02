module ocn_fortrilinos_imp_mod
!---------------------------
#include "ForTrilinosInterface_config.hpp"
#include "ForTrilinos.h"
  use forteuchos
  use fortpetra
  use fortrilinos
  use ocn_model_evaluator_mod
!---------------------------
  use, intrinsic :: iso_c_binding
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar
!---------------------------

  implicit none
  private  
  save


  public :: ocn_time_integration_imp_btrmode
! public :: noxinit, noxsolve, noxfinish
  public :: noxsolve,noxfinish
  class(ForModelEvaluator),private,allocatable :: model_evaluator

!*********************************************************************
                                contains
!*********************************************************************

  subroutine ocn_time_integration_imp_btrmode(domain,dt,xstate,nvec,n_tot_vec)
   
    implicit none
    !--------------
    type (TeuchosComm) :: comm
    type (ParameterList) :: params
    !--------------
    type (domain_type) :: domain
    type (dm_info) :: dminfo
    real (kind=RKIND)  :: dt
    real (c_double), dimension(nvec) :: xstate
    integer (c_int) :: ierr,nvec,n_tot_vec
    integer(global_size_type) :: n_global_vec,n_global2_vec,n_total_vec
    integer(global_size_type) :: nCells1_vec,nCells4_vec
    logical :: init_nox = .true.
    ! -------------------------------------------------------------------
   
    n_global_vec = nvec
     n_total_vec = n_tot_vec

    !--- Only for initial
    if ( init_nox ) then
   
!     print*, 'PRINT in init'

      init_nox = .false.

      dminfo = domain % dminfo

      comm = TeuchosComm(dminfo % comm)
      params = ParameterList("ocnModelEvaluator")
      call load_from_xml(params,'nox_params.xml')

      allocate(model_evaluator, source=ocnModelEvaluator(comm,domain,xstate,dt,n_global_vec,n_total_vec))
      call init_ForModelEvaluator(model_evaluator)
      call model_evaluator%setup(params)
      call params%release()
    endif ! init_nox
    

    call noxsolve(comm,domain,dt,params,xstate,nvec,ierr)

  end subroutine ocn_time_integration_imp_btrmode

!------------------------------------------------------------------------

  subroutine noxsolve(comm,domain,dt,params,xstate,nvec,ierr)
    ! -------------------------------------------------------------------
    use,intrinsic :: iso_c_binding
    implicit none

    type (TeuchosComm) :: comm
    type (domain_type) :: domain
    real (kind=RKIND)  :: dt
    integer(c_int) :: vector_size,ierr,nvec
    real(c_double),dimension(nvec) :: xstate

    type(NOXSolver) :: nox_solver
    type(ocnModelEvaluator) :: self
    integer(kind(NOXStatusType)) :: status
    type(ParameterList) :: params
    integer(global_size_type) :: n_global_vec
    ! -------------------------------------------------------------------

    ierr = 0

    n_global_vec = nvec

!     print*, 'PRINT in noxsolve'
    
    ! ------------------------------------------------- !
    call mpas_timer_start("si fortrilinos params")
    params = ParameterList("ocnModelEvaluator")
    call load_from_xml(params,'nox_params.xml')
    call mpas_timer_stop ("si fortrilinos params")
    ! ------------------------------------------------- !

    ! ------------------------------------------------- !
    call mpas_timer_start("si fortrilinos evalse")
    call model_evaluator%setup(params)
    call mpas_timer_stop ("si fortrilinos evalse")
    ! ------------------------------------------------- !

    nox_solver = NOXSolver(model_evaluator)

    ! ------------------------------------------------- !
    call mpas_timer_start("si fortrilinos nox setup")
    call nox_solver%setup(params)
    call mpas_timer_stop ("si fortrilinos nox setup")
    ! ------------------------------------------------- !

    ! ------------------------------------------------- !
    call mpas_timer_start("si fortrilinos nox solver")
    status = nox_solver%solve()
    call mpas_timer_stop ("si fortrilinos nox solver")
    ! ------------------------------------------------- !

    ! ------------------------------------------------- !
    call mpas_timer_start("si fortrilinos nox release")
    call nox_solver%release()
    call mpas_timer_stop ("si fortrilinos nox release")
    ! ------------------------------------------------- !

    if (status /= NOXConverged) ierr = 1
!   print*, "Status@@",status,NOXConverged

    !--- Get the final solution
    call getFinalSolution(self,xstate,n_global_vec)

    call params%release()

  end subroutine noxsolve

!------------------------------------------------------------------------

  subroutine noxfinish()
    call model_evaluator%release()
    deallocate(model_evaluator)
  end subroutine noxfinish

end module
