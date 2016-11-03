module mcmc_progress_class

  implicit none

  type :: mcmc_progress
    logical      :: verbose
    integer      :: burnin
    integer      :: i
    integer      :: steps(20)
    character(6) :: perct(20)
  contains
    procedure, public :: init => init_mcmc_progress
    procedure, public :: show => show_mcmc_progress
  end type mcmc_progress


contains


  !-----------------------------------------------------------------------------

  subroutine init_mcmc_progress(this, burnin, iter, verbose)
    implicit none
    class(mcmc_progress) :: this
    integer, intent(in)  :: iter, burnin
    logical, intent(in)  :: verbose
    integer              :: i

    this%verbose = verbose
    this%burnin = burnin
    this%i = 1

    this%steps = (/(i, i=iter/20, iter, iter/20)/)
    this%steps(20) = iter   ! make sure "100%" displayed after last iteration

    this%perct = ["    5%", "   10%", "   15%", "   20%", "   25%", "   30%", &
                  "   35%", "   40%", "   45%", "   50%", "   55%", "   60%", &
                  "   65%", "   70%", "   75%", "   80%", "   85%", "   90%", &
                  "   95%", "  100%"]

  end subroutine init_mcmc_progress


  !-----------------------------------------------------------------------------

  subroutine show_mcmc_progress(this, rep)
    implicit none
    class(mcmc_progress) :: this
    integer, intent(in)  :: rep

    if(modulo(rep, 100) == 0) call rchkusr()

    if(this%verbose) then
      if(rep == this%burnin) then
        call intpr("done with burn-in period", 24, 0, 0)
      end if
      if(rep == this%steps(this%i)) then
        call intpr(this%perct(this%i), 6, 0, 0)
        this%i = this%i+1
      end if
    end if

  end subroutine show_mcmc_progress


end module mcmc_progress_class

