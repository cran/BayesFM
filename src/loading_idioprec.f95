module loading_idioprec_class

  use global
  use probability, only : rnorm, rgamma
  implicit none

  private
  public :: loading_idioprec


  type :: loading_idioprec
    logical  :: cont            ! TRUE if continuous case, FALSE if categorical
    real(r8) :: alpha           ! factor loading
    real(r8) :: alpha_mu0       ! prior:
    real(r8) :: alpha_prec0     !   alpha ~ N(alpha_mu0, sig2/alpha_prec0)
    real(r8) :: var             ! idiosyncratic variance
    real(r8) :: prec            ! idiosyncratic precision
    real(r8) :: prec_a0         ! prior:
    real(r8) :: prec_b0         !   sig2  ~ IGamma(prec_a0, prec_b0)
    real(r8) :: prec_a_post     ! posterior shape for prec
    ! back up values
    real(r8) :: alpha_bak
    real(r8) :: var_bak
    real(r8) :: prec_bak
  contains
    procedure, public :: init    => init_loading_idioprec
    procedure, public :: update  => update_loading_idioprec
    procedure, public :: backup  => backup_loading_idioprec
    procedure, public :: restore => restore_loading_idioprec
  end type loading_idioprec


contains


  !-----------------------------------------------------------------------------

  subroutine init_loading_idioprec(this, nobs, cont, prior, start)
    implicit none
    class(loading_idioprec) :: this
    integer,  intent(in)    :: nobs
    logical,  intent(in)    :: cont
    real(r8), intent(in)    :: prior(3)
    real(r8), intent(in)    :: start(2)

    this%cont        = cont
    this%alpha_mu0   = 0._r8
    this%alpha_prec0 = prior(1)
    this%prec_a0     = prior(2)
    this%prec_b0     = prior(3)

    this%prec_a_post = this%prec_a0 + .5_r8*nobs

    this%alpha = start(1)
    this%var   = start(2)
    this%prec  = 1._r8/this%var

    this%alpha_bak = this%alpha
    this%var_bak   = this%var
    this%prec_bak  = this%prec

  end subroutine init_loading_idioprec


  !-----------------------------------------------------------------------------

  subroutine update_loading_idioprec(this, Yaux, dedic, fac)
    implicit none
    class(loading_idioprec)       :: this
    real(r8), intent(in)          :: Yaux(:)
    integer,  intent(in)          :: dedic
    real(r8), intent(in), target  :: fac(:,:)
    real(r8),             pointer :: fp(:) => null()
    real(r8)                      :: aN, AAN, prec_b_post

    !----- 'null model'
    if(dedic == 0) then

      if(.not.this%cont) return  ! nothing to do in categorical case

      ! sample idiosyncratic precision in continuous case
      prec_b_post = this%prec_b0 + .5_r8*sum(Yaux**2)
      this%prec = rgamma(this%prec_a_post, 1._r8/prec_b_post)
      this%var = 1._r8/this%prec

    !----- general case
    else

      fp => fac(:, dedic)

      aN  = dot_product(Yaux, fp)
      AAN = 1._r8/(sum(fp**2) + this%alpha_prec0)

      if(this%cont) then    ! sample precision in continuous case
        prec_b_post = this%prec_b0 + .5_r8*(sum(Yaux**2) - AAN*(aN**2))
        this%prec = rgamma(this%prec_a_post, 1._r8/prec_b_post)
        this%var = 1._r8/this%prec
      end if

      this%alpha = rnorm(AAn*aN, AAN*this%var)

    end if

  end subroutine update_loading_idioprec


  !-----------------------------------------------------------------------------

  subroutine backup_loading_idioprec(this)
    implicit none
    class(loading_idioprec) :: this

    this%alpha_bak = this%alpha
    this%var_bak   = this%var
    this%prec_bak  = this%prec

  end subroutine backup_loading_idioprec


  !-----------------------------------------------------------------------------

  subroutine restore_loading_idioprec(this)
    implicit none
    class(loading_idioprec) :: this

    this%alpha = this%alpha_bak
    this%var   = this%var_bak
    this%prec  = this%prec_bak

  end subroutine restore_loading_idioprec


end module loading_idioprec_class

