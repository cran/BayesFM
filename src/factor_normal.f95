module factor_normal_class

  use global
  use covmat_block_invwishart_class
  use probability, only : rnorm
  use matrix,      only : matinv, chol
  implicit none

  private
  public :: factor_normal


  type :: factor_normal
    integer               :: nobs
    integer               :: nmeas
    integer               :: nfac
    real(r8), allocatable :: theta(:,:)
    ! back up
    real(r8), allocatable :: theta_bak(:,:)
  contains
    procedure, public :: init    => init_factor_normal
    procedure, public :: update  => update_factor_normal
    procedure, public :: backup  => backup_factor_normal
    procedure, public :: restore => restore_factor_normal
  end type factor_normal


contains


  !-----------------------------------------------------------------------------

  subroutine init_factor_normal(this, nobs, nmeas, nfac, start)
    implicit none
    class(factor_normal) :: this
    integer,  intent(in) :: nobs
    integer,  intent(in) :: nmeas
    integer,  intent(in) :: nfac
    real(r8), intent(in) :: start(nobs,nfac)

    this%nobs  = nobs
    this%nmeas = nmeas
    this%nfac  = nfac

    allocate(this%theta(nobs, nfac))
    allocate(this%theta_bak(nobs, nfac))
    this%theta = start
    this%theta_bak = start

  end subroutine init_factor_normal


  !-----------------------------------------------------------------------------

  subroutine update_factor_normal(this, Y, alpha, dedic, idioprec, fdist)
    implicit none
    class(factor_normal)                       :: this
    real(r8),                       intent(in) :: Y(this%nobs, this%nmeas)
    real(r8),                       intent(in) :: idioprec(this%nmeas)
    real(r8),                       intent(in) :: alpha(this%nmeas)
    integer,                        intent(in) :: dedic(this%nmeas)
    class(covmat_block_invwishart), intent(in) :: fdist
    real(r8) :: ap2(this%nmeas)
    real(r8) :: ing(this%nmeas, this%nfac)
    real(r8) :: mean_post(this%nobs, this%nfac)
    real(r8) :: var_post(this%nfac, this%nfac)
    integer  :: i, j, k

    ! posterior covariance matrix
    var_post = fdist%prec
    ap2 = idioprec * (alpha**2)
    forall(k = 1:this%nfac) var_post(k,k) = var_post(k,k) + sum(ap2, dedic==k)
    var_post = matinv(var_post)

    ! posterior mean
    do j = 1, this%nmeas
      if(dedic(j) == 0) then
        ing(j,:) = 0._r8
      else
        forall(k = 1:this%nfac)
          ing(j,k) = var_post(dedic(j), k) * alpha(j) * idioprec(j)
        end forall
      end if
    end do
    mean_post = matmul(Y, ing)

    ! sample factors
    do k = 1, this%nfac
      do i = 1, this%nobs
        this%theta(i,k) = rnorm()
      end do
    end do
    this%theta = mean_post + matmul(this%theta, transpose(chol(var_post)))

  end subroutine update_factor_normal


  !-----------------------------------------------------------------------------

  subroutine backup_factor_normal(this)
    implicit none
    class(factor_normal) :: this

    this%theta_bak = this%theta

  end subroutine backup_factor_normal


  !-----------------------------------------------------------------------------

  subroutine restore_factor_normal(this)
    implicit none
    class(factor_normal) :: this

    this%theta = this%theta_bak

  end subroutine restore_factor_normal


end module factor_normal_class

