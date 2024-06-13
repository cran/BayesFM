module covariates_class

  use global
  use probability, only : rnorm
  use matrix,      only : chol, solvl, solvu
  implicit none


  type :: covariates
    integer               :: nobs
    integer               :: npar
    real(r8), allocatable :: beta(:)
    real(r8), allocatable :: X(:,:)
    real(r8), allocatable :: XX(:,:)
    real(r8), allocatable :: Xbeta(:)
    real(r8)              :: prec0
    ! back up
    real(r8), allocatable :: beta_bak(:)
    real(r8), allocatable :: Xbeta_bak(:)
  end type covariates


contains


  !-----------------------------------------------------------------------------

  subroutine init_covariates(this, nobs, npar, X, prior, start, mask)
    implicit none
    type(covariates), intent(out)          :: this
    integer,          intent(in)           :: nobs
    integer,          intent(in)           :: npar
    real(r8),         intent(in)           :: X(nobs, npar)
    real(r8),         intent(in)           :: prior
    real(r8),         intent(in)           :: start(npar)
    logical,          intent(in), optional :: mask(npar)
    integer                                :: i

    this%nobs = nobs

    if(present(mask)) then
      this%npar = count(mask)
    else
      this%npar = npar
    end if

    allocate(this%Xbeta(nobs))
    this%Xbeta = 0._r8

    if(this%npar == 0) return    ! nothing more to allocate if no this

    allocate(this%beta(this%npar))
    allocate(this%X(nobs, this%npar))
    allocate(this%XX(this%npar, this%npar))
    allocate(this%beta_bak(this%npar))
    allocate(this%Xbeta_bak(nobs))

    if(present(mask)) then
      do i = 1, nobs
        this%X(i,:) = pack(X(i,:), mask)
      end do
      this%beta = pack(start, mask)
    else
      this%X    = X
      this%beta = start
    end if

    this%XX    = matmul(transpose(this%X), this%X)
    this%Xbeta = matmul(this%X, this%beta)
    this%prec0 = prior

    this%beta_bak  = this%beta
    this%Xbeta_bak = this%Xbeta

  end subroutine init_covariates


  !-----------------------------------------------------------------------------

  subroutine update_covariates(this, Y, prec)
    implicit none
    type(covariates),                intent(inout) :: this
    real(r8), dimension(this%nobs), intent(in)    :: Y
    real(r8),                        intent(in)    :: prec
    real(r8), dimension(this%npar, this%npar)    :: beta_B, beta_L
    real(r8), dimension(this%npar)                :: beta_m, beta_z, beta_e
    integer                                        :: i

    if(this%npar == 0) return

    beta_m = prec * matmul(transpose(this%X), Y)
    beta_B = prec * this%XX
    do i = 1, this%npar
      beta_B(i,i) = beta_B(i,i) + this%prec0
    end do
    beta_L = chol(beta_B)
    beta_z = solvl(beta_L, beta_m)

    do i = 1, this%npar
      beta_e(i) = rnorm()
    end do
    this%beta  = solvu(transpose(beta_L), beta_e + beta_z)
    this%Xbeta = matmul(this%X, this%beta)

  end subroutine update_covariates


  !-----------------------------------------------------------------------------

  subroutine backup_covariates(this)
    implicit none
    type(covariates), intent(inout) :: this

    if(this%npar == 0) return

    this%beta_bak  = this%beta
    this%Xbeta_bak = this%Xbeta

  end subroutine backup_covariates


  !-----------------------------------------------------------------------------

  subroutine restore_covariates(this)
    implicit none
    type(covariates), intent(inout) :: this

    if(this%npar == 0) return

    this%beta  = this%beta_bak
    this%Xbeta = this%Xbeta_bak

  end subroutine restore_covariates


  !-----------------------------------------------------------------------------

  function get_all_covariates(X) result(par)
    implicit none
    type(covariates) :: X(:)
    real(r8)         :: par(sum(X%npar))
    integer          :: i, j

    i = 0
    do j = 1, size(X)
      if(X(j)%npar == 0) cycle
      par(i+1:i+X(j)%npar) = X(j)%beta
      i = i + X(j)%npar
    end do

  end function get_all_covariates


end module covariates_class

