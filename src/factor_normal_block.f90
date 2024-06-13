module factor_normal_block_class

  use global
  use indicators_dedic_class
  use covmat_block_invwishart_class
  use probability, only : rnorm
  use matrix,      only : matinv, chol
  implicit none


  type :: factor_normal_block
    integer               :: nobs
    integer               :: nmeas
    integer               :: nfac
    real(r8), allocatable :: theta(:,:)
    integer,  allocatable :: indi(:)
    integer,  allocatable :: mi(:)
    ! back up
    real(r8), allocatable :: theta_bak(:,:)
  end type factor_normal_block


contains


  !-----------------------------------------------------------------------------

  subroutine init_factor_normal_block(this, nobs, nmeas, nfac, start)
    implicit none
    type(factor_normal_block), intent(out) :: this
    integer,                   intent(in)  :: nobs
    integer,                   intent(in)  :: nmeas
    integer,                   intent(in)  :: nfac
    real(r8),                  intent(in)  :: start(nobs,nfac)
    integer                                :: i

    this%nobs  = nobs
    this%nmeas = nmeas
    this%nfac  = nfac

    allocate(this%theta(nobs, nfac))
    allocate(this%theta_bak(nobs, nfac))
    this%theta = start
    this%theta_bak = start

    allocate(this%indi(nfac))
    allocate(this%mi(nmeas))
    this%indi = (/(i, i=1, this%nfac)/)
    this%mi = (/(i, i=1, this%nmeas)/)

  end subroutine init_factor_normal_block


  !-----------------------------------------------------------------------------

  subroutine update_factor_normal_block_acti(this, Y, alpha, dedic, idioprec, fdist)
    implicit none
    type(factor_normal_block),     intent(inout) :: this
    real(r8),                      intent(in)    :: Y(this%nobs,this%nmeas)
    real(r8),                      intent(in)    :: idioprec(this%nmeas)
    real(r8),                      intent(in)    :: alpha(this%nmeas)
    type(indic_dedic),             intent(in)    :: dedic
    type(covmat_block_invwishart), intent(in)    :: fdist
    logical                                      :: acti(this%nfac)
    integer                                      :: actind(dedic%K1)
    real(r8)                                     :: mean_post(this%nobs,dedic%K1)
    real(r8)                                     :: var_post(dedic%K1,dedic%K1)
    integer                                      :: i, k

    if(dedic%K1 == 0) return    ! no active factors to sample

    acti   = dedic%active
    actind = pack(this%indi, acti)

    ! compute posterior moments for active factors
    var_post = matinv(fdist%var(actind,actind))
    do k = 1, dedic%K1
      var_post(k,k) = var_post(k,k) &
                    + sum(idioprec*(alpha**2), mask=(dedic%group==actind(k)))
      mean_post(:,k) = matmul(Y(:, pack(this%mi, dedic%group==actind(k))), &
                              pack(idioprec*alpha, dedic%group==actind(k)))
    end do
    var_post = matinv(var_post)

    ! sample active thiss factors
    do k = 1, this%nfac
      if(.not.acti(k)) cycle
      do i = 1, this%nobs
        this%theta(i,k) = rnorm()
      end do
    end do
    this%theta(:,actind) = matmul(mean_post, var_post) &
                           + matmul(this%theta(:,actind), &
                                    transpose(chol(var_post)))

  end subroutine update_factor_normal_block_acti


  !-----------------------------------------------------------------------------

  subroutine update_factor_normal_block_inacti(this, fdist, dedic)
    implicit none
    type(factor_normal_block),     intent(inout) :: this
    type(covmat_block_invwishart), intent(in)    :: fdist
    type(indic_dedic),             intent(in)    :: dedic
    integer                                      :: actind(dedic%K1)
    integer                                      :: inactind(dedic%K2)
    integer                                      :: i, k

    if(dedic%K2 == 0) return   ! no inactive factors to sample

    actind   = pack(this%indi, dedic%active)
    inactind = pack(this%indi, .not.dedic%active)

    do k = 1, dedic%K2
      do i = 1, this%nobs
        this%theta(i, inactind(k)) = rnorm()
      end do
    end do

    if(dedic%K1 == 0) then

      this%theta = matmul(this%theta, transpose(fdist%L221))

    else

      this%theta(:,inactind) = matmul(this%theta(:,inactind), &
                                        transpose(fdist%L221)) &
                               + matmul(this%theta(:,actind), fdist%G1112)

    end if

  end subroutine update_factor_normal_block_inacti


  !-----------------------------------------------------------------------------

  subroutine update_factor_normal_block(this, Y, alpha, dedic, idioprec, fdist)
    implicit none
    type(factor_normal_block),     intent(inout) :: this
    real(r8),                      intent(in)    :: Y(this%nobs, this%nmeas)
    real(r8),                      intent(in)    :: idioprec(this%nmeas)
    real(r8),                      intent(in)    :: alpha(this%nmeas)
    integer,                       intent(in)    :: dedic(this%nmeas)
    type(covmat_block_invwishart), intent(in)    :: fdist
    real(r8) :: ap2(this%nmeas)
    real(r8) :: ing(this%nmeas, this%nfac)
    real(r8) :: mean_post(this%nobs, this%nfac)
    real(r8) :: var_post(this%nfac, this%nfac)
    integer  :: i, j, k

    ! posterior covariance matrix
    var_post = fdist%prec
    ap2 = idioprec * (alpha**2)
    do k = 1, this%nfac
      var_post(k,k) = var_post(k,k) + sum(ap2, dedic==k)
    end do
    var_post = matinv(var_post)

    ! posterior mean
    do j = 1, this%nmeas
      if(dedic(j) == 0) then
        ing(j,:) = 0._r8
      else
        do k = 1, this%nfac
          ing(j,k) = var_post(dedic(j), k) * alpha(j) * idioprec(j)
        end do
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

  end subroutine update_factor_normal_block


  !-----------------------------------------------------------------------------

  subroutine backup_factor_normal_block(this)
    implicit none
    type(factor_normal_block), intent(inout) :: this

    this%theta_bak = this%theta

  end subroutine backup_factor_normal_block


  !-----------------------------------------------------------------------------

  subroutine restore_factor_normal_block(this)
    implicit none
    type(factor_normal_block), intent(inout):: this

    this%theta = this%theta_bak

  end subroutine restore_factor_normal_block


end module factor_normal_block_class

