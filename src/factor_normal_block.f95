module factor_normal_block_class

  use global
  use factor_normal_class
  use indicators_dedic_class
  use covmat_block_invwishart_class
  use probability, only : rnorm
  use matrix,      only : matinv, chol
  implicit none

  private
  public :: factor_normal_block


  type, extends(factor_normal) :: factor_normal_block
    integer,  allocatable :: indi(:)
    integer,  allocatable :: mi(:)
  contains
    procedure, public :: init          => init_factor_normal
    procedure, public :: update_acti   => update_factor_block_acti
    procedure, public :: update_inacti => update_factor_block_inacti
  end type factor_normal_block


contains


  !-----------------------------------------------------------------------------

  subroutine init_factor_normal(this, nobs, nmeas, nfac, start)
    implicit none
    class(factor_normal_block) :: this
    integer,  intent(in) :: nobs
    integer,  intent(in) :: nmeas
    integer,  intent(in) :: nfac
    real(r8), intent(in) :: start(nobs,nfac)
    integer              :: i

    call this%factor_normal%init(nobs, nmeas, nfac, start)

    allocate(this%indi(nfac))
    allocate(this%mi(nmeas))
    this%indi = (/(i, i=1, this%nfac)/)
    this%mi = (/(i, i=1, this%nmeas)/)

  end subroutine init_factor_normal


  !-----------------------------------------------------------------------------

  subroutine update_factor_block_acti(this, Y, alpha, dedic, idioprec, fdist)
    implicit none
    class(factor_normal_block)                 :: this
    real(r8),                       intent(in) :: Y(this%nobs,this%nmeas)
    real(r8),                       intent(in) :: idioprec(this%nmeas)
    real(r8),                       intent(in) :: alpha(this%nmeas)
    class(indic_dedic),             intent(in) :: dedic
    class(covmat_block_invwishart), intent(in) :: fdist
    logical                                    :: acti(this%nfac)
    integer                                    :: actind(dedic%K1)
    real(r8)                                   :: mean_post(this%nobs,dedic%K1)
    real(r8)                                   :: var_post(dedic%K1,dedic%K1)
    integer                                    :: i, k

    if(dedic%K1 == 0) return    ! no active factors to sample

    acti   = dedic%active
    actind = pack(this%indi, acti)

    ! compute posterior moments for active factors
    var_post = matinv(fdist%var(actind,actind))
    forall(k = 1:dedic%K1)
      var_post(k,k) = var_post(k,k) &
                    + sum(idioprec*(alpha**2), mask=(dedic%group==actind(k)))
      mean_post(:,k) = matmul(Y(:, pack(this%mi, dedic%group==actind(k))), &
                              pack(idioprec*alpha, dedic%group==actind(k)))
    end forall
    var_post = matinv(var_post)

    ! sample active factors factors
    do k = 1, this%nfac
      if(.not.acti(k)) cycle
      do i = 1, this%nobs
        this%theta(i,k) = rnorm()
      end do
    end do
    this%theta(:,actind) = matmul(mean_post, var_post) &
                         + matmul(this%theta(:,actind), &
                                  transpose(chol(var_post)))

  end subroutine update_factor_block_acti


  !-----------------------------------------------------------------------------

  subroutine update_factor_block_inacti(this, fdist, dedic)
    implicit none
    class(factor_normal_block)                 :: this
    class(covmat_block_invwishart), intent(in) :: fdist
    class(indic_dedic),             intent(in) :: dedic
    integer                                    :: actind(dedic%K1)
    integer                                    :: inactind(dedic%K2)
    integer                                    :: i, k

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

  end subroutine update_factor_block_inacti


end module factor_normal_block_class

