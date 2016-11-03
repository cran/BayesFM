module covmat_block_invwishart_class

  use global
  use indicators_dedic_class
  use matrix,      only : crossprod, matinv, chol
  use probability, only : rnorm, rgamma, rinvwishart
  implicit none

  private
  public :: covmat_block_invwishart
  public :: covmat_block_HuangWand


  type :: covmat_block_invwishart
    integer               :: nfac
    integer               :: npar
    real(r8), allocatable :: prec(:,:)
    real(r8), allocatable :: var(:,:)
    logical,  allocatable :: var_mask(:,:)  ! mask for lower triangular elements
    real(r8)              :: nu0            ! prior in expanded model:
    real(r8), allocatable :: S0(:)          !   var ~ IWish(nu0, diag(S0))
    real(r8)              :: df_post
    real(r8), allocatable :: L221(:,:)      ! ingredients used to sample
    real(r8), allocatable :: G1112(:,:)     ! inactive factors later
    integer,  allocatable :: ind(:)
    ! back up
    real(r8), allocatable :: var_bak(:,:)
    real(r8), allocatable :: prec_bak(:,:)
  contains
    procedure, public :: init    => init_covmat_block_invwishart
    procedure, public :: update  => update_covmat_block_invwishart
    procedure, public :: get     => get_covmat_block_invwishart
    procedure, public :: backup  => backup_covmat_block_invwishart
    procedure, public :: restore => restore_covmat_block_invwishart
  end type covmat_block_invwishart


  type, extends(covmat_block_invwishart) :: covmat_block_HuangWand
    real(r8), allocatable :: A2k(:)
    real(r8)              :: nus
    real(r8)              :: Rmat_b_shape
  end type covmat_block_HuangWand


contains


  !-----------------------------------------------------------------------------

  subroutine init_covmat_block_invwishart(this, nobs, nfac, prior, start)
    implicit none
    class(covmat_block_invwishart) :: this
    integer,  intent(in)           :: nobs
    integer,  intent(in)           :: nfac
    real(r8), intent(in)           :: prior(1:nfac+1)
    real(r8), intent(in)           :: start(nfac, nfac)
    integer                        :: i, j

    allocate(this%prec(nfac, nfac))
    allocate(this%var(nfac, nfac))
    allocate(this%var_mask(nfac, nfac))
    allocate(this%S0(nfac))
    allocate(this%var_bak(nfac, nfac))
    allocate(this%prec_bak(nfac, nfac))
    this%nfac = nfac
    this%npar = nfac*(nfac-1)/2

    ! mask for lower triangular elements
    this%var_mask = .false.
    forall(i = 1:nfac, j = 1:nfac, j < i)
      this%var_mask(i,j) = .true.
    end forall

    ! starting values
    this%var = start
    this%prec = matinv(this%var)
    this%var_bak = start
    this%prec_bak = this%prec

    ! prior
    this%nu0 = prior(1)
    this%S0  = prior(2:)
    this%df_post = this%nu0 + dble(nobs)

    ! additional parameters for Huang-wand prior
    select type (this)
      type is (covmat_block_HuangWand)
        allocate(this%A2k(nfac))
        this%A2k = prior(2:)
        this%nus = this%nu0 - this%nfac + 1._r8
        this%Rmat_b_shape = .5_r8*(this%nu0 + 1._r8)
    end select

    ! factor indicators 1, ..., nfac
    allocate(this%ind(nfac))
    this%ind = (/ (i, i=1, nfac) /)

  end subroutine init_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  subroutine update_covmat_block_invwishart(this, fac, dedic)
    implicit none
    class(covmat_block_invwishart)        :: this
    real(r8), dimension(:,:), intent(in)  :: fac
    class(indic_dedic),       intent(in)  :: dedic
    integer,  dimension(:),   allocatable :: acti, inacti
    real(r8), dimension(:,:), allocatable :: S_post, G221, G12
    real(r8)                              :: df
    real(r8)                              :: Rmat_b_rate
    integer                               :: K1, K2, k, l

    K1 = dedic%K1       ! active factors
    K2 = dedic%K2       ! inactive factors
    allocate(acti(K1))
    acti = pack(this%ind, dedic%active)

    !----- update scale parameters for Huang-Wand prior

    select type (this)
      type is (covmat_block_HuangWand)
        do k = 1, this%nfac
          Rmat_b_rate = .5_r8*(this%prec(k,k) + 1._r8/(this%nus*this%A2k(k)))
          this%S0(k) = rgamma(this%Rmat_b_shape, 1._r8/Rmat_b_rate)
        end do
    end select

    !----- sample submatrix for active factors

    if(K1 > 0) then

      allocate(S_post(K1, K1))

      S_post = crossprod(fac(:,acti))
      forall(k = 1:K1) S_post(k,k) = S_post(k,k) + this%S0(acti(k))
      df = this%df_post - dble(K2)
      this%var(acti,acti) = rinvwishart(df, S_post)

      deallocate(S_post)

    end if

    !----- sample remaining block corresponding to inactive factors

    if(allocated(this%L221)) deallocate(this%L221)

    if(K1 == 0) then      ! only inactive factors, sample matrix from prior

      allocate(this%L221(K2,K2), S_post(K2,K2))
      S_post = 0._r8
      forall(k = 1:K2) S_post(k,k) = this%S0(k)
      this%var = rinvwishart(this%nu0, S_post)
      this%L221 = chol(this%var)

    else if(K2 > 0) then

      if(allocated(this%G1112)) deallocate(this%G1112)
      allocate(this%L221(K2,K2), this%G1112(K1,K2))
      allocate(inacti(K2), S_post(K2,K2), G221(K2,K2), G12(K1,K2))

      inacti = pack(this%ind, mask=.not.dedic%active)

      S_post = 0._r8
      forall(k = 1:K2) S_post(k,k) = this%S0(inacti(k))
      G221 = rinvwishart(this%nu0, S_post)    ! G221 = Omega_22.1
      this%L221  = chol(G221)
      do k = 1, K1                            ! G1112 = inv(Omega_11) * Omega_12
        do l = 1, K2
          this%G1112(k,l) = rnorm()
        end do
      end do
      this%G1112 = matmul(this%G1112, transpose(this%L221))
      forall(k = 1:K1) this%G1112(k,:) = this%G1112(k,:)/sqrt(this%S0(acti(k)))
      G12 = matmul(this%var(acti,acti), this%G1112)      ! G12 = Omega_12

      this%var(acti, inacti) = G12
      this%var(inacti, acti) = transpose(G12)
      this%var(inacti, inacti) = G221 + matmul(transpose(G12), this%G1112)

    end if

    this%prec = matinv(this%var)

  end subroutine update_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  function get_covmat_block_invwishart(this) result(par)
    implicit none
    class(covmat_block_invwishart) :: this
    real(r8)                       :: par(this%npar)

    par = pack(this%var, this%var_mask)

  end function get_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  subroutine backup_covmat_block_invwishart(this)
    implicit none
    class(covmat_block_invwishart) :: this

    this%var_bak  = this%var
    this%prec_bak = this%prec

  end subroutine backup_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  subroutine restore_covmat_block_invwishart(this)
    implicit none
    class(covmat_block_invwishart) :: this

    this%var  = this%var_bak
    this%prec = this%prec_bak

  end subroutine restore_covmat_block_invwishart


end module covmat_block_invwishart_class

