module covmat_block_invwishart_class

  use global
  use indicators_dedic_class
  use matrix,      only : crossprod, matinv, chol
  use probability, only : rnorm, rgamma, rinvwishart
  implicit none


  type :: covmat_block_invwishart
    logical               :: use_HuangWand
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
    ! Huang-Wand
    real(r8), allocatable :: A2k(:)
    real(r8)              :: nus
    real(r8)              :: Rmat_b_shape
    ! back up
    real(r8), allocatable :: var_bak(:,:)
    real(r8), allocatable :: prec_bak(:,:)
  end type covmat_block_invwishart


contains


  !-----------------------------------------------------------------------------

  subroutine init_covmat_block_invwishart(this, nobs, nfac, use_HuangWand, &
                                          prior, start)
    implicit none
    type(covmat_block_invwishart), intent(out) :: this
    integer,                       intent(in)  :: nobs
    integer,                       intent(in)  :: nfac
    logical,                       intent(in)  :: use_HuangWand
    real(r8),                      intent(in)  :: prior(1:nfac+1)
    real(r8),                      intent(in)  :: start(nfac, nfac)
    integer                                    :: i, j

    this%use_HuangWand = use_HuangWand

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
    do i = 1, nfac
      do j = 1, nfac
        if(j < i) then
          this%var_mask(i,j) = .true.
        end if
      end do
    end do

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
    if(this%use_HuangWand) then
      allocate(this%A2k(nfac))
      this%A2k = prior(2:)
      this%nus = this%nu0 - this%nfac + 1._r8
      this%Rmat_b_shape = .5_r8*(this%nu0 + 1._r8)
    end if

    ! factor indicators 1, ..., nfac
    allocate(this%ind(nfac))
    this%ind = (/ (i, i=1, nfac) /)

  end subroutine init_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  subroutine update_covmat_block_invwishart(this, fac, dedic)
    implicit none
    type(covmat_block_invwishart), intent(inout) :: this
    real(r8), dimension(:,:),      intent(in)    :: fac
    type(indic_dedic),             intent(in)    :: dedic
    integer,  dimension(:),        allocatable   :: acti, inacti
    real(r8), dimension(:,:),      allocatable   :: S_post, G221, G12
    real(r8)                                     :: df
    real(r8)                                     :: Rmat_b_rate
    integer                                      :: K1, K2, k, l

    K1 = dedic%K1       ! active factors
    K2 = dedic%K2       ! inactive factors

    !----- update scale parameters for Huang-Wand prior

    if(this%use_HuangWand) then
      do k = 1, this%nfac
        Rmat_b_rate = .5_r8*(this%prec(k,k) + 1._r8/(this%nus*this%A2k(k)))
        this%S0(k) = rgamma(this%Rmat_b_shape, 1._r8/Rmat_b_rate)
      end do
    end if

    !----- sample submatrix for active factors

    if(K1 > 0) then

      if(allocated(acti)) deallocate(acti)
      if(allocated(S_post)) deallocate(S_post)
      allocate(acti(K1))
      allocate(S_post(K1, K1))

      acti = pack(this%ind, dedic%active)
      S_post = crossprod(fac(:,acti))
      do k = 1, K1
        S_post(k,k) = S_post(k,k) + this%S0(acti(k))
      end do
      df = this%df_post - dble(K2)
      this%var(acti,acti) = rinvwishart(df, S_post)

    end if

    !----- sample remaining block corresponding to inactive factors


    if(K1 == 0) then      ! only inactive factors, sample matrix from prior

      if(allocated(this%L221)) deallocate(this%L221)
      if(allocated(S_post)) deallocate(S_post)
      allocate(this%L221(K2,K2))
      allocate(S_post(K2,K2))

      S_post = 0._r8
      do k = 1, K2
        S_post(k,k) = this%S0(k)
      end do
      this%var = rinvwishart(this%nu0, S_post)
      this%L221 = chol(this%var)

    else if(K2 > 0) then

      if(allocated(this%L221))  deallocate(this%L221)
      if(allocated(this%G1112)) deallocate(this%G1112)
      if(allocated(inacti))     deallocate(inacti)
      if(allocated(S_post))     deallocate(S_post)
      if(allocated(G221))       deallocate(G221)
      if(allocated(G12))        deallocate(G12)
      allocate(this%L221(K2,K2))
      allocate(this%G1112(K1,K2))
      allocate(inacti(K2))
      allocate(S_post(K2,K2))
      allocate(G221(K2,K2))
      allocate(G12(K1,K2))

      inacti = pack(this%ind, mask=.not.dedic%active)

      S_post = 0._r8
      do k = 1, K2
        S_post(k,k) = this%S0(inacti(k))
      end do
      G221 = rinvwishart(this%nu0, S_post)    ! G221 = Omega_22.1
      this%L221  = chol(G221)
      do k = 1, K1                            ! G1112 = inv(Omega_11) * Omega_12
        do l = 1, K2
          this%G1112(k,l) = rnorm()
        end do
      end do
      this%G1112 = matmul(this%G1112, transpose(this%L221))
      do k = 1, K1
        this%G1112(k,:) = this%G1112(k,:)/sqrt(this%S0(acti(k)))
      end do
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
    type(covmat_block_invwishart) :: this
    real(r8)                      :: par(this%npar)

    par = pack(this%var, this%var_mask)

  end function get_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  subroutine backup_covmat_block_invwishart(this)
    implicit none
    type(covmat_block_invwishart), intent(inout) :: this

    this%var_bak  = this%var
    this%prec_bak = this%prec

  end subroutine backup_covmat_block_invwishart


  !-----------------------------------------------------------------------------

  subroutine restore_covmat_block_invwishart(this)
    implicit none
    type(covmat_block_invwishart), intent(inout) :: this

    this%var  = this%var_bak
    this%prec = this%prec_bak

  end subroutine restore_covmat_block_invwishart


end module covmat_block_invwishart_class

