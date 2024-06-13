module indicators_dedic_class

  use global
  use probability, only : runif
  implicit none


  !-----------------------------------------------------------------------------

  type :: ratio_marglik
    logical  :: is_categorical  ! TRUE = categorical, FALSE = continuous
    integer  :: nfac
    real(r8) :: alpha_prec0
    real(r8) :: prec_scale0
    real(r8) :: prec_shape
  end type ratio_marglik


  !-----------------------------------------------------------------------------

  type :: param_tau
    logical               :: use_alt  ! TRUE = use alternative parametrization
    integer               :: nmeas
    integer               :: nfac
    real(r8)              :: xi0
    real(r8)              :: kappa0_xi0
    real(r8), allocatable :: kappa0(:)
    real(r8)              :: skappa0
    real(r8), allocatable :: atau(:,:)
    real(r8), allocatable :: btau(:)
    real(r8), allocatable :: stau(:)
  end type param_tau


  !-----------------------------------------------------------------------------

  type :: indic_dedic
    integer                          :: nmeas
    integer                          :: nobs
    integer                          :: nfac
    integer,             allocatable :: group(:)
    integer,             allocatable :: ngroup(:)
    logical,             allocatable :: active(:)   ! TRUE for active factors
    integer                          :: K1          ! number of active factors
    integer                          :: K2          ! number of inactive factors
    type(ratio_marglik), allocatable :: mlik(:)
    type(param_tau)                  :: ltau
    ! back up
    integer,             allocatable :: group_bak(:)
    integer,             allocatable :: ngroup_bak(:)
  end type indic_dedic


contains


  !-----------------------------------------------------------------------------

  subroutine init_ratio_marglik(this, is_categorical, nobs, nfac, prior)
    implicit none
    type(ratio_marglik), intent(inout) :: this
    logical,             intent(in)    :: is_categorical
    integer,             intent(in)    :: nobs, nfac
    real(r8),            intent(in)    :: prior(3)

    this%is_categorical = is_categorical
    this%nfac = nfac
    this%alpha_prec0 = prior(1)

    if(.not. this%is_categorical) then
      this%prec_scale0 = prior(3)
      this%prec_shape = prior(2) + .5_r8*dble(nobs)
    end if

  end subroutine init_ratio_marglik


  !-----------------------------------------------------------------------------

  function get_ratio_marglik(this, Y2, fac2, Yfac2) result(lmlik)
    implicit none
    type(ratio_marglik), intent(inout) :: this
    real(r8),            intent(in)    :: Y2, fac2(this%nfac), Yfac2(this%nfac)
    real(r8)                           :: Pm(this%nfac), Qm(this%nfac), Rm(this%nfac)
    real(r8)                           :: CNn, lmlik(0:this%nfac, 0:this%nfac)

    if(this%is_categorical) then

      Pm = 1._r8 + fac2/this%alpha_prec0
      Qm = .5_r8*Yfac2/(this%alpha_prec0+fac2)
      Rm = -.5_r8*log(Pm) + Qm

      lmlik(0, 1:) = Rm

    else

      Pm = 1._r8 + fac2/this%alpha_prec0
      Qm = .5_r8*Yfac2/(this%alpha_prec0+fac2)

      CNn = this%prec_scale0 + .5_r8*Y2
      Rm = -.5_r8*log(Pm) - this%prec_shape*log(CNn - Qm)

      lmlik(0, 1:) = Rm + this%prec_shape*log(CNn)

    end if

    lmlik(0,0) = 0._r8
    lmlik(1:, 0) = -lmlik(0, 1:)
    lmlik(1:,1:) = spread(Rm, dim=1, ncopies=this%nfac) &
                 - spread(Rm, dim=2, ncopies=this%nfac)

  end function get_ratio_marglik


  !-----------------------------------------------------------------------------

  subroutine init_param_tau(this, use_alt, nmeas, nfac, prior)
    implicit none
    type(param_tau), intent(out) :: this
    logical,         intent(in)  :: use_alt
    integer,         intent(in)  :: nmeas, nfac
    real(r8),        intent(in)  :: prior(0:nfac+1)
    integer                      :: j, k

    this%use_alt = use_alt
    this%nmeas = nmeas
    this%nfac = nfac

    allocate(this%kappa0(0:this%nfac))
    this%xi0 = prior(0)
    this%kappa0 = prior(1:nfac+1)
    this%skappa0 = sum(this%kappa0(1:))

    if(this%use_alt) then
      this%kappa0_xi0 = log(this%kappa0(0)) - log(this%xi0)
    end if

    ! tables to avoid recomputing log function at each MCMC iteration
    allocate(this%atau(0:nmeas, 0:nfac))
    allocate(this%btau(0:nmeas))
    allocate(this%stau(0:nmeas))
    do j = 0, nmeas
      do k = 0, nfac
        this%atau(j,k) = log(dble(j) + this%kappa0(k))
      end do
    end do
    do j = 0, nmeas
      this%btau(j) = log(dble(j) + this%xi0)
      this%stau(j) = log(dble(j) + this%skappa0)
    end do

  end subroutine init_param_tau


  !-----------------------------------------------------------------------------

  function get_param_tau(this, kold, ngroup) result(ltau)
    implicit none
    type(param_tau), intent(inout) :: this
    integer,         intent(in)    :: kold, ngroup(this%nfac)
    real(r8)                       :: ltau(0:this%nfac, 0:this%nfac)
    real(r8)                       :: lngk(this%nfac)
    integer                        :: ng(this%nfac)
    integer                        :: k, sng

    ng = ngroup
    if(kold > 0) ng(kold) = ng(kold) - 1
    sng = sum(ng)

    do k = 1, this%nfac
      lngk(k) = this%atau(ng(k), k)
    end do
    ltau(0, 0) = 0._r8
    ltau(1:, 1:) = spread(lngk, dim=1, ncopies=this%nfac) &
                 - spread(lngk, dim=2, ncopies=this%nfac)

    if(this%use_alt) then
      ltau(1:, 0) = this%stau(sng) - lngk + this%kappa0_xi0
    else
      ltau(1:, 0) = this%atau(this%nmeas-sng-1, 0) &
                  + this%stau(sng) &
                  - this%btau(sng) &
                  - lngk
    end if

    ltau(0, 1:) = -ltau(1:, 0)

  end function get_param_tau


  !-----------------------------------------------------------------------------

  subroutine init_indic_dedic(this, nobs, nmeas, nfac, Ycat, prior, start)
    implicit none
    type(indic_dedic), intent(out) :: this
    integer,           intent(in)  :: nobs
    integer,           intent(in)  :: nmeas
    integer,           intent(in)  :: nfac
    integer,           intent(in)  :: Ycat(nmeas)
    real(r8),          intent(in)  :: prior(0:3*nmeas+nfac+2)
    integer,           intent(in)  :: start(nmeas)
    integer                        :: j, k
    real(r8)                       :: prior_marglik(nmeas, 3)

    allocate(this%group(nmeas))
    allocate(this%ngroup(nfac))
    allocate(this%active(nfac))
    allocate(this%group_bak(nmeas))
    allocate(this%ngroup_bak(nfac))

    this%nmeas = nmeas
    this%nobs = nobs
    this%nfac = nfac
    this%group = start
    this%group_bak = start

    do k = 1, this%nfac
      this%ngroup(k) = count(this%group == k)
    end do
    this%ngroup_bak = this%ngroup
    this%active = this%ngroup > 0
    this%K1 = count(this%active)
    this%K2 = nfac - this%K1

    prior_marglik = reshape(prior(1:3*nmeas), shape=[nmeas,3])
    allocate(this%mlik(nmeas))
    do j = 1, nmeas
      call init_ratio_marglik(this%mlik(j), Ycat(j) > 0, nobs, nfac, &
                              prior_marglik(j,:))
    end do

    if(int(prior(0)) == 0) then
      call init_param_tau(this%ltau, .false., nmeas, nfac, prior(3*nmeas+1:))
    else
      call init_param_tau(this%ltau, .true., nmeas, nfac, prior(3*nmeas+1:))
    end if

  end subroutine init_indic_dedic


  !-----------------------------------------------------------------------------

  subroutine update_indic_dedic(this, Yaux, Fac)
    implicit none
    type(indic_dedic), intent(inout) :: this
    real(r8),          intent(in)    :: Yaux(this%nobs, this%nmeas)
    real(r8),          intent(in)    :: Fac(this%nobs, this%nfac)
    real(r8)                         :: Fac2(this%nfac)
    real(r8)                         :: FacYXb2(this%nfac, this%nmeas)
    real(r8)                         :: YXb2(this%nmeas)
    real(r8)                         :: prob_sub(0:this%nfac, 0:this%nfac)
    real(r8)                         :: prob(0:this%nfac)
    real(r8)                         :: sprob, e
    integer                          :: j, k, kold

    ! cross-products used to compute posterior log odds
    Fac2 = sum(Fac**2, dim=1)
    YXb2 = sum(Yaux**2, dim=1)
    FacYXb2 = matmul(transpose(Fac), Yaux)**2

    do j = 1, this%nmeas

      ! compute posterior probabilities
      kold = this%group(j)
      prob_sub = get_ratio_marglik(this%mlik(j), YXb2(j), Fac2, FacYXb2(:,j)) &
               + get_param_tau(this%ltau, kold, this%ngroup)
      prob = 1._r8/sum(exp(prob_sub), dim=2)

      ! sample indicator
      e = runif()
      sprob = 0._r8
      do k = 0, this%nfac
        sprob = sprob + prob(k)
        if(e <= sprob) then
          this%group(j) = k
          if(kold > 0) this%ngroup(kold) = this%ngroup(kold) - 1
          if(k > 0)    this%ngroup(k)    = this%ngroup(k)    + 1
          exit
        end if
      end do

    end do

    ! update indicators for active factors
    this%active = this%ngroup > 0
    this%K1 = count(this%active)
    this%K2 = this%nfac - this%K1

  end subroutine update_indic_dedic


  !-----------------------------------------------------------------------------

  subroutine backup_indic_dedic(this)
    implicit none
    type(indic_dedic), intent(inout) :: this

    this%group_bak  = this%group
    this%ngroup_bak = this%ngroup

  end subroutine backup_indic_dedic


  !-----------------------------------------------------------------------------

  subroutine restore_indic_dedic(this)
    implicit none
    type(indic_dedic), intent(inout) :: this

    this%group  = this%group_bak
    this%ngroup = this%ngroup_bak

  end subroutine restore_indic_dedic


end module indicators_dedic_class

