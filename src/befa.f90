subroutine befa &
    (nmeas, nobs, kmax, nid, Yobs, Ycat, Ymiss, nX, Xobs, Xloc, iter, burnin, &
     search_delay, Rmat_delay, step_rnd, step_lambda, seed, prior_loadprec, &
     prior_beta, prior_dedic, prior_facdist, start_loadprec, start_beta, &
     start_dedic, start_factor, start_facdist, verbose, npar, MCMCdraws, &
     MCMCdedic, MHacc)

  use global
  use probability, only : set_seed, rpoisson
  use measurement_class
  use loading_idioprec_class
  use covariates_class
  use indicators_dedic_class
  use factor_normal_block_class
  use covmat_block_invwishart_class
  use mda_class
  use mcmc_progress_class
  implicit none

  !----- input
  integer,  intent(in) :: nmeas
  integer,  intent(in) :: nobs
  integer,  intent(in) :: kmax
  integer,  intent(in) :: nid
  real(r8), intent(in) :: Yobs(nobs,nmeas)
  integer,  intent(in) :: Ycat(nmeas)
  logical,  intent(in) :: Ymiss(nobs,nmeas)
  integer,  intent(in) :: nX
  real(r8), intent(in) :: Xobs(nobs,nX)
  logical,  intent(in) :: Xloc(nmeas,nX)
  integer,  intent(in) :: iter
  integer,  intent(in) :: burnin
  integer,  intent(in) :: search_delay
  integer,  intent(in) :: Rmat_delay
  logical,  intent(in) :: step_rnd
  real(r8), intent(in) :: step_lambda
  integer,  intent(in) :: seed
  real(r8), intent(in) :: prior_loadprec(nmeas,3)
  real(r8), intent(in) :: prior_beta(nmeas)
  real(r8), intent(in) :: prior_dedic(0:3*nmeas+kmax+2)
  real(r8), intent(in) :: prior_facdist(0:kmax+1)
  real(r8), intent(in) :: start_loadprec(nmeas,2)
  real(r8), intent(in) :: start_beta(nmeas,nX)
  integer,  intent(in) :: start_dedic(nmeas)
  real(r8), intent(in) :: start_factor(nobs,kmax)
  real(r8), intent(in) :: start_facdist(kmax,kmax)
  logical,  intent(in) :: verbose

  !----- output
  integer,  intent(in)  :: npar
  real(r8), intent(out) :: MCMCdraws(burnin+1:iter,npar)
  integer,  intent(out) :: MCMCdedic(burnin+1:iter,nmeas)
  logical,  intent(out) :: MHacc(burnin+1:iter)

  !----- model ingredients
  type(measurement)             :: Ylat(nmeas)
  type(covariates)              :: Xcov(nmeas)
  type(loading_idioprec)        :: alpha_prec(nmeas)
  type(indic_dedic)             :: dedic
  type(factor_normal_block)     :: factors
  type(covmat_block_invwishart) :: facdist
  type(workpar)                 :: mda

  logical :: norestr
  integer :: i, ii, j, k, step

  type(mcmc_progress) :: mcmc_prog
  call init_mcmc_progress(mcmc_prog, burnin, iter, verbose)


  !=============================================================================
  ! initialization

  call set_seed(seed)

  do j = 1, nmeas

    ! latent variables underlying measurements
    call init_measurement(Ylat(j), nobs, Ycat(j) > 0, Yobs(:,j), Ymiss(:,j))

    ! slope parameters for covariates
    call init_covariates(Xcov(j), nobs, nX, Xobs, prior_beta(j), &
                         start_beta(j,:), Xloc(j,:))

    ! factor loadings and idiosyncratic precisions
    call init_loading_idioprec(alpha_prec(j), nobs, Ycat(j) > 0, &
                               prior_loadprec(j,:), start_loadprec(j,:))

  end do

  ! indicators
  call init_indic_dedic(dedic, nobs, nmeas, kmax, Ycat, prior_dedic, &
                        start_dedic)
  where(dedic%group == 0) alpha_prec%alpha = 0._r8

  ! latent factors
  call init_factor_normal_block(factors, nobs, nmeas, kmax, start_factor)
  if(int(prior_facdist(0)) == 0) then
    call init_covmat_block_invwishart(facdist, nobs, kmax, .false., &
                                      prior_facdist(1:), start_facdist)
  else
    call init_covmat_block_invwishart(facdist, nobs, kmax, .true., &
                                      prior_facdist(1:), start_facdist)
  end if

  ! working parameters for MDA
  call init_workpar(mda, kmax, nmeas, nobs)

  norestr = nid==1
  step = int(step_lambda)
  MHacc = .false.


  !=============================================================================
  ! MCMC sampling

  do i = 1, iter

    if(i <= search_delay) then

      call gibbs_sweep_forward(search=.false., sample_Rmat = i>Rmat_delay)
      if(i == search_delay) call backup_draws()  ! backup before starting search

    else if(norestr) then

      call gibbs_sweep_forward(search=.true., sample_Rmat = i>Rmat_delay)

    else

      ! number of intermediate steps in expanded model
      if(step_rnd) step = 1 + rpoisson(step_lambda)

      ! intermediate steps forward
      do ii = 1, step
        call gibbs_sweep_forward(search=.true., sample_Rmat = i>Rmat_delay)
      end do

      ! intermediate steps backward
      do ii = 1, step
        call gibbs_sweep_backward(sample_Rmat = i>Rmat_delay)
      end do

      ! M-H step: check identification restrictions
      if(any(dedic%ngroup<nid .and. dedic%ngroup>0)) then
        ! model is not identified, restore previous draws
        call restore_draws()
      else
        ! model is identified, back up current draws for next iteration
        call backup_draws()
        if(i > burnin) MHacc(i) = .true.
      end if

    end if

    ! save current draws after burn-in
    if(i > burnin) then
      MCMCdedic(i,:) = dedic%group
      MCMCdraws(i,:) = [ alpha_prec%alpha, &
                         alpha_prec%var, &
                         get_covmat_block_invwishart(facdist), &
                         get_all_covariates(Xcov) ]
    end if

    ! show MCMC progress
    call show_mcmc_progress(mcmc_prog, i)

  end do


contains


  !=============================================================================

  subroutine gibbs_sweep_forward(search, sample_Rmat)
    implicit none
    logical, intent(in) :: search, sample_Rmat
    real(r8)            :: Yaux(nobs, nmeas)

    do j = 1, nmeas
      Yaux(:,j) = Ylat(j)%Y - Xcov(j)%Xbeta
    end do

    ! sample indicators
    if(search) then
      call update_indic_dedic(dedic, Yaux, factors%theta)
      where(dedic%group == 0) alpha_prec%alpha = 0._r8
    end if

    do j = 1, nmeas

      k = dedic%group(j)

      ! sample factor loading and idiosyncratic precision
      call update_loading_idioprec(alpha_prec(j), Yaux(:,j), k, factors%theta)

      ! sample slope parameters for covariates
      Yaux(:,j) = Ylat(j)%Y
      if(k > 0) Yaux(:,j) = Yaux(:,j) - alpha_prec(j)%alpha*factors%theta(:,k)
      call update_covariates(Xcov(j), Yaux(:,j), alpha_prec(j)%prec)

      ! sample latent variable underlying measurement
      if(Ycat(j)>0 .or. any(Ymiss(:,j))) then
        Yaux(:,j) = Xcov(j)%Xbeta
        if(k > 0) Yaux(:,j) = Yaux(:,j) + alpha_prec(j)%alpha*factors%theta(:,k)
        call update_measurement(Ylat(j), Yaux(:,j), alpha_prec(j)%prec)
      end if

    end do

    if(sample_Rmat) then

      ! MDA: expand model
      call expand_workpar(mda, dedic%group, alpha_prec%alpha, facdist)

      ! sample active factors
      do j = 1, nmeas
        Yaux(:,j) = Ylat(j)%Y - Xcov(j)%Xbeta
      end do
      call update_factor_normal_block_acti &
             (factors, Yaux, alpha_prec%alpha, dedic, alpha_prec%prec, facdist)

      ! sample factor covariance matrix
      call update_covmat_block_invwishart(facdist, factors%theta, dedic)

      ! sample inactive factors
      call update_factor_normal_block_inacti(factors, facdist, dedic)

      ! MDA: transform back to identified model
      call transform_back_workpar(mda, dedic%group, alpha_prec%alpha, facdist, &
                                  factors%theta)

    else

      do j = 1, nmeas
        Yaux(:,j) = Ylat(j)%Y - Xcov(j)%Xbeta
      end do
      call update_factor_normal_block(factors, Yaux, alpha_prec%alpha, &
                                      dedic%group, alpha_prec%prec, facdist)

    end if

  end subroutine gibbs_sweep_forward


  !=============================================================================

  subroutine gibbs_sweep_backward(sample_Rmat)
    implicit none
    logical, intent(in) :: sample_Rmat
    real(r8)            :: Yaux(nobs, nmeas)

    if(sample_Rmat) then

      ! MDA: expand model
      call expand_workpar(mda, dedic%group, alpha_prec%alpha, facdist)

      ! sample active factors
      do j = 1, nmeas
        Yaux(:,j) = Ylat(j)%Y - Xcov(j)%Xbeta
      end do
      call update_factor_normal_block_acti &
             (factors, Yaux, alpha_prec%alpha, dedic, alpha_prec%prec, facdist)

      ! sample factor covariance matrix
      call update_covmat_block_invwishart(facdist, factors%theta, dedic)

      ! sample inactive factors
      call update_factor_normal_block_inacti(factors, facdist, dedic)

      ! MDA: transform back to identified model
      call transform_back_workpar(mda, dedic%group, alpha_prec%alpha, facdist, &
                                  factors%theta)

    else

      do j = 1, nmeas
        Yaux(:,j) = Ylat(j)%Y - Xcov(j)%Xbeta
      end do
      call update_factor_normal_block(factors, Yaux, alpha_prec%alpha, &
                                      dedic%group, alpha_prec%prec, facdist)

    end if

    do j = 1, nmeas

      k = dedic%group(j)

      ! sample latent variable underlying measurement
      if(Ycat(j)>0 .or. any(Ymiss(:,j))) then
        Yaux(:,j) = Xcov(j)%Xbeta
        if(k > 0) Yaux(:,j) = Yaux(:,j) + alpha_prec(j)%alpha*factors%theta(:,k)
        call update_measurement(Ylat(j), Yaux(:,j), alpha_prec(j)%prec)
      end if

      ! sample slope parameters for covariates
      Yaux(:,j) = Ylat(j)%Y
      if(k > 0) Yaux(:,j) = Yaux(:,j) - alpha_prec(j)%alpha*factors%theta(:,k)
      call update_covariates(Xcov(j), Yaux(:,j), alpha_prec(j)%prec)

    end do

    ! sample indicators
    do j = 1, nmeas
      Yaux(:,j) = Ylat(j)%Y - Xcov(j)%Xbeta
    end do
    call update_indic_dedic(dedic, Yaux, factors%theta)
    where(dedic%group == 0) alpha_prec%alpha = 0._r8

    ! sample factor loading and idiosyncratic precision
    do j = 1, nmeas
      call update_loading_idioprec(alpha_prec(j), Yaux(:,j), dedic%group(j), &
                                   factors%theta)
    end do

  end subroutine gibbs_sweep_backward


  !=============================================================================

  subroutine backup_draws()
    implicit none
    integer :: j

    call backup_indic_dedic(dedic)
    call backup_factor_normal_block(factors)
    call backup_covmat_block_invwishart(facdist)
    do j = 1, nmeas
      call backup_loading_idioprec(alpha_prec(j))
      call backup_covariates(Xcov(j))
      call backup_measurement(Ylat(j))
    end do

  end subroutine backup_draws


  !=============================================================================

  subroutine restore_draws()
    implicit none
    integer :: j

    call restore_indic_dedic(dedic)
    call restore_factor_normal_block(factors)
    call restore_covmat_block_invwishart(facdist)
    do j = 1, nmeas
      call restore_loading_idioprec(alpha_prec(j))
      call restore_covariates(Xcov(j))
      call restore_measurement(Ylat(j))
    end do

  end subroutine restore_draws


end subroutine befa

