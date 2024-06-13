module mda_class

  use global
  use covmat_block_invwishart_class
  use probability, only : rgamma, rinvgamma
  implicit none


  type :: workpar
    integer               :: nfac
    integer               :: nmeas
    integer               :: nobs
    real(r8), allocatable :: D(:)
  end type workpar


contains


  !-----------------------------------------------------------------------------

  subroutine init_workpar(this, nfac, nmeas, nobs)
    implicit none
    type(workpar), intent(out) :: this
    integer,       intent(in)  :: nfac
    integer,       intent(in)  :: nmeas
    integer,       intent(in)  :: nobs

    allocate(this%D(nfac))

    this%nfac  = nfac
    this%nmeas = nmeas
    this%nobs  = nobs

  end subroutine init_workpar


  !-----------------------------------------------------------------------------

  subroutine expand_workpar(this, dedic, alpha, covmat)
    implicit none
    type(workpar),                 intent(inout) :: this
    integer,                       intent(in)    :: dedic(this%nmeas)
    real(r8),                      intent(inout) :: alpha(this%nmeas)
    type(covmat_block_invwishart), intent(inout) :: covmat
    integer                                      :: k, l

    if(covmat%use_HuangWand) then
      ! sample scale parameters from prior for Huang-Wand prior
      do k = 1, this%nfac
        covmat%S0(k) = rgamma(.5_r8, 2._r8*covmat%nus*covmat%A2k(k))
      end do
    end if

    ! sample working parameters from conditional prior
    do k = 1, this%nfac
      this%D(k) = rinvgamma(.5_r8*covmat%nu0, &
                            .5_r8*covmat%S0(k)*covmat%prec(k,k))
    end do
    this%D = sqrt(this%D)

    ! transform model parameters
    ! (note: factors not tansformed, as they are updated right in the next step)
    do k = 1, this%nfac
      where(dedic == k) alpha = alpha / this%D(k)
    end do
    do l = 1, this%nfac  ! k <= l
      do k = 1, l
        covmat%var(k,l)  = covmat%var(k,l)  * this%D(k) * this%D(l)
        covmat%prec(k,l) = covmat%prec(k,l) / this%D(k) / this%D(l)
        covmat%var(l,k)  = covmat%var(k,l)
        covmat%prec(l,k) = covmat%prec(k,l)
      end do
    end do

  end subroutine expand_workpar



  !-----------------------------------------------------------------------------

  subroutine transform_back_workpar(this, dedic, alpha, covmat, fac)
    implicit none
    type(workpar),                  intent(inout) :: this
    integer,                        intent(in)    :: dedic(this%nmeas)
    real(r8),                       intent(inout) :: alpha(this%nmeas)
    type(covmat_block_invwishart),  intent(inout) :: covmat
    real(r8),                       intent(inout) :: fac(this%nobs,this%nfac)
    integer                                       :: k, l

    ! retrieve updated working parameters
    do k = 1, this%nfac
      this%D(k) = sqrt(covmat%var(k,k))
    end do

    ! transform back
    do k = 1, this%nfac
      fac(:,k) = fac(:,k) / this%D(k)
      where(dedic == k) alpha = alpha * this%D(k)
    end do
    do l = 1, this%nfac  ! k <= l
      do k = 1, l
        covmat%var(k,l)  = covmat%var(k,l)  / this%D(k) / this%D(l)
        covmat%prec(k,l) = covmat%prec(k,l) * this%D(k) * this%D(l)
        covmat%var(l,k)  = covmat%var(k,l)
        covmat%prec(l,k) = covmat%prec(k,l)
      end do
    end do

  end subroutine transform_back_workpar


end module mda_class

