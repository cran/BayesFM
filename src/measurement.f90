module measurement_class

  use global
  use probability, only : rnorm, rtnorm
  implicit none


  type :: measurement
    logical               :: is_binary  ! TRUE = binary, FALSE = continuous
    integer               :: nobs
    real(r8), allocatable :: Y(:)
    logical,  allocatable :: Ybin(:)
    logical,  allocatable :: Ymiss(:)
    ! back up
    real(r8), allocatable :: Y_bak(:)
  end type measurement


contains


  !-----------------------------------------------------------------------------

  subroutine init_measurement(this, nobs, is_binary, Y, Ymiss)
    implicit none
    type(measurement), intent(out)          :: this
    integer,           intent(in)           :: nobs
    logical,           intent(in)           :: is_binary
    real(r8),          intent(in)           :: Y(nobs)
    logical,           intent(in), optional :: Ymiss(nobs)
    integer                                 :: i

    this%is_binary = is_binary

    this%nobs = nobs
    allocate(this%Y(nobs))

    if(present(Ymiss)) then
      if(any(Ymiss)) then
        allocate(this%Ymiss(nobs))
        this%Ymiss = Ymiss
      end if
    end if

    if(this%is_binary) then
      allocate(this%Ybin(nobs))
      allocate(this%Y_bak(nobs))
      this%Ybin = int(Y) == 1
      do i = 1, nobs
        this%Y(i) = abs(rnorm())
      end do
      if(allocated(this%Ymiss)) then
        do i = 1, nobs
          if(.not.this%Ymiss(i) .and. .not.this%Ybin(i)) then
            this%Y(i) = -this%Y(i)
          end if
        end do
      else
        do i = 1, nobs
          if(.not.this%Ybin(i)) then
            this%Y(i) = -this%Y(i)
          end if
        end do
      end if
    else
      this%Y = Y
      if(allocated(this%Ymiss)) then
        allocate(this%Y_bak(nobs))
        do i = 1, nobs
          if(.not.this%Ymiss(i)) cycle
          this%Y(i) = rnorm()
        end do
      end if
    end if

    if(allocated(this%Y_bak)) this%Y_bak = this%Y

  end subroutine init_measurement


  !-----------------------------------------------------------------------------

  subroutine update_measurement(this, mean, var)
    implicit none
    type(measurement), intent(inout) :: this
    real(r8),          intent(in)    :: mean(this%nobs)
    real(r8),          intent(in)    :: var
    integer                          :: i

    if(this%is_binary) then

      if(allocated(this%Ymiss)) then
  
        do i = 1, this%nobs
          if(this%Ymiss(i)) then
            this%Y(i) = rnorm(mean(i), var)
          else
            this%Y(i) = rtnorm(mean(i), var, 0._r8, this%Ybin(i))
          end if
        end do
  
      else
  
        do i = 1, this%nobs
          this%Y(i) = rtnorm(mean(i), var, 0._r8, this%Ybin(i))
        end do
  
      end if

    else

      if(.not.allocated(this%Ymiss)) return
  
      do i = 1, this%nobs
        if(this%Ymiss(i)) then
          this%Y(i) = rnorm(mean(i), var)
        end if
      end do

    end if

  end subroutine update_measurement


  !-----------------------------------------------------------------------------

  subroutine backup_measurement(this)
    implicit none
    type(measurement), intent(inout) :: this

    if(allocated(this%Y_bak)) then
      this%Y_bak = this%Y
    end if

  end subroutine backup_measurement


  !-----------------------------------------------------------------------------

  subroutine restore_measurement(this)
    implicit none
    type(measurement), intent(inout) :: this

    if(allocated(this%Y_bak)) then
      this%Y = this%Y_bak
    end if

  end subroutine restore_measurement


end module measurement_class

