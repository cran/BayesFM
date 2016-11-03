module measurement_class

  use global
  use probability, only : rnorm, rtnorm
  implicit none


  private
  public :: measurement
  public :: measurement_cont
  public :: measurement_bin


  type :: measurement_cont
    integer               :: nobs
    real(r8), allocatable :: Y(:)
    logical,  allocatable :: Ymiss(:)
    ! back up
    real(r8), allocatable :: Y_bak(:)
  contains
    procedure, public :: init    => init_measurement
    procedure, public :: update  => update_measurement_cont
    procedure, public :: backup  => backup_measurement
    procedure, public :: restore => restore_measurement
  end type measurement_cont


  type, extends(measurement_cont) :: measurement_bin
    logical, allocatable :: Ybin(:)
  contains
    procedure, public :: update => update_measurement_bin
  end type measurement_bin


  type :: measurement
    class(measurement_cont), allocatable :: p
  end type measurement


contains


  !-----------------------------------------------------------------------------

  subroutine init_measurement(this, nobs, Y, Ymiss)
    implicit none
    class(measurement_cont)        :: this
    integer,  intent(in)           :: nobs
    real(r8), intent(in)           :: Y(nobs)
    logical,  intent(in), optional :: Ymiss(nobs)
    integer                        :: i

    this%nobs = nobs
    allocate(this%Y(nobs))

    if(present(Ymiss)) then
      if(any(Ymiss)) then
        allocate(this%Ymiss(nobs))
        this%Ymiss = Ymiss
      end if
    end if

    select type (this)

      type is (measurement_cont)
        this%Y = Y
        if(allocated(this%Ymiss)) then
          allocate(this%Y_bak(nobs))
          do i = 1, nobs
            if(.not.this%Ymiss(i)) cycle
            this%Y(i) = rnorm()
          end do
        end if

      type is (measurement_bin)
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

    end select

    if(allocated(this%Y_bak)) this%Y_bak = this%Y

  end subroutine init_measurement


  !-----------------------------------------------------------------------------

  subroutine update_measurement_cont(this, mean, var)
    implicit none
    class(measurement_cont) :: this
    real(r8), intent(in)    :: mean(this%nobs)
    real(r8), intent(in)    :: var
    integer                 :: i

    if(.not.allocated(this%Ymiss)) return

    do i = 1, this%nobs
      if(this%Ymiss(i)) then
        this%Y(i) = rnorm(mean(i), var)
      end if
    end do

  end subroutine update_measurement_cont


  !-----------------------------------------------------------------------------

  subroutine update_measurement_bin(this, mean, var)
    implicit none
    class(measurement_bin) :: this
    real(r8), intent(in)   :: mean(this%nobs)
    real(r8), intent(in)   :: var
    integer                :: i

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

  end subroutine update_measurement_bin


  !-----------------------------------------------------------------------------

  subroutine backup_measurement(this)
    implicit none
    class(measurement_cont) :: this

    if(.not.allocated(this%Y_bak)) return

    this%Y_bak = this%Y

  end subroutine backup_measurement


  !-----------------------------------------------------------------------------

  subroutine restore_measurement(this)
    implicit none
    class(measurement_cont) :: this

    if(.not.allocated(this%Y_bak)) return

    this%Y = this%Y_bak

  end subroutine restore_measurement


end module measurement_class

