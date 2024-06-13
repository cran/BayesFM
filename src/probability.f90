module probability

  use global
  use matrix, only : chol, matinv
  implicit none

  private

  public :: set_seed
  public :: runif
  public :: rexpon
  public :: rnorm
  public :: rtnorm
  public :: rgamma
  public :: rinvgamma
  public :: rdirich
  public :: rwishart
  public :: rinvwishart
  public :: rpoisson

  real(r8),    parameter :: pi    = 3.141592653589793238462643383276_r8
  real(r8),    parameter :: half  = .5_r8
  real(r8),    parameter :: third = 1._r8/3._r8


  !----- variables for 64-bit Mersenne Twister algorithm
  integer(i8), parameter :: nn      = 312_i8
  integer(i8), parameter :: mm      = 156_i8
  integer(i8), parameter :: matx_a  = -5403634167711393303_i8
  integer(i8), parameter :: um      = -2147483648_i8 ! most significant 33 bits
  integer(i8), parameter :: lm      =  2147483647_i8 ! least significant 31 bits
  real(r8),    parameter :: pi253_1 = 1._r8/(2._r8**53 - 1._r8)

  integer(i8) :: mt(nn)       ! array for the state vector
  integer     :: mti = nn+1   ! mti==nn+1 means mt(nn) is not initialized


  interface runif
   module procedure runif_01
   module procedure runif_ab
  end interface runif

 interface rnorm
   module procedure rnorm_01
   module procedure rnorm_mu_var
  end interface rnorm


contains


  !-----------------------------------------------------------------------------
  ! 64-bit Mersenne Twister algorithm
  !
  ! Fortran translation from C-program for MT19937-64 (2004/9/29 version)
  ! originally coded by Takuji Nishimura and Makoto Matsumoto
  ! see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
  ! translation by Rémi Piatek, University of Copenhagen

  !----- initializes mt(nn) with a seed

  subroutine set_seed(seed)
    implicit none
    integer, intent(in) :: seed
    integer :: i

    mt(1) = int(seed, kind=i8)
    do i = 1, nn-1
      mt(i+1) = 6364136223846793005_i8 * ieor(mt(i), ishft(mt(i), -62)) + i
    end do
    mti = nn

  end subroutine set_seed


  !----- generates a random number on [-2^63, 2^63-1]-interval

  integer(r8) function genrand64_int64()
    implicit none
    integer(i8) :: mag01(0:1) = (/0_i8, matx_a/)
    integer(i8) :: x
    integer     :: i

    if(mti >= nn) then ! generate nn words at one time

      ! if set_seed() has not been called, a default initial seed is used
      if(mti == nn+1) call set_seed(5489)

      do i = 1, nn-mm
        x = ior(iand(mt(i),um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do

      do i = nn-mm+1, nn-1
        x = ior(iand(mt(i), um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm-nn), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do

      x = ior(iand(mt(nn), um), iand(mt(1), lm))
      mt(nn) = ieor(ieor(mt(mm), ishft(x, -1)), mag01(iand(x, 1_i8)))

      mti = 0

    end if

    mti = mti + 1
    x = mt(mti)

    x = ieor(x, iand(ishft(x,-29), 6148914691236517205_i8))
    x = ieor(x, iand(ishft(x, 17), 8202884508482404352_i8))
    x = ieor(x, iand(ishft(x, 37),   -2270628950310912_i8))
    x = ieor(x, ishft(x, -43))

    genrand64_int64 = x

  end function genrand64_int64


  !-----------------------------------------------------------------------------
  ! generates random number from uniform distribution (genrand64_real1)

  real(r8) function runif_01()
    implicit none

    runif_01 = real(ishft(genrand64_int64(), -11), kind=r8) * pi253_1

  end function runif_01


  real(r8) function runif_ab(a, b)
    implicit none
    real(r8), intent(in) :: a, b

    if(b <= a) call rexit('*** ERROR: a should be < b (runif) ***')

    runif_ab = runif_01() * (b-a) + a

  end function runif_ab


  !-----------------------------------------------------------------------------
  ! generates random number from exponential distribution
  ! with scale parameter b such that mean = b and variance = b^2

  real(r8) function rexpon(b)
    implicit none
    real(r8), intent(in) :: b

    if(b <= 0._r8) then
      call rexit('*** ERROR: rate parameter should be > 0 (rexpon) ***')
    end if

    rexpon = -log(runif()) / b

  end function rexpon


  !-----------------------------------------------------------------------------
  ! generates random number from normal distribution
  !
  ! source:
  !   Joseph L. Leva
  !   ``A Fast Normal Random Number Generator''
  !   ACM Transactions on Mathematical Software
  !   Vol. 18, No. 4, December 1992, pages 449-453

  real(r8) function rnorm_01()
    implicit none
    real(r8) :: u, v, x, y, q

    do
      u = runif()
      v = 1.7156_r8*(runif() - 0.5_r8)
      x = u - 0.449871_r8
      y = abs(v) + 0.386595_r8
      q = x**2 + y*(0.19600_r8*y - 0.25472_r8*x)
      if(q < 0.27597_r8) exit
      if(q > 0.27846_r8) cycle
      if(v**2 <= -4._r8*(u**2)*log(u)) exit    ! rarely evaluated
    end do
    rnorm_01 = v/u

  end function rnorm_01


  real(r8) function rnorm_mu_var(mu, var)
    implicit none
    real(r8), intent(in) :: mu, var

    if(var <= 0._r8) then
      call rexit('*** ERROR: var should be positive (rnorm) ***')
    end if

    rnorm_mu_var = mu + sqrt(var)*rnorm_01()

  end function rnorm_mu_var


  !-----------------------------------------------------------------------------
  ! generates random number from truncated normal distribution
  ! with mean mu and standard deviation sd with truncation
  !   (a, +infty)    if left = true
  !   (-infty, a)    if left = false
  !
  ! source:
  !   John Geweke (1991)
  !   ``Efficient Simulation from the Multivariate Normaland Student-t
  !     Distributions Subject to Linear Constraintsand the Evaluation
  !     of Constraint Probabilities''
  !   Computing Science and Statistics:
  !   Proceedings of the 23rd Symposium on the Interface
  !   Ed. E. Keramidas and S. Kaufman, pages 571-578
  !   Fairfax Station, VA: InterfaceFoundation of North America.

  real(r8) function rtnorm(mu, var, a, left)
    implicit none
    real(r8), intent(in) :: mu, var, a
    logical,  intent(in) :: left
    real(r8) :: sd, c, u, z

    if(var <= 0._r8) then
      call rexit('*** ERROR: var should be positive (rtnorm) ***')
    end if
    sd = sqrt(var)

    c = (a - mu)/sd
    if(.not.left) c = -c

    if(c <= .45_r8) then   ! normal rejection sampling
      do
        z = rnorm()
        if(z > c) exit
      end do
    else                   ! exponential rejection sampling
      do
        z = rexpon(c)
        u = runif()
        if(u < exp(-.5_r8*(z**2))) exit
      end do
      z = z + c
    end if

    if(left) then
       rtnorm = mu + z*sd
    else
       rtnorm = mu - z*sd
    end if

  end function rtnorm


  !-----------------------------------------------------------------------------
  ! generates random number from gamma distribution
  ! with shape parameter a and scale parameter b
  ! such that mean=a*b and var=a*b^2
  !
  ! source:
  !   George Marsaglia and Wai Wan Tsang
  !   ``A Simple Method for Generating Gamma Variables''
  !   ACM Transactions on Mathematical Software
  !   Vol. 26, No. 3, September 2000, pages 363-372

  real(r8) function rgamma(a, b)
    implicit none
    real(r8), intent(in) :: a, b
    real(r8) :: a1, c, d, u, v, x

    if(a <= 0._r8) call rexit('*** ERROR: a should be positive (rgamma) ***')
    if(b <= 0._r8) call rexit('*** ERROR: b should be positive (rgamma) ***')

    a1 = a
    if(a < 1._r8) a1 = a + 1._r8
    d = a1 - third
    c = 1._r8/sqrt(9._r8*d)

    do
      do
        x = rnorm()
        v = (1._r8 + c*x)
        if(v > 0._r8) exit
      end do
      v = v**3
      u = runif()
      if(u < 1._r8 - 0.0331_r8*(x**4)) exit
      if(log(u) < .5_r8*(x**2) + d*(1._r8-v+log(v))) exit
    end do
    rgamma = d*v*b

    ! case where a < 1 (see note p.371)
    if(a < 1._r8) then
      do
        u = runif()
        if(u > 0._r8) then   ! cycle if u = 0
          rgamma = rgamma * u**(1._r8/a)
          exit
        end if
      end do
    end if

  end function rgamma


  !-----------------------------------------------------------------------------
  ! generates random number from inverse-gamma distribution
  ! with shape parameter a and rate parameter b

  real(r8) function rinvgamma(a, b)
    implicit none
    real(r8), intent(in) :: a, b

    if(a <= 0._r8) call rexit('*** ERROR: a should be positive (rinvgamma) ***')
    if(b <= 0._r8) call rexit('*** ERROR: b should be positive (rinvgamma) ***')

    rinvgamma = 1._r8/rgamma(a, 1._r8/b)

  end function rinvgamma


  !-----------------------------------------------------------------------------
  ! generates vector from Dirichlet distribution
  ! with concentration parameters alpha

  function rdirich(alpha)
    implicit none
    real(r8), intent(in) :: alpha(:)
    real(r8)             :: rdirich(size(alpha))
    integer              :: i

    if(any(alpha <= 0._r8)) then
      call rexit('*** ERROR: alpha should be strictly positive (rdirich) ***')
    end if

    do i = 1, size(alpha)
      rdirich(i) = rgamma(alpha(i), 1._r8)
    end do
    rdirich = rdirich / sum(rdirich)

  end function rdirich


  !-----------------------------------------------------------------------------
  ! generates random matrix from Wishart distribution
  ! with df degrees of freedom and scale matrix S:
  !
  !   X ~ Wishart(df, S)
  !   p(X) \propto |X|^{(df-p-1)/2} exp{ -tr(X * S^-1)/2 }
  !
  ! reference:
  !   A. K. Gupta & D. K. Nagar (2000)
  !   ``Matrix Variate Distributions''
  !   Chapman & Hall/CRC
  !   Monographs and Surveys in Pure and Applied Mathematics, No. 104
  !   Theorems 3.3.1 and 3.3.4, pp.90-91

  function rwishart(df, S)
    implicit none
    real(r8), intent(in) :: df, S(:,:)
    real(r8), dimension(size(S,1),size(S,2)) :: rwishart, B, BB, L
    integer :: i, j, p

    p = size(S,1)
    if(df < dble(p)) then
      call rexit('*** ERROR: degrees of freedom should be > p-1 (rwishart) ***')
    end if

    B = 0._r8
    do i = 1, p
      do j = 1, i
        if (i == j) then
          B(i,i) = sqrt(rgamma(.5_r8*(df-dble(i)+1._r8), 2._r8))
        else
          B(i,j) = rnorm()
        end if
      end do
    end do

    BB = matmul(B, transpose(B))
    L = chol(S)
    rwishart = matmul(L, matmul(BB, transpose(L)))

  end function rwishart


  !-----------------------------------------------------------------------------
  ! generates random matrix from inverse-Wishart distribution
  !
  !   X ~ inv-Wishart(df, W)
  !   p(X) \propto |X|^{-(df+p+1)/2} exp{ -tr(W * X^-1)/2 }

  function rinvwishart(df, W)
    implicit none
    real(r8), intent(in) :: df, W(:,:)
    real(r8) :: rinvwishart(size(W,1), size(W,1))

    rinvwishart = matinv(rwishart(df, matinv(W)))

  end function rinvwishart


  !-----------------------------------------------------------------------------
  ! generates random number from Poisson distribution

  integer function rpoisson(b)
    implicit none
    real(r8), intent(in) :: b
    real(r8)             :: em, t
    real(r8), save       :: g, oldb=-1._r8

    if(b <= 0._r8) then
      call rexit('*** ERROR: b should be > 0 (rpoisson) ***')
    end if

    ! Knuth's algorithm
    if(abs(b - oldb) > 0._r8) then
      oldb = b
      g = exp(-b)
    end if
    em = -1._r8
    t = 1._r8
    do
      em = em + 1._r8
      t = t * runif()
      if(t <= g) exit
    end do

    rpoisson = int(em)

  end function rpoisson


end module probability

