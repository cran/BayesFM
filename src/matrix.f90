module matrix

  use global
  implicit none

  private
  public :: trace
  public :: chol
  public :: det
  public :: matinv
  public :: solvl
  public :: solvu
  public :: crossprod
  public :: outerprod


  interface outerprod
    module procedure outerprod1
    module procedure outerprod2
  end interface outerprod


contains


  !-----------------------------------------------------------------------------
  ! returns the trace of the matrix A

  real(r8) function trace(A) result(t)
    implicit none
    real(r8), intent(in) :: A(:,:)
    integer :: i, n

    n = size(A, 1)
    if(n /= size(A, 2)) then
      call rexit('### ERROR: trace only for square matrices')
    end if

    t = 0._r8
    do i = 1, n
      t = t + A(i,i)
    end do

  end function trace


  !-----------------------------------------------------------------------------
  ! Cholesky decomposition of a positive-definite symmetric matrix A:
  !     A = L * L'
  ! Input:
  !     A        positive-definite symmetric matrix to be decomposed
  !              (only upper triangular part of A is used)
  ! Output:
  !     Achol    Cholesky factor L is returned

  function chol(A) result(Achol)
    implicit none
    real(r8), intent(in) :: A(:,:)
    real(r8) :: Achol(size(A,1), size(A,2))
    real(r8) :: p(size(A,1))
    real(r8) :: asum
    integer  :: i, j, n

    n = size(A,1)
    if(n .ne. size(A,2)) then
      call rexit('*** ERROR: matrix is not square (chol) ***')
    end if

    Achol = A
    do i = 1, n
      asum = Achol(i,i) - sum(Achol(i,1:i-1)**2)
      if(asum <= 0.0_r8) call rexit('*** ERROR: chol failed')
      p(i) = sqrt(asum)
      Achol(i,i) = p(i)
      Achol(i+1:n,i) = (Achol(i,i+1:n) - matmul(Achol(i+1:n,1:i-1), &
                                                Achol(i,1:i-1))) / p(i)
    end do

    do i = 1, n
      do j = 1, n
        if(i < j) then
          Achol(i,j) = 0._r8
        end if
      end do
    end do

  end function chol


  !------------------------------------------------------------------------------
  ! solves the set of linear equations R * x = b
  ! where R is a triangular matrix that is either upper triangular (solvu)
  ! or lower triangular (solvl)

  function solvu(U, b) result(x)
    ! U is upper triangular
    implicit none
    real(r8), intent(in) :: U(:,:)
    real(r8), intent(in) :: b(:)
    real(r8) :: x(size(b))
    integer  :: i, n

    ! check diagonal elements are not zero
    n = size(b)
    do i = 1, n
      if(abs(U(i,i)) > 0._r8) cycle
      call rexit('*** ERROR: zero diagonal element(s) (solvu) ***')
    end do

    ! solve equations
    x(n) = b(n) / U(n,n)
    do i = n-1, 1, -1
      x(i) = ( b(i) - dot_product(U(i,i+1:n), x(i+1:n)) ) / U(i,i)
    end do

  end function solvu

  function solvl(L, b) result(x)
    ! L is lower triangular
    implicit none
    real(r8), intent(in) :: L(:,:)
    real(r8), intent(in) :: b(:)
    real(r8)             :: x(size(b))
    integer              :: i, n

    ! check diagonal elements are not zero
    n = size(b)
    do i = 1, n
      if(abs(L(i,i)) > 0._r8) cycle
      call rexit('*** ERROR: zero diagonal element(s) (solvl) ***')
    end do

    ! solve equations
    x(1) = b(1) / L(1,1)
    do i = 2, n
      x(i) = ( b(i) - dot_product(L(i,1:i-1), x(1:i-1)) ) / L(i,i)
    end do

  end function solvl


  !-----------------------------------------------------------------------------
  ! computes the inverse of the matrix A using the LU decomposition
  ! depends on LAPACK
  ! see http://fortranwiki.org/fortran/show/Matrix+inversion

  function matinv(A) result(Ainv)
    real(r8), intent(in) :: A(:,:)
    real(r8) :: Ainv(size(A,1), size(A,2))
    real(r8) :: work(size(A,1))
    integer  :: piv(size(A,1))
    integer  :: n, info

    ! external LAPACK procedures
    external DGETRF
    external DGETRI

    ! check matrix is square
    n = size(A,1)
    if(n .ne. size(A,2)) then
      call rexit('*** ERROR: matrix is not square (matinv) ***')
    end if

    ! compute LU factorization with LAPACK subroutine DGETRF
    ! using partial pivoting with interchange rows
    Ainv = A
    call DGETRF(n, n, Ainv, n, piv, info)

    if(info /= 0) then
      call rexit('*** ERROR: singular matrix (matinv) ***')
    end if

    ! compute the inverse of the matrix using LAPACK subroutine DGETRI,
    ! using LU factorization computed by DGETRF
    call DGETRI(n, Ainv, n, piv, work, n, info)

    if(info /= 0) then
      call rexit('*** ERROR: matrix inversion failed (matinv) ***')
    end if

  end function matinv


  !-----------------------------------------------------------------------------
  ! computes the determinant of the matrix A using the LU decomposition

  real(r8) function det(A)
    implicit none
    real(r8), intent(in) :: A(:,:)
    real(r8) :: Ain(size(A,1), size(A,2))
    integer  :: piv(size(A,1))
    integer  :: i, n, info

    ! external LAPACK procedures
    external DGETRF

    ! check matrix is square
    n = size(A,1)
    if(n .ne. size(A,2)) then
      call rexit('*** ERROR: matrix is not square (matinv) ***')
    end if

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    Ain = A
    call DGETRF(n, n, Ain, n, piv, info)

    if(info /= 0) then
      call rexit('*** ERROR: LU decomposition failed (det) ***')
    end if

    ! compute determinant
    det = 1._r8
    do i = 1, n
      if(piv(i) /= i) then
        det = -det * Ain(i,i)
      else
        det = det * Ain(i,i)
      end if
    end do

  end function det


  !-----------------------------------------------------------------------------
  ! returns the cross-product A' * A of the matrix A

  function crossprod(A) result(AA)
    implicit none
    real(r8), intent(in) :: A(:,:)
    real(r8) :: AA(size(A,2), size(A,2))
    integer  :: i, j, n

    n = size(A,2)
    do j = 1, n  ! i <= j
      do i = 1, j
        AA(i,j) = dot_product(A(:,i), A(:,j))
        AA(j,i) = AA(i,j)
      end do
    end do

  end function crossprod


  !-----------------------------------------------------------------------------
  ! compute outer product p = a * a' or p = a * b'

  function outerprod1(a) result(p)
    real(r8), intent(in) :: a(:)
    real(r8)             :: p(size(a), size(a))

    p = spread(a, dim=2, ncopies=size(a))
    p = p * transpose(p)

  end function outerprod1

  function outerprod2(a, b) result(p)
    real(r8), intent(in) :: a(:)
    real(r8), intent(in) :: b(:)
    real(r8)             :: p(size(a), size(b))

    p = spread(a, dim=2, ncopies=size(b)) * &
        spread(b, dim=1, ncopies=size(a))

  end function outerprod2


end module matrix

