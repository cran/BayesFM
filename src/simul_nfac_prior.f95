subroutine simnfacprior(nmeas, Kmax, Nid, kappa, nrep, seed, nfac, restrid)

  use global
  use probability, only : set_seed, rdirich, runif
  implicit none

  integer,  intent(in)  :: nmeas
  integer,  intent(in)  :: Kmax
  integer,  intent(in)  :: Nid
  real(r8), intent(in)  :: kappa(Kmax)
  integer,  intent(in)  :: nrep
  integer,  intent(in)  :: seed
  integer,  intent(out) :: nfac(nrep)
  logical,  intent(out) :: restrid(nrep)

  real(r8) :: prob(Kmax)
  integer  :: dedic(nmeas)
  integer  :: ndedic(Kmax)
  real(r8) :: e, csum
  logical  :: checkid
  integer  :: i, j, k

  checkid = Nid > 1
  restrid = .true.

  call set_seed(seed)

  do i = 1, nrep

    ! sample indicator probabilities
    prob = rdirich(kappa)

    ! sample indicators
    ndedic = 0
    do j = 1, nmeas
      e = runif()
      csum = 0._r8
      do k = 1, Kmax
        csum = csum + prob(k)
        if(e <= csum) exit
      end do
      dedic(j) = k
      ndedic(k) = ndedic(k) + 1
    end do

    ! count number of factors
    nfac(i) = count(ndedic > 0)

    ! check identification restriction
    if(checkid) restrid(i) = all(ndedic >= Nid .or. ndedic == 0)

  end do

end subroutine simnfacprior
