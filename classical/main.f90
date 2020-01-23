include "../external/luxury.f90"
module global_variables
  implicit none

  integer :: num_particle
  real(8),allocatable :: xt(:)
  real(8),allocatable :: vt(:),vt_n(:),vt_o(:)

  real(8) :: KbT, gamma, omega

! propagation
  integer :: nt
  real(8) :: Tprop, dt


end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call initialize_random_number_generator
  call initialize_parameters
  call dynamics

end program main
!-------------------------------------------------------------------------------
subroutine initialize_parameters
  use global_variables
  implicit none

  num_particle = 64
  allocate(xt(num_particle), vt(num_particle))
  allocate(vt_n(num_particle), vt_o(num_particle))

  KbT = 1d0
  omega = 1d0
  gamma = 0.05d0*omega

  xt = 0d0
  vt_o = 0d0
  vt = 0d0


  Tprop = 1000d0/omega
  dt = 0.1d0
  nt = aint(Tprop/dt)+1


end subroutine initialize_parameters
!-------------------------------------------------------------------------------
subroutine dynamics
  use global_variables
  implicit none
  integer :: it, ip, ibin
  real(8),allocatable :: psi(:)
  real(8) :: x1, x2
  integer :: nbin
  real(8) :: xmin, xmax, dx, vmin, vmax, dv, ss
  real(8),allocatable :: rho_x(:), rho_v(:)
  character(512) :: cit, cfilename

  nbin = 128
  xmin = -10d0
  xmax = 10d0
  vmin = -10d0
  vmax = 10d0
  dx = (xmax-xmin)/nbin
  dv = (vmax-vmin)/nbin

  allocate(psi(num_particle))
  allocate(rho_x(nbin), rho_v(nbin))
  rho_x = 0d0
  rho_v = 0d0

  do it = 0, nt

    do ip = 1, num_particle
      call gaussian_random_number(x1,x2)
      psi(ip) = x1
    end do
    vt_n = vt_o*exp(-gamma*dt) -omega*xt*dt + sqrt(2d0*KbT*gamma*dt)*psi
    vt = 0.5d0*(vt_n + vt_o)

    do ip = 1, num_particle
      ibin = aint((xt(ip)-xmin)/dx)
      ibin = max(ibin,1); ibin = min(ibin,nbin)
      rho_x(ibin) = rho_x(ibin) + 1d0

      ibin = aint((vt(ip)-vmin)/dv)
      ibin = max(ibin,1); ibin = min(ibin,nbin)
      rho_v(ibin) = rho_v(ibin) + 1d0

    end do

    xt = xt + dt*vt_n
    vt_o = vt_n
    
    if(mod(it, 1000)==0)then

      write(cit,"(I10.10)")it
      cfilename = 'rho_x_'//trim(cit)//'.out'
      open(30,file=cfilename)
      ss = sum(rho_x)*dx
      do ibin = 1, nbin
        write(30,"(999e16.6e3)")xmin+ibin*dx+0.5d0*dx,rho_x(ibin)/ss
      end do
      close(30)

      cfilename = 'rho_v_'//trim(cit)//'.out'
      open(30,file=cfilename)
      ss = sum(rho_v)*dv
      do ibin = 1, nbin
        write(30,"(999e16.6e3)")vmin+ibin*dv+0.5d0*dv,rho_v(ibin)/ss
      end do
      close(30)

    end if



  end do

end subroutine dynamics
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine initialize_random_number_generator
  use luxury
  implicit none
! parameters for random_number generator
  integer :: lux_ran = 3, K1_ran = 0, K2_ran = 0
  integer :: INT_ran

!  INT_ran = myrank**2 + 1000*myrank + 1
!  INT_ran = myrank**2 + 1000*myrank + 10001
!  INT_ran = 19900126 + 1234567
  INT_ran = 1234567
  CALL RLUXGO(lux_ran,int_ran,K1_ran,K2_ran)


end subroutine initialize_random_number_generator
!-------------------------------------------------------------------------------
subroutine ranlux_double(rvec,len)
  use luxury
  implicit none
  integer,intent(in) :: len
  real(8),intent(out) :: rvec(len)
  integer :: len4 = 2
  real :: rvec4(2)
  integer(8) :: int1,i


  do i = 1,len

    CALL RANLUX (rvec4, len4)

    int1 = aint(rvec4(1)*1d6)*10000000 + aint(rvec4(2)*1d7)
    rvec(i) =dble(int1)*1d-13

  end do

end subroutine ranlux_double
!-------------------------------------------------------------------------------
subroutine gaussian_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  integer :: len = 2
  real(8) :: rvec(2),tmp


  CALL ranlux_double (rvec, len)
!
  if(rvec(1) == 0d0)then
    x1 = 0d0
    x2 = 0d0
  else 
    tmp = sqrt(-2d0*log(rvec(1)))
    x1 = tmp*cos(2d0*pi*rvec(2))
    x2 = tmp*sin(2d0*pi*rvec(2))
  end if

end subroutine gaussian_random_number
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
