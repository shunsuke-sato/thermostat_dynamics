include "external/luxury.f90"
module global_variables
  implicit none
! math parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical system
  integer :: nsite, nelec
  real(8) :: t_hop, delta_gap
  real(8) :: KbT
  real(8),allocatable :: ham(:,:), ham0(:,:)
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: psi(:,:)
  real(8),allocatable :: rho_e(:)
  real(8) :: omega0, g_couple, gamma_damp
  real(8),allocatable :: xt(:), xt_n(:), vt(:), vt_n(:), vt_o(:)

! time-propagation
  real(8) :: tprop, dt
  integer :: nt


end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call initialize

  call calc_quantum_classical_ground_state


end program main
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: i,j

! set parameters
  nelec = 64
  nsite = 2*nelec

  t_hop     =  1d0
  delta_gap =  0d0
  
  omega0     =  0.1d0
  g_couple   =  0.2d0
  gamma_damp =  0.2d0*omega0

  KbT = 0.5d0

  tprop      =  2d0*pi*10000d0/omega0
  dt = 0.1d0
  nt = aint(tprop/dt)+1


  allocate(psi(nsite,nelec),zpsi(nsite,nelec))
  allocate(rho_e(nsite))
  allocate(xt(nsite), xt_n(nsite), vt(nsite), vt_n(nsite), vt_o(nsite))

  xt   = 0d0
  vt   = 0d0
  vt_o = 0d0
  

  allocate(ham(nsite,nsite), ham0(nsite, nsite))

  ham0 = 0d0
  do i = 1, nsite
    j = mod(i + 1 -1 + nsite, nsite) + 1
    ham0(i,j) = -t_hop
    j = mod(i - 1 -1 + nsite, nsite) + 1
    ham0(i,j) = -t_hop
  end do

  call initialize_random_number_generator



end subroutine initialize
!-------------------------------------------------------------------------------
subroutine calc_quantum_classical_ground_state
  use global_variables
  use luxury
  implicit none
  integer :: iscf, nscf
  integer :: istate
  integer :: i
  real(8) :: Eelec, Eion, Etot
!LAPACK ==
  real(8), allocatable :: amat(:,:)
  integer :: lwork
  real(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  lwork = 6*nsite**2
  allocate(work_lp(lwork))
  allocate(rwork(3*nsite-2))
  allocate(w(nsite))
!LAPACK ==

! initialize
!  xt = 0d0
  call ranlux_double(xt,nsite); xt = xt-0.5d0

  nscf = 300
  open(30,file='gs_conv.out')
  do iscf = 1, nscf

! calc ham
    ham = ham0
    do i = 1, nsite
      ham(i,i) = ham(i,i) - g_couple*xt(i)
    end do

    amat = ham
    call dsyev('V', 'U', nsite, amat(:,:), nsite &
      , w(:), work_lp(:), lwork, info)

    psi(1:nsite,1:nelec) = amat(1:nsite,1:nelec)

    Eelec = sum(w(1:nelec))
    Eion = 0.5d0*omega0**2*sum(xt**2)
    Etot = Eelec + Eion
    
    write(30,"(I7,2x,999e26.16e3)")iscf,Eelec,Eion,Etot
    write(*,"(I7,2x,999e26.16e3)")iscf,Eelec,Eion,Etot

    rho_e = 0d0
    do istate =1, nelec
      rho_e(:) = rho_e(:) + psi(:,istate)**2
    end do

    xt = g_couple*rho_e/omega0**2
    
  end do
  close(30)

  open(30,file='gs_data.out')
  do i=1,nsite
    write(30,"(I7,2x,999e26.16e3)")i,rho_e(i),xt(i)
  end do
  close(30)

end subroutine calc_quantum_classical_ground_state
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
subroutine gaussian_random_number_vec(rvec,nvec)
  implicit none
  integer,intent(in)  :: nvec
  real(8),intent(out) :: rvec(nvec)
  real(8) :: x1, x2
  integer :: i

  do i = 1, nvec/2
    call gaussian_random_number(x1,x2)
    rvec(2*i-1) = x1
    rvec(2*i)   = x2
  end do

  if(mod(nvec, 2) /= 0)then
    call gaussian_random_number(x1,x2)
    rvec(nvec) = x1
  end if

end subroutine gaussian_random_number_vec
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  
