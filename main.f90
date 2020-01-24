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
  call initialize_random_number_generator

  call Langevin_dynamics


end program main
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: i,j

! set parameters
  nelec = 8
  nsite = 2*nelec

  t_hop     =  1d0
  delta_gap =  0d0
  
  omega0     =  0.1d0
  g_couple   =  0.1d0
  gamma_damp =  0.1d0*omega0

  KbT = 0.5d0

  tprop      =  2d0*pi*10000d0/omega0
  dt = 0.1d0
  nt = aint(tprop/dt)+1


  allocate(zpsi(nsite,nelec))
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

  call set_gs_wavefunction


end subroutine initialize
!-------------------------------------------------------------------------------
subroutine set_gs_wavefunction
  use global_variables
  implicit none
  integer :: i
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

  allocate(amat(nsite,nsite))
  amat = ham0

  call dsyev('V', 'U', nsite, amat(:,:), nsite &
    , w(:), work_lp(:), lwork, info) 

  open(20,file='eigen_val_gs.out')
  do i = 1, nsite
    write(20,"(I7,2x,999e26.16e3)")i,w(i)
  end do
  close(20)

  open(20,file='eigen_wf_gs.out')
  do i = 1, nsite
    write(20,"(I7,2x,999e26.16e3)")i,amat(i, 1:nsite)
  end do
  close(20)

  zpsi(:,1:nelec) = amat(:,1:nelec)


end subroutine set_gs_wavefunction
!-------------------------------------------------------------------------------
subroutine Langevin_dynamics
  use global_variables
  implicit none
  integer :: it
  real(8),allocatable :: xi(:), force(:)
  integer :: icount
  real(8) :: Etot, Eelec, Eion, Ecoup
  real(8) :: Etot_t, Eelec_t, Eion_t, Ecoup_t
  
  allocate(xi(nsite), force(nsite))

  Etot=0d0
  Eelec=0d0
  Eion=0d0
  Ecoup=0d0

  open(30,file='energy_t.out')

  icount = 0
  do it = 0, nt

    if(mod(it, max(nt/200, 1)) == 0)write(*,*)"it/nt=",it,dble(it)/dble(nt)
    call gaussian_random_number_vec(xi,nsite)
    call calc_density
    force(:) = g_couple*rho_e(:) - omega0*xt
    vt_n = vt_o*exp(-gamma_damp*dt) + force*dt + sqrt(2d0*KbT*gamma_damp*dt)*xi
    vt = 0.5d0*(vt_n + vt_o)

! start: calculate observables
    icount = icount + 1
    Eion_t = 0.5d0*sum(vt**2+omega0**2*xt**2)/nsite
    Eelec_t = sum(conjg(zpsi)*matmul(ham0,zpsi))/nsite
    Ecoup_t = - g_couple*sum(rho_e(:)*xt(:))/nsite
    Etot_t = Eelec_t + Ecoup_t + Eion_t

    Eelec = Eelec + Eelec_t
    Ecoup = Ecoup + Ecoup_t
    Eion  = Eion  + Eion_t
    Etot  = Etot  + Etot_t

    write(30,"(999e26.16e3)")dt*it,Eelec/icount,Ecoup/icount,Eion/icount,Etot/icount

! end  : calculate observables



    xt_n = xt + dt*vt_n
    call dt_evolve_elec
    xt = xt_n
    vt_o = vt_n

  end do

  close(30)


end subroutine Langevin_dynamics
!-------------------------------------------------------------------------------
subroutine calc_density
  use global_variables
  implicit none
  integer :: istate

  rho_e = 0d0
  do istate = 1, nelec
    rho_e(:) = rho_e(:) + abs(zpsi(:,istate))**2
  end do

!  write(*,*)'num_elec. = ', sum(rho_e)

end subroutine calc_density
!-------------------------------------------------------------------------------
subroutine dt_evolve_elec
  use global_variables
  implicit none
  integer :: i,j
  complex(8),allocatable :: zUprop(:,:)
!LAPACK ==
  real(8), allocatable :: amat(:,:)
  integer :: lwork
  real(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  allocate(amat(nsite,nsite))
  lwork = 6*nsite**2
  allocate(work_lp(lwork))
  allocate(rwork(3*nsite-2))
  allocate(w(nsite))
!LAPACK ==

  allocate(zUprop(nsite,nsite))

  ham = ham0
  do i = 1, nsite
    ham(i,i) = ham(i,i) - g_couple*0.5d0*(xt(i) + xt_n(i))
  end do

  amat = ham
  call dsyev('V', 'U', nsite, amat(:,:), nsite &
    , w(:), work_lp(:), lwork, info) 


  zUprop(:,:) = 0d0
  do i = 1, nsite
    zUprop(i,i) = exp(-zI*dt*w(i))
  end do
  
  zUprop = matmul(amat, matmul(zUprop, transpose(amat)))
  zpsi = matmul(zUprop,zpsi)

end subroutine dt_evolve_elec
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

