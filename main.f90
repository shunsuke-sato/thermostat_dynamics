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

! block average
  integer :: nblock
  real(8) :: t_relax


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
  nelec = 64
  nsite = 2*nelec

  t_hop     =  1d0
  delta_gap =  1d0
  
  omega0     =  0.1d0
  g_couple   =  0.1d0*omega0  ! debug
  gamma_damp =  0.1d0*omega0 ! debug

  KbT = 3d0
  open(30,file='inp_tmp')
  read(30,*)KbT
  close(30)

  tprop      =  10000d0 !2d0*pi*10000d0/omega0
  dt = 0.1d0
  nt = aint(tprop/dt)+1

  nblock = 10


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
    ham0(i,i) = ham0(i,i) + (-1)**(i+1)*0.5d0*delta_gap
  end do

  call initialize_random_number_generator

  call set_gs_wavefunction


end subroutine initialize
!-------------------------------------------------------------------------------
subroutine set_gs_wavefunction
  use global_variables
  implicit none
  integer :: i
  complex(8) :: zhpsi_t(nsite,nelec), zhpsi2_t(nsite,nelec)
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

!  zhpsi_t = matmul(ham0,zpsi)
!  zhpsi2_t = matmul(ham0,zhpsi_t)
!  write(*,*)sum(conjg(zpsi)*zhpsi_t),sum(conjg(zpsi)*zhpsi2_t)
!  write(*,*)sum(conjg(zpsi)*zhpsi2_t)-sum(conjg(zpsi)*zhpsi_t)**2
!  do i = 1, nelec
!    write(*,*)i,w(i),sum(conjg(zpsi(:,i))*zhpsi_t(:,i))
!    write(*,*)i,w(i)**2,sum(conjg(zpsi(:,i))*zhpsi2_t(:,i))
!  end do
!
!  do i = 1, nelec
!    zhpsi2_t(:,i) = zhpsi_t(:,i) - w(i)*zpsi(:,i)
!  end do
!  write(*,*)sum(abs(zhpsi2_t)**2)
!  stop

end subroutine set_gs_wavefunction
!-------------------------------------------------------------------------------
subroutine Langevin_dynamics
  use global_variables
  implicit none
  integer :: it, i,j
  real(8),allocatable :: xi(:), force(:)
  integer :: icount
  real(8) :: Etot, Eelec, Eion, Ecoup
  real(8) :: Etot_t, Eelec_t, Eion_t
  real(8) :: Etot2
  real(8) :: Etot2_t, Eelec2_t
  real(8) :: rvec1(nsite),rvec2(nsite)
  real(8) :: ss
  complex(8) :: zs
  logical :: if_file_exists

  integer :: iblock
  real(8) :: ss_ave, ss_sigma
  real(8),allocatable :: Eelec_bave(:), Eion_bave(:)
  real(8),allocatable :: Etot_bave(:), cv_bave(:)
  real(8) :: results_data(99)

  complex(8),allocatable :: zhpsi_t(:,:),zhpsi2_t(:,:)
  complex(8),allocatable :: zham2(:,:)
  real(8),allocatable :: eps_sp_t(:)
  
  allocate(xi(nsite), force(nsite))
  allocate(zham2(nelec,nelec))
  allocate(zhpsi_t(nsite,nelec))
  allocate(zhpsi2_t(nsite,nelec))
  allocate(eps_sp_t(nelec))

  allocate(Eelec_bave(0:nblock), Eion_bave(0:nblock), Etot_bave(0:nblock), cv_bave(0:nblock))

  inquire(file="checkpoint.out",exist=if_file_exists)

  if(if_file_exists)then
    open(40,file="checkpoint.out",form='unformatted')
    read(40)zpsi
    read(40)xt_n,xt
    read(40)vt_n,vt_o
    close(40)

  else
!    do i = 1, nelec
!      call ranlux_double(rvec1,nsite)
!      call ranlux_double(rvec2,nsite)
!      zpsi(:,i)=rvec1(:)*exp(zi*2d0*pi*rvec2(:))
!    end do

! set to GS state
    vt = 0d0
    vt_o = 0d0
    call calc_ground_state_for_whole_system

  end if

  do i = 1, nelec
    ss = sum(abs(zpsi(:,i))**2)
    zpsi(:,i) = zpsi(:,i)/sqrt(ss)
    do j = 1, i-1
      zs = sum(conjg(zpsi(:,j))*zpsi(:,i))
      zpsi(:,i) = zpsi(:,i) -zs*zpsi(:,j)
    end do
  end do

  block_ave: do iblock = 0, nblock

  Etot=0d0
  Etot2 = 0d0
  Eelec=0d0
  Eion=0d0



  open(30,file='energy_t.out')

  icount = 0
  do it = 0, nt

    if(mod(it, max(nt/200, 1)) == 0)write(*,*)"it/nt=",it,dble(it)/dble(nt)
    call gaussian_random_number_vec(xi,nsite)
    call calc_density
    force(:) = g_couple*rho_e(:) - omega0**2*xt
    vt_n = vt_o*exp(-gamma_damp*dt) + force*dt + sqrt(2d0*KbT*gamma_damp*dt)*xi
    vt = 0.5d0*(vt_n + vt_o)

! start: calculate observables
    icount = icount + 1
    ham = ham0
    do i = 1, nsite
      ham(i,i) = ham(i,i) - g_couple*xt_n(i)
    end do

    Eion_t = 0.5d0*sum(vt**2+omega0**2*xt**2)

    zhpsi_t = matmul(ham,zpsi)
    zhpsi2_t = matmul(ham,zhpsi_t)

    Eelec_t = 0d0
    do i = 1,nelec
      zham2(i,i) = sum(conjg(zpsi(:,i))*zhpsi_t(:,i))
      Eelec_t = Eelec_t + zham2(i,i)
      do j = i+1,nelec
        zham2(j,i) = sum(conjg(zpsi(:,j))*zhpsi_t(:,i))
        zham2(i,j) = conjg(zham2(j,i))
      end do
    end do


!    Eelec2_t = sum(conjg(zpsi)*zhpsi2_t)
!    do i = 1,nelec
!      do j = i+1,nelec
!        Eelec2_t = Eelec2_t + 2d0*real(zham2(i,i))*real(zham2(j,j)) &
!          -2d0*abs(zham2(i,j))**2
!      end do
!    end do

    Eelec2_t = 0d0
    do i = 1,nelec
      do j = i+1,nelec
        Eelec2_t = Eelec2_t + 2d0*real(zham2(i,i))*real(zham2(j,j)) &
          -2d0*abs(zham2(i,j))**2
      end do
   end do
   Eelec2_t = 2d0*Eelec2_t/dble(nelec-1)
    Eelec2_t = Eelec2_t + sum(conjg(zpsi)*zhpsi2_t)    

!    write(*,*)Eelec_t,Eelec2_t,Eelec2_t-Eelec_t**2


    Etot_t = Eelec_t + Eion_t
    Etot2_t = Eelec2_t + 2d0*Eelec_t*Eion_t + Eion_t**2
!    Etot_t =  Eion_t ! debug
!    Etot2_t = Eion_t**2 ! debug

    Eelec = Eelec + Eelec_t
    Eion  = Eion  + Eion_t
    Etot  = Etot  + Etot_t
    Etot2  = Etot2  + Etot2_t

    write(30,"(999e26.16e3)")dt*it,Eelec/icount,Eion/icount,Etot/icount&
      ,((Etot2/icount-(Etot/icount)**2)/nsite)/KbT**2

! end  : calculate observables



    xt_n = xt + dt*vt_n
    call dt_evolve_elec
    xt = xt_n
    vt_o = vt_n

  end do

  close(30)

  Eelec_bave(iblock)= Eelec/icount
  Eion_bave(iblock) = Eion/icount
  Etot_bave(iblock) = Etot/icount
  cv_bave(iblock) = ((Etot2/icount-(Etot/icount)**2)/nsite)/KbT**2

  end do block_ave

  
  ss_ave = sum(Eelec_bave(1:nblock))/nblock
  ss_sigma = sqrt(sum((Eelec_bave(1:nblock)-ss_ave)**2  )/(nblock - 1))
  write(*,"(A,2x,999e26.16e3)")"Eelec              =",ss_ave
  write(*,"(A,2x,999e26.16e3)")"standard deviation =",ss_sigma
  write(*,"(A,2x,999e26.16e3)")"standard error     =",ss_sigma/sqrt(dble(nblock))
  results_data(1) = kbt
  results_data(2) = ss_ave
  results_data(3) = ss_sigma/sqrt(dble(nblock))

  ss_ave = sum(Eion_bave(1:nblock))/nblock
  ss_sigma = sqrt(sum((Eion_bave(1:nblock)-ss_ave)**2  )/(nblock - 1))
  write(*,"(A,2x,999e26.16e3)")"Eion               =",ss_ave
  write(*,"(A,2x,999e26.16e3)")"standard deviation =",ss_sigma
  write(*,"(A,2x,999e26.16e3)")"standard error     =",ss_sigma/sqrt(dble(nblock))
  results_data(4) = ss_ave
  results_data(5) = ss_sigma/sqrt(dble(nblock))

  ss_ave = sum(Etot_bave(1:nblock))/nblock
  ss_sigma = sqrt(sum((Etot_bave(1:nblock)-ss_ave)**2  )/(nblock - 1))
  write(*,"(A,2x,999e26.16e3)")"Etot               =",ss_ave
  write(*,"(A,2x,999e26.16e3)")"standard deviation =",ss_sigma
  write(*,"(A,2x,999e26.16e3)")"standard error     =",ss_sigma/sqrt(dble(nblock))
  results_data(6) = ss_ave
  results_data(7) = ss_sigma/sqrt(dble(nblock))

  ss_ave = sum(cv_bave(1:nblock))/nblock
  ss_sigma = sqrt(sum((cv_bave(1:nblock)-ss_ave)**2  )/(nblock - 1))
  write(*,"(A,2x,999e26.16e3)")"cv                 =",ss_ave
  write(*,"(A,2x,999e26.16e3)")"standard deviation =",ss_sigma
  write(*,"(A,2x,999e26.16e3)")"standard error     =",ss_sigma/sqrt(dble(nblock))
  results_data(8) = ss_ave
  results_data(9) = ss_sigma/sqrt(dble(nblock))


  open(40,file='results_data.out')
  write(40,"(999e26.16e3)")results_data(1:9)
  close(40)

  open(40,file="checkpoint.out",form='unformatted')
  write(40)zpsi
  write(40)xt_n,xt
  write(40)vt_n,vt_o
  close(40)


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
subroutine calc_ground_state_for_whole_system
  use global_variables
  implicit none
  integer :: i, iscf, nscf
  complex(8) :: zhpsi_t(nsite,nelec), zhpsi2_t(nsite,nelec)
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
  
  nscf = 100

  do i = 1, nsite
    xt(i) = sin(dble(i))
  end do

  write(*,"(A)")"Compute the ground state"
  do iscf = 1, nscf
    write(*,*)"iscf=",iscf
    ham = ham0
    do i = 1, nsite
      ham(i,i) = ham(i,i) - g_couple*xt(i)
    end do

    amat = ham

    call dsyev('V', 'U', nsite, amat(:,:), nsite &
      , w(:), work_lp(:), lwork, info) 

    zpsi(:,1:nelec) = amat(:,1:nelec)
    call calc_density
    write(*,*)"scf-error",sum((xt-g_couple*rho_e/omega0**2)**2)/nsite
    xt = g_couple*rho_e/omega0**2

  end do

  open(30,file='gs_rho_xn.out')
  do i = 1, nsite
    write(30,"(I7,2x,999e26.16e3)")i,rho_e(i),xt(i)
  end do
  close(30)


end subroutine calc_ground_state_for_whole_system
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

