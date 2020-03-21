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
  real(8),allocatable :: lambda_sp(:)
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: psi(:,:)
  real(8),allocatable :: rho_e(:)
  real(8) :: omega0, g_couple, gamma_damp
  real(8),allocatable :: xt(:), xt_n(:), vt(:), vt_n(:), vt_o(:)
  real(8) :: Egs_tot
 
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
!  call calc_electronic_canonical_ensemble

  call calc_quantum_classical_canonical_ensemble

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
  g_couple   =  0.1d0*omega0**2
  gamma_damp =  0.1d0*omega0

  KbT = 0.5d0 !0.5d0

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
  allocate(lambda_sp(nsite))

  ham0 = 0d0
  do i = 1, nsite
    j = mod(i + 1 -1 + nsite, nsite) + 1
    ham0(i,j) = -t_hop
    j = mod(i - 1 -1 + nsite, nsite) + 1
    ham0(i,j) = -t_hop
    ham0(i,i) = ham0(i,i) + (-1)**(i+1)*0.5d0*delta_gap
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

    lambda_sp(:) = w(:)
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

  Egs_tot = Etot

  open(30,file='gs_data.out')
  do i=1,nsite
    write(30,"(I7,2x,999e26.16e3)")i,rho_e(i),xt(i)
  end do
  close(30)

end subroutine calc_quantum_classical_ground_state
!-------------------------------------------------------------------------------
subroutine calc_electronic_canonical_ensemble
  use global_variables
  use luxury
  implicit  none
  integer :: nocc_dist(nsite)
  integer :: nlist_occ(nelec), nlist_unocc(nsite-nelec)
  real(8),allocatable :: prob(:)
  real(8) :: energy_diff, ss, rvec(1)
  integer :: nchoise, ichoise
  integer :: iter, niter, i, j, k
  real(8) :: Eelec_t, Eelec2_t, Egs
  real(8) :: Eelec_mb, Cv_mb
  integer :: njump
  real(8) :: occ_dist(nsite), mu_chem, mu_chem_max, mu_chem_min
  real(8) :: occ_dist_mb(nsite)
  real(8) :: tmp_nelec, eps_chk
  real(8) :: Eelec_sp, Cv_sp, Eelec2_sp

  niter = 1000000
  eps_chk = 1d-14

  njump = 0


  nchoise = nelec*(nsite-nelec)
  allocate(prob(0:nchoise))
  nocc_dist = 0
  nocc_dist(1:nelec) = 1
  Egs = sum(nocc_dist*lambda_sp)

  occ_dist = 0d0

  call pre_thermalization

  iter = 0
  open(20,file='scf_chk.out')
  do
    iter = iter + 1

    j = 0; k = 0
    do i = 1, nsite
      if(nocc_dist(i) == 0)then
        j = j + 1
        nlist_unocc(j) = i
      else if(nocc_dist(i) == 1)then
        k = k + 1
        nlist_occ(k) = i
      else
        stop 'error'
      end if
    end do

    prob(0) = 0d0
    k = 0
    do i = 1, nelec
      do j = 1, nsite-nelec
        k = k + 1
        energy_diff = lambda_sp(nlist_unocc(j))-lambda_sp(nlist_occ(i))
        prob(k) = prob(k-1) + 1d0 ! debug
!        prob(k) = prob(k-1) + exp(-0.5d0*energy_diff/KbT) ! debug
!        prob(k) = prob(k-1) + min(1d0, exp(-energy_diff/KbT)) ! debug

      end do
    end do
    if(k /= nchoise)stop 'error'


!    if(iter == 2)then
!    open(30,file='prob.out')
!    k = 0
!    write(30,*)k,prob(k),0d0
!    do i = 1, nelec
!      do j = 1, nsite-nelec
!        k = k + 1
!        energy_diff = lambda_sp(nlist_unocc(j))-lambda_sp(nlist_occ(i))
!        write(30,*)k,prob(k),energy_diff
!      end do
!    end do
!
!    close(30)
!    stop
!    end if

    ss = prob(nchoise)
    prob = prob/ss

    call  ranlux_double (rvec, 1)
    do ichoise = 0, nchoise
      if(prob(ichoise)>= rvec(1))exit
    end do
    if(ichoise > nchoise)stop 'error0'
!    write(*,*)ichoise,rvec(1),prob(ichoise)


    if(ichoise /= 0)then      
      njump = njump + 1
      k = 0
      do i = 1, nelec
        do j = 1, nsite-nelec
          k = k +1
          if(k==ichoise)then
            if(nocc_dist(nlist_occ(i)) /= 1) stop 'error1'
            if(nocc_dist(nlist_unocc(j)) /= 0) stop 'error2'
            call  ranlux_double (rvec, 1)
            energy_diff = lambda_sp(nlist_unocc(j))-lambda_sp(nlist_occ(i))
            if(rvec(1) < exp(-energy_diff/KbT))then
              nocc_dist(nlist_occ(i)) = 0
              nocc_dist(nlist_unocc(j)) = 1
            end if
          end if
        end do
      end do
    end if

    ss = sum(nocc_dist*lambda_sp)
    Eelec_t  = Eelec_t  + ss
    Eelec2_t = Eelec2_t + ss**2
    write(20,"(I7,2x,999e26.16e3)")iter,(Eelec_t/iter-Egs)/nsite,Eelec2_t/iter &
      ,((Eelec2_t/iter)-(Eelec_t/iter)**2)/(KbT**2)

    occ_dist = occ_dist + nocc_dist
    if(niter == iter)exit
  end do
  close(20)
  write(*,*)"njump=",njump
  write(*,"(A,2x,999e26.16e3)")"Eelec =",(Eelec_t/niter-Egs)/nsite
  write(*,"(A,2x,999e26.16e3)")"Eelec2=",Eelec2_t/niter
  write(*,"(A,2x,999e26.16e3)")"C_t   =",((Eelec2_t/iter)-(Eelec_t/iter)**2)/(KbT**2)
  Eelec_mb = Eelec_t/niter-Egs
  Cv_mb    = ((Eelec2_t/iter)-(Eelec_t/iter)**2)/(KbT**2)

  occ_dist_mb = occ_dist/niter

! single-particle
  mu_chem_min = lambda_sp(1)
  mu_chem_max = lambda_sp(nsite)
  do
    mu_chem = 0.5d0*(mu_chem_max + mu_chem_min)
    do i = 1, nsite
      occ_dist(i) = 1d0/(exp((lambda_sp(i)-mu_chem)/KbT)+1d0)
    end do
    tmp_nelec = sum(occ_dist)
    if(tmp_nelec > dble(nelec))then
      mu_chem_max = mu_chem
    else
      mu_chem_min = mu_chem
    end if
    if((mu_chem_max - mu_chem_min) < eps_chk)exit

  end do
  mu_chem = 0.5d0*(mu_chem_max + mu_chem_min)
  do i = 1, nsite
    occ_dist(i) = 1d0/(exp((lambda_sp(i)-mu_chem)/KbT)+1d0)
  end do

  open(30,file="occ_dist.out")
  do i = 1, nsite
    write(30,"(I7,2x,99e26.16e3)")i,lambda_sp(i),occ_dist_mb(i),occ_dist(i)
  end do
  close(30)  

  Eelec_sp  = sum(occ_dist*lambda_sp)
  Eelec2_sp = sum(occ_dist*lambda_sp**2)
  Cv_sp = (Eelec2_sp-Eelec_sp**2)/(kbT**2)
  write(*,"(A)")"Single-particle evaluation"
  write(*,"(A,2x,e26.16e3,2x,I7)")"num elecc =", sum(occ_dist),nelec
  write(*,"(A,2x,999e26.16e3)")"Eelec =",(Eelec_sp-Egs)/nsite
  write(*,"(A,2x,999e26.16e3)")"Eelec2=",Eelec2_sp
  write(*,"(A,2x,999e26.16e3)")"C_t   =",Cv_sp

  contains
    subroutine pre_thermalization
      implicit none
      
      do iter = 1,nsite

        j = 0; k = 0
        do i = 1, nsite
          if(nocc_dist(i) == 0)then
            j = j + 1
            nlist_unocc(j) = i
          else if(nocc_dist(i) == 1)then
            k = k + 1
            nlist_occ(k) = i
          else
            stop 'error'
          end if
        end do

        prob(0) = 1d0
        k = 0
        do i = 1, nelec
          do j = 1, nsite-nelec
            k = k + 1
            energy_diff = lambda_sp(nlist_unocc(j))-lambda_sp(nlist_occ(i))
            prob(k) = prob(k-1) + exp(-0.5d0*energy_diff/KbT) ! debug
!            prob(k) = prob(k-1) + min(1d0, exp(-energy_diff/KbT)) ! debug

          end do
        end do
        if(k /= nchoise)stop 'error'


!    if(iter == 2)then
!    open(30,file='prob.out')
!    k = 0
!    write(30,*)k,prob(k),0d0
!    do i = 1, nelec
!      do j = 1, nsite-nelec
!        k = k + 1
!        energy_diff = lambda_sp(nlist_unocc(j))-lambda_sp(nlist_occ(i))
!        write(30,*)k,prob(k),energy_diff
!      end do
!    end do
!
!    close(30)
!    stop
!    end if

        ss = prob(nchoise)
        prob = prob/ss

        call  ranlux_double (rvec, 1)
        do ichoise = 0, nchoise
          if(prob(ichoise)>= rvec(1))exit
        end do
        if(ichoise > nchoise)stop 'error0'
!    write(*,*)ichoise,rvec(1),prob(ichoise)


        if(ichoise /= 0)then      
          njump = njump + 1
          k = 0
          do i = 1, nelec
            do j = 1, nsite-nelec
              k = k +1
              if(k==ichoise)then
                if(nocc_dist(nlist_occ(i)) /= 1) stop 'error1'
                if(nocc_dist(nlist_unocc(j)) /= 0) stop 'error2'
                nocc_dist(nlist_occ(i)) = 0
                nocc_dist(nlist_unocc(j)) = 1
              end if
            end do
          end do
        end if

      end do


    end subroutine pre_thermalization

end subroutine calc_electronic_canonical_ensemble
!-------------------------------------------------------------------------------
subroutine calc_quantum_classical_canonical_ensemble
  use global_variables
  implicit none
  real(8) :: x0(nsite)
  real(8) :: dx(nsite),pt(nsite)
  real(8) :: ss
  integer :: ncount
  integer :: i,j,k
  integer :: i0,j0,k0
  integer :: nocc_dist(nsite)
  integer :: nlist_occ(nelec), nlist_unocc(nsite-nelec)
  integer :: isample,nsample
  integer :: iter
  real(8) :: Eelec_t, Eph_t, Etot_t
  real(8) :: Eelec, Eph, Etot, Etot2, cv
  real(8) :: rvec(1), energy_diff
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
  allocate(amat(nsite,nsite))
!LAPACK ==

  nsample = 1024

  nocc_dist = 0
  nocc_dist(1:nelec) = 1

  x0 = g_couple*rho_e/omega0**2

  Eelec = 0d0
  Eph = 0d0
  Etot =0d0
  Etot2 = 0d0
  ncount = 0
  open(30,file='convergence.out')
  do isample = 1, nsample
    call gaussian_random_number_vec(dx,nsite)
    call gaussian_random_number_vec(pt,nsite)
    dx = sqrt(KbT)/omega0*dx ! debug
    pt = sqrt(KbT)*pt ! debug
    xt = x0 + dx

    ss = g_couple*sum(rho_e*xt)/nelec
    ham = ham0
    do i = 1, nsite
      ham(i,i) = ham(i,i) - g_couple*xt(i) + ss
    end do


    amat = ham
    call dsyev('V', 'U', nsite, amat(:,:), nsite &
      , w(:), work_lp(:), lwork, info)

    lambda_sp(:) = w(:)
    psi(1:nsite,1:nelec) = amat(1:nsite,1:nelec)

    do iter = 1, nsite*128

      j = 0; k = 0
      do i = 1, nsite
        if(nocc_dist(i) == 0)then
          j = j + 1
          nlist_unocc(j) = i
        else if(nocc_dist(i) == 1)then
          k = k + 1
          nlist_occ(k) = i
        else
          stop 'error'
        end if
      end do

      call  ranlux_double (rvec, 1)
      rvec(1) = rvec(1)*nelec
      k0 = aint(rvec(1))
      k0 = mod(k0,nelec)+1
      call  ranlux_double (rvec, 1)
      rvec(1) = rvec(1)*(nsite-nelec)
      j0 =aint(rvec(1))
      j0 = mod(j0,(nsite-nelec))+1
      energy_diff = lambda_sp(nlist_unocc(j0))-lambda_sp(nlist_occ(k0))
      call  ranlux_double (rvec, 1)

      if(exp(-energy_diff/KbT) > rvec(1))then
        if(nocc_dist(nlist_occ(k0)) /= 1) stop 'error1'
        if(nocc_dist(nlist_unocc(j0)) /= 0) stop 'error2'
        nocc_dist(nlist_occ(k0)) = 0
        nocc_dist(nlist_unocc(j0)) = 1
      end if

      if(mod(iter,nsite)==0)then
        ncount = ncount + 1

        Eelec_t = sum(nocc_dist*lambda_sp) 
        Eph_t   = sum(0.5d0*pt**2+0.5d0*omega0**2*dx**2+0.5d0*omega0**2*x0**2 &
          -g_couple*x0*rho_e)
!        Eph_t   = sum(0.5d0*pt**2+0.5d0*omega0**2*dx**2) ! debug
        Etot_t = Eelec_t + Eph_t - Egs_tot

        Eelec = Eelec + Eelec_t
        Eph   = Eph   + Eph_t
        Etot  = Etot  + Etot_t
        Etot2 = Etot2 + Etot_t**2
        write(30,"(I7,2x,99e26.16e3)")ncount,Eelec/ncount,Eph/ncount,Etot/ncount &
          ,(Etot2/ncount-(Etot/ncount)**2)/KbT**2
      end if

      
    end do

  end do
  close(30)

  Eelec = Eelec/ncount
  Eph   = Eph/ncount
  Etot  = Etot/ncount
  Etot2 = Etot2/ncount
  cv = (Etot2-Etot**2)/KbT**2
  write(*,"(A,2x,999e26.16e3)")'Eelec=',Eelec/nsite
  write(*,"(A,2x,999e26.16e3)")'Eph  =',Eph/nsite
  write(*,"(A,2x,999e26.16e3)")'Etot =',Etot/nsite
  write(*,"(A,2x,999e26.16e3)")'cv   =',cv/nsite

end subroutine calc_quantum_classical_canonical_ensemble
!-------------------------------------------------------------------------------
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
  
