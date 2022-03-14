
module ranmod1
    integer,parameter::i4b=selected_int_kind(9)
    integer,parameter::k15=selected_int_kind(15)
    integer,parameter::dp=selected_real_kind(8)
end module ranmod1
   
module ranmod2
    use ranmod1
    real(dp)::U(97),cc,cd,cm
    integer(k15)::i97,j97
end module ranmod2
   
module ranutil
    use ranmod1
    integer(i4b)::rand_inst=0
    integer(i4b)::rand_feedback=1
   contains
    subroutine initrandom(i,i2)
     use ranmod1
     implicit none
     integer(k15),optional,intent(in)::i,i2
     integer(k15)::seed_in,kl,ij
     character(len=10)::fred
     real(dp)::klr
     
     if (present(i)) then
      seed_in = i
     else
      seed_in = -1
     end if
     if (seed_in /= -1) then
      if (present(i2)) then
       kl = i2
       if (i2 > 30081) stop 'i2 > 30081'
      else
       kl = 9373
      end if
      ij = i
     else
      call system_clock(count=ij)
      ij = mod(ij+rand_inst*100,31328)
      call date_and_time(time=fred)
      read(fred,'(e10.3)')klr
      kl = mod(int(klr*1000),30081)
     end if
     call rmarin(ij,kl)
    end subroutine initrandom
    
    subroutine rmarin(ij,kl)
     use ranmod1
     use ranmod2
     implicit none
     integer(k15),intent(in)::ij,kl
     integer(k15)::i,j,k,l,ii,jj,m
     real(dp)::s,t
     
     if (ij<0 .or. ij>31328 .or. kl<0 .or. kl>30081) then
      stop 'Erro ao iniciar o gerador de números aleatórios'
     end if
     i = mod(ij/177,177)+2
     j = mod(ij,177)+2
     k = mod(kl/169,178)+1
     l = mod(kl,169)
     do ii = 1,97
      s = 0.0d0
      t = 0.5d0
      do jj = 1,24
       m = mod(mod(i*j,179)*k,179)
       i = j
       j = k
       k = m
       l = mod(53*l+1,169)
       if (mod(l*m,64)>=32) then
        s = s + t
       end if
       t = 0.5d0*t
      end do
      u(ii) = s
     end do
     cc = 362436.0d0/16777216.0d0
     cd = 7654321.0d0/16777216.0d0
     cm = 16777213.0d0/16777216.0d0
     i97 = 97
     j97 = 33
     
     return
    end subroutine rmarin
    
    function ranmar() result(fn_val)
     use ranmod1
     use ranmod2
     implicit none
     real(dp) :: fn_val
     real(dp) :: uni
     
     uni =u(i97) -u(j97)
     if ( uni < 0.0d0 ) uni = uni + 1.0d0
     u(i97) = uni
     i97 = i97 - 1
     if (i97 == 0) i97 = 97
     j97 = j97 - 1
     if (j97 == 0) j97 = 97
     cc = cc - cd
     if ( cc < 0.0d0 ) cc = cc + cm
     uni = uni - cc
     if ( uni < 0.0d0 ) uni = uni + 1.0d0
     fn_val = uni
     
     return
    end function ranmar
    
    real(dp) function rand_real(a,b)
     use ranmod1
     implicit none
     real(dp),intent(in)::a,b
     
     rand_real = (b-a)*ranmar() + a
    end function rand_real
    
    integer(i4b) function rand_int(i,j)
     use ranmod1
     implicit none
     integer(i4b),intent(in)::i,j
     
     rand_int = int(j*ranmar()) + i
    end function rand_int
    
end module ranutil
 
   !-----------------------------------------------------------------!
   !                          MARSAGLIA                              !
   ! Usage example:                                                  !
   !    use ranutil                                                  !
   !    ...                                                          !
   !    integer :: i                                                 !
   !    real(rk) :: x,g                                              !
   !    ...                                                          !                                           !
   !    call initrandom()                                            !
   !    x=ranmar()                                                   !
   !    x=rand_real(25.,59.)                                         !
   !    i=rand_int(0,100)                                            !
   !    ...                                                          !
   !                                                                 !
   !-----------------------------------------------------------------!
 

module var_FS
    integer :: Ns,N_temp,N_mc
    real(8), dimension(:), allocatable :: beta,Zk,En,Mx,My,M
    real(8), dimension(:), allocatable :: betan,Zb,Ub1,Ub2
    real(8), dimension(:), allocatable :: Mx1,My1,M1
    real(8), parameter :: eps = 1.0d-7
end module var_FS

real(8) function log1p(x)
    implicit none
    real(8), intent(in) :: x
    real(8) :: y,z

    y = 1.0d0 + x
    z = y - 1.0d0
    if (z == 0.0d0) then
        log1p = x
    else
        log1p = x*dlog(y)/z
    end if
    return
end function log1p

real(8) function soma_log(l1,l2)
    implicit none
    real(8), intent(in) :: l1,l2
    real(8) :: log1p

    if (l1 >= l2) then
        soma_log = l1 + log1p(dexp(l2-l1))
    else
        soma_log = l2 + log1p(dexp(l1-l2))
    end if
    return
end function soma_log

subroutine inicia
    use var_FS, only : Ns,N_temp,N_mc,beta,En,Zk,Mx,My,M
    implicit none
    integer :: i,j,k

    open(1,file='energia.dat')
    read(1,*) Ns,N_mc,N_temp
    N_temp = N_temp - 0
    allocate(beta(N_temp))
    allocate(En(N_temp*N_mc),Mx(N_temp*N_mc),My(N_temp*N_mc),M(N_temp*N_mc))
    allocate(Zk(N_temp))
    k = 0
    Mx = 0.0d0
    My = 0.0d0
    M = 0.0d0
    En = 0.0d0
    do i = 1,N_temp
        read(1,*) beta(i)
        do j = 1,N_mc
            k = k + 1
            read(1,*) En(k),Mx(k),My(k)
            M(k) = sqrt(Mx(k)**2 + My(k)**2)
        end do
    end do
    close(1)
    beta = 1.0d0/beta
    Zk = 1.0d0
    return
end subroutine inicia

subroutine Zk_calc
    use var_FS, only: N_temp,N_mc,beta,Zk,En,eps
    implicit none
    integer :: i,j,k,s,ie
    real(8), dimension(:), allocatable :: Zj
    real(8) :: diff,logden,lognum,Aden,soma_log

    diff = 1000.0d0

    allocate(Zj(N_temp))
    do while(diff > eps)
        do k = 1,N_temp
            Zj(k) = -10000.0d0
            ie = 0
            do i = 1,N_temp
                do s = 1,N_mc
                    ie = ie + 1
                    logden = -10000.0d0
                    lognum = -10000.0d0
                    do j = 1,N_temp
                        Aden = (beta(k)-beta(j))*En(ie)-Zk(j)+dlog(real(N_mc,8))
                        logden = soma_log(logden,Aden)
                        lognum = soma_log(lognum,0.0d0)
                    end do
                    Zj(k) = soma_log(Zj(k),lognum-logden)
                end do
            end do
        end do
        Zj = Zj - 0.5d0*(maxval(Zj)+minval(Zj))
        diff = 0.0d0
        do k = 1,N_temp
            diff = diff + abs(1.0d0 - Zj(k)/Zk(k))
        end do
        Zk = Zj
    end do
    deallocate(Zj)

    return
end subroutine Zk_calc

subroutine FS_interpol
    use var_FS
    implicit none
    integer :: N_tn
    integer :: it,i,j,s,ie
    real(8) :: Ti,Tf,dT,T,Emin,Mxmin,Mymin
    real(8) :: Zb_den,Zb_num,U1_num,U2_num
    real(8) :: Mx_num,My_num,M_num
    real(8) :: soma_log,Zd
    real(8) :: E,E2,Cv,Mxx,Myy,MM

    Emin = minval(En)
    Mxmin = minval(Mx)
    Mymin = minval(My)

    Ti = 1.0d0/beta(N_temp)
    Tf = 1.0d0/beta(1)
    write(*,'("Ti:",1x,f20.5)') Ti
    write(*,'("Tf:",1x,f20.5)') Tf
    print*, 'Qual o espaçamento dT:'
    read(*,*) dT
    N_tn = ceiling((Tf-Ti)/dT)+1

    allocate(betan(N_tn),Zb(N_tn),Ub1(N_tn),Ub2(N_tn))
    allocate(Mx1(N_tn),My1(N_tn),M1(N_tn))

    open(2,file='saida_hist.dat')
    write(2,*) ' Temp ',' En ',' Cv ',' M ',' Mx ',' My '

    do it = 1,N_tn
        T = Ti + real((it-1),8)*dT
        betan(it) = 1.0d0/T
        Zb(it) = -1000.0d0
        Ub1(it) = -1000.0d0
        Ub2(it) = -1000.0d0
        Mx1(it) = -1000.0d0
        My1(it) = -1000.0d0
        M1(it) = -1000.0d0
        ie = 0
        do i = 1,N_temp
            do s = 1,N_mc
                ie = ie + 1
                Zb_den = -1000.0d0
                Zb_num = -1000.0d0
                U1_num = -1000.0d0
                U2_num = -1000.0d0
                Mx_num = -1000.0d0
                My_num = -1000.0d0
                M_num = -1000.0d0
                do j = 1,N_temp
                    Zd = (betan(it)-beta(j))*En(ie)-Zk(j)+dlog(real(N_mc,8))
                    Zb_den = soma_log(Zb_den,Zd)
                    Zb_num = soma_log(Zb_num,0.0d0)
                    U1_num = soma_log(U1_num,dlog(En(ie)-Emin))
                    U2_num = soma_log(U2_num,dlog(En(ie)*En(ie)))
                    Mx_num = soma_log(Mx_num,dlog(Mx(ie)-Mxmin))
                    My_num = soma_log(My_num,dlog(My(ie)-Mymin))
                    M_num = soma_log(M_num,dlog(M(ie)))
                end do
                Zb(it) = soma_log(Zb(it),Zb_num-Zb_den)
                Ub1(it) = soma_log(Ub1(it),U1_num-Zb_den)
                Ub2(it) = soma_log(Ub2(it),U2_num-Zb_den)
                Mx1(it) = soma_log(Mx1(it),Mx_num-Zb_den)
                My1(it) = soma_log(My1(it),My_num-Zb_den)
                M1(it) = soma_log(M1(it),M_num-Zb_den)
            end do
        end do
        E = (dexp(Ub1(it)-Zb(it))+Emin)/real(Ns,8)
        E2 = (dexp(Ub2(it)-Zb(it)))/real(Ns*Ns,8)
        Cv = (betan(it)**2)*(E2 - E*E)*real(Ns,8)
        Mxx = (dexp(Mx1(it)-Zb(it))+Mxmin)/real(Ns)
        Myy = (dexp(My1(it)-Zb(it))+Mymin)/real(Ns)
        MM = (dexp(M1(it)-Zb(it)))/real(Ns)
        write(2,'(6(1x,f20.6))') T,E,Cv,MM,Mxx,Myy
    end do
    close(2)

    deallocate(beta,Zk,En)
    deallocate(Mx,My,M)
    deallocate(betan,Zb,Ub1,Ub2)
    deallocate(Mx1,My1,M1)

end subroutine FS_interpol

subroutine momentos
    use var_FS, only: N_temp,N_mc,beta,En,M!,Mx,My,M
    implicit none
    integer :: i,j,k
    real(8), dimension(:), allocatable :: d1,d4!,d2,d3,d4
    real(8) :: E,errE,Cvmed,errCv,bt
    real(8) :: M1,errM,Sucmed,errSuc

    allocate(d1(N_mc),d4(N_mc))!,d2(N_mc),d3(N_mc),d4(N_mc))

    k = 0
    open(3,file='Energia.dat')
    open(4,file='Magnetizacao.dat')
    write(3,*) ' Temp ',' En ',' dE ',' Cv ', ' dCv '
    write(4,*) ' Temp ',' M ',' dM ',' Suc ', ' dSuc '
    do i = 1,N_temp
        bt = beta(i)
        do j = 1,N_mc
            k = k + 1
            d1(j) = En(k)
            !d2(j) = Mx(k)
            !d3(j) = My(k)
            d4(j) = M(k)
        end do
        call jackknife1(d1,E,errE,Cvmed,errCv,bt)

        call jackknife2(d4,M1,errM,Sucmed,errSuc,bt)
        write(3,*) 1.0d0/bt,E,errE,Cvmed,errCv
        write(4,*) 1.0d0/bt,M1,errM,Sucmed,errSuc
    end do
    close(3)

    return

end subroutine momentos

subroutine jackknife1(d1,E,errE,Cvmed,errCv,bt)
    use var_FS, only: Ns,N_mc
    implicit none
    real(8), dimension(N_mc), intent(in) :: d1
    real(8), intent(in) :: bt
    real(8), intent(out) :: E,errE,Cvmed,errCv
    integer :: i,j,k
    real(8), dimension(N_mc-1) :: dn
    real(8) :: E1,E2,Cvi,dE

    call avevar(N_mc,d1,E,dE)
    E2 = sum(d1*d1)/real(N_mc,8)
    Cvmed = bt*bt*(E2 - E**2)/real(Ns,8)

    errE = 0.0d0
    errCv = 0.0d0
    do i = 1,N_mc
        k = 0
        dn = 0.0d0
        do j = 1,i-1
            k = k + 1
            dn(k) = d1(j)
        end do
        do j = i+1,N_mc
            k = k + 1
            dn(k) = d1(j)
        end do
        E1 = sum(dn)/real(N_mc-1,8)
        E2 = sum(dn**2)/real(N_mc-1,8)

        Cvi = bt*bt*(E2 - E1**2)/real(Ns,8)
        errCv = errCv + (Cvi - Cvmed)**2
        errE = errE + (E1 - E)**2

    end do
    E = E/real(Ns,8)
    errCv = sqrt(errCv)
    errE = sqrt(errE)/real(Ns,8)
    return
end subroutine jackknife1

subroutine jackknife2(d4,M,errM,Sucmed,errSuc,bt)
    use var_FS, only: Ns,N_mc
    implicit none
    real(8), dimension(N_mc), intent(in) :: d4
    real(8), intent(in) :: bt
    real(8), intent(out) :: M,errM,Sucmed,errSuc
    integer :: i,j,k
    real(8), dimension(N_mc-1) :: dn
    real(8) :: M1,M2,Sui,dM

    call avevar(N_mc,d4,M,dM)
    M2 = sum(d4*d4)/real(N_mc,8)
    Sucmed = bt*(M2 - M**2)/real(Ns,8)

    errM = 0.0d0
    errSuc = 0.0d0
    do i = 1,N_mc
        k = 0
        dn = 0.0d0
        do j = 1,i-1
            k = k + 1
            dn(k) = d4(j)
        end do
        do j = i+1,N_mc
            k = k + 1
            dn(k) = d4(j)
        end do
        M1 = sum(dn)/real(N_mc-1,8)
        M2 = sum(dn**2)/real(N_mc-1,8)

        Sui = bt*(M2 - M1**2)/real(Ns,8)
        errSuc = errSuc + (Sui - Sucmed)**2
        errM = errM + (M1 - M)**2

    end do
    M = M/real(Ns,8)
    errSuc = sqrt(errSuc)
    errM = sqrt(errM)/real(Ns,8)
    return
end subroutine jackknife2

! subroutine bootstrap
!     use ranutil
!     use var_FS
!     implicit none
!     integer, parameter :: N_BS = 500
!     integer :: i,j,k,k1,k2,o
!     real(8) :: E,errE,Cvmed,errCv,bt
!     real(8), dimension(:), allocatable :: BS_sample,E1,E2,Cv
    
!     call initrandom()
!     allocate(BS_sample(N_mc),E1(N_BS),E2(N_BS))
    
!     open(4,file='bootstrap.dat')
!     write(4,*) ' Temp ',' En ',' dE ',' Cv ', ' dCv '

!     do i = 1,N_temp
!         bt = beta(i)
!         k1 = (i-1)*N_mc
!         k2 = i*N_mc
!         do j = 1,N_BS
!             do o = 1,N_mc
!                 k = rand_int(k1,k2)
!                 BS_sample(o) = En(k)
!             end do
!             E1(j) = sum(BS_sample)/real(N_mc,8)
!             E2(j) = sum(BS_sample)/real(N_mc,8)
!         end do
!         do k = 



!         write(4,*) 1.0d0/bt,E,errE,Cvmed,errCv
!     end do


! end subroutine bootstrap

subroutine avevar(N,d1,ave,var)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N), intent(in) :: d1
    real(8), intent(out) :: ave,var
    real(8), dimension(N) :: s

    ave = sum(d1)/real(N,8)
    s = d1 - ave
    var = dot_product(s,s)
    var = sqrt((var - (sum(s)**2)/real(N,8))/real(N-1,8))
    return

end subroutine avevar

program Ferrenberg_Swendsen
    call inicia
    !call Zk_calc
    !call FS_interpol
    call momentos
end program Ferrenberg_Swendsen
        
