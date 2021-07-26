module var_FS
    integer :: Ns,N_temp,N_mc
    real(8), dimension(:), allocatable :: beta,Zk,En
    real(8), dimension(:), allocatable :: betan,Zb,Ub1,Ub2
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
    use var_FS, only : Ns,N_temp,N_mc,beta,En,Zk
    implicit none
    integer :: i,j,k

    open(1,file='energia.dat')
    read(1,*) Ns,N_mc,N_temp
    !N_temp = N_temp - 1
    allocate(beta(N_temp))
    allocate(En(N_temp*N_mc))
    allocate(Zk(N_temp))
    k = 0
    do i = 1,N_temp
        read(1,*) beta(i)
        do j = 1,N_mc
            k = k + 1
            read(1,*) En(k)
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
    real(8) :: Ti,Tf,dT,T,Emin
    real(8) :: Zb_den,Zb_num,U1_num,U2_num
    real(8) :: soma_log,Zd
    real(8) :: E,E2,Cv
    real(8), dimension(:), allocatable :: Eshift

    allocate(Eshift(size(En)))
    Emin = minval(En)
    Eshift = En - Emin
    Ti = 1.0d0/beta(N_temp)
    Tf = 1.0d0/beta(1)
    write(*,'("Ti:",1x,f20.5)') Ti
    write(*,'("Tf:",1x,f20.5)') Tf
    print*, 'Qual o espaçamento dT:'
    read(*,*) dT
    N_tn = ceiling((Tf-Ti)/dT)+1
    allocate(betan(N_tn),Zb(N_tn),Ub1(N_tn),Ub2(N_tn))

    open(2,file='saida_hist.dat')
    write(2,*) 'Temp','En','En2','Cv','Zb'

    do it = 1,N_tn
        T = Ti + real((it-1),8)*dT
        betan(it) = 1.0d0/T
        Zb(it) = -10000.0d0
        Ub1(it) = -10000.0d0
        Ub2(it) = -10000.0d0
        ie = 0
        do i = 1,N_temp
            do s = 1,N_mc
                ie = ie + 1
                Zb_den = -10000.0d0
                Zb_num = -10000.0d0
                U1_num = -10000.0d0
                U2_num = -10000.0d0
                do j = 1,N_temp
                    Zd = (betan(it)-beta(j))*En(ie)-Zk(j)+dlog(real(N_mc,8))
                    Zb_den = soma_log(Zb_den,Zd)
                    Zb_num = soma_log(Zb_num,0.0d0)
                    U1_num = soma_log(U1_num,Eshift(ie))
                    U2_num = soma_log(U2_num,En(ie**2))
                end do
                Zb(it) = soma_log(Zb(it),Zb_num-Zb_den)
                Ub1(it) = soma_log(Ub1(it),U1_num-Zb_den)
                Ub2(it) = soma_log(Ub2(it),U2_num-Zb_den)
            end do
        end do
        E = (dexp(Ub1(it)-Zb(it))+Emin)/real(Ns,8)
        E2 = (dexp(Ub2(it)-Zb(it)))/real(Ns**2,8)
        Cv = (T*T)*(E2 - E*E)*real(Ns,8)
        write(2,'(5(1x,f20.6))') T,E,E2,Cv,Zb(it)
    end do
    close(2)

    deallocate(Eshift)
    deallocate(beta,Zk,En)
    deallocate(betan,Zb,Ub1,Ub2)

end subroutine FS_interpol

        