!   Código para gerar o Fator de Estrutura Magnético de um spin ice
!
! Estrutura do arquivo de input:
! Ns,Nmssf --> Número de spin na rede e o número total de configurações para média
! rx(1),ry(1),mx(1),my(1)
!...
! rx(Ns),ry(Ns),mx(Ns),my(N)
!
!1: S(1) .... S(Ns)
!2: S(1) .... S(Ns)
!...
!Nmssf: S(1) .... S(Ns)
!---------------------------------------------------------!
module var
    integer :: Ns,Nmssf
    integer :: Nq
    integer, dimension(:), allocatable :: S
    real(8), parameter :: pi = 4.0d0*atan(1.0d0) 
    real(8), dimension(:), allocatable :: rx,ry,mx,my
    real(8), dimension(:), allocatable :: qx,qy,qmod
    real(8), dimension(:), allocatable :: Ax,Ay,Bx,By
    real(8), dimension(:), allocatable :: Intensidade
    real(8) :: qmax
end module var

subroutine ler_input
    use var, only : Ns,Nmssf,rx,ry,mx,my,S,Nq,qmax
    implicit none
    integer :: i

    open(10,file='mssf.dat')
    
    read(10,*) Ns, Nmssf
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns),S(Ns))
    S = 0
    !print*, 'Entre com N_sample'
    !read(*,*) Nmssf

    do i = 1,Ns
        read(10,*) rx(i),ry(i),mx(i),my(i)
    end do

    print*, "Entre com o número de pontos pra dq"
    read(*,*) Nq
    print*, "Entre com o módulo máximo de q"
    read(*,*) qmax

    return
end subroutine ler_input

subroutine inicia_q
    use var, only : Nq,qx,qy,qmod,pi,qmax
    implicit none
    integer :: i,j,iq
    real(8) :: dq

    allocate(qx(Nq*Nq),qy(Nq*Nq),qmod(Nq*Nq))

    dq = 2.0d0*qmax*pi/real(Nq-1,8)
    iq = 0

    do i = 1,Nq
        do j = 1,Nq
            iq = iq + 1
            qx(iq) = -pi*qmax + (i-1)*dq
            qy(iq) = -pi*qmax + (j-1)*dq
            if (qx(iq)==0.0d0 .and. qy(iq)==0.0d0) then
                qmod(iq) = 1.0d0
            else
                qmod(iq) = 1.0d0/(qx(iq)**2 + qy(iq)**2)
            end if
        end do
    end do

    return
end subroutine inicia_q

subroutine inicia_A
    use var, only : Nq,Ns,rx,ry,mx,my,qx,qy,qmod,Ax,Ay,Bx,By,Intensidade
    implicit none
    integer :: i,j,k,l,iq
    real(8) :: spx,spy,a,theta

    allocate(Ax(Nq*Nq*Ns),Ay(Nq*Nq*Ns),Bx(Nq*Nq*Ns),By(Nq*Nq*Ns))
    allocate(Intensidade(Nq*Nq))
    Intensidade = 0.0d0
    
    l = 0
    iq = 0
    do j = 1,Nq
        do k = 1,Nq
            iq = iq + 1
            do i = 1,Ns
                l = l + 1
                a = qx(iq)*mx(i) + qy(iq)*my(i)

                theta = qx(iq)*rx(i) + qy(iq)*ry(i)

                spx = mx(i) - a*qx(iq)*qmod(iq)
                spy = my(i) - a*qy(iq)*qmod(iq)

                Ax(l) = spx*cos(theta)
                Ay(l) = spy*cos(theta)

                Bx(l) = spx*sin(theta)
                By(l) = spy*sin(theta)

            end do
        end do
    end do
    deallocate(qmod)

    return
end subroutine inicia_A

subroutine inicia_var
    implicit none

    call ler_input 
    call inicia_q       ! Vetores de onda
    call inicia_A   

    return
end subroutine inicia_var

subroutine ler_S(soma)
    use var, only : Ns,S
    implicit none
    real(8), dimension(Ns), intent(inout) :: soma
    integer :: i,si

    S = 0
    do i = 1,Ns
        read(10,*) si
        S(i) = 2*si - 1
        soma(i) = soma(i) + 1.0d0*S(i)
    end do

    return
end subroutine ler_S

subroutine intensity(soma)
    use var, only : Ns,Nq,Ax,Ay,Bx,By,S
    implicit none
    real(8), dimension(Nq*Nq), intent(out) :: soma
    integer :: i,j,k,l,iq
    real(8) :: iAx,iAy,iBx,iBy

    soma = 0.0d0

    iq = 0
    l = 0

    do j = 1,Nq
        do k = 1,Nq
            iq = iq + 1
            iAx = 0.0d0
            iAy = 0.0d0
            iBx = 0.0d0
            iBy = 0.0d0
            do i = 1,Ns
                l = l + 1
                iAx = iAx + S(i)*Ax(l)
                iAy = iAy + S(i)*Ay(l)
                iBx = iBx + S(i)*Bx(l)
                iBy = iBy + S(i)*By(l)
            end do
            soma(iq) = (iAx**2 + iAy**2) + (iBx**2 + iBy**2)
        end do
    end do

    !soma = soma / real(Ns,8)

    return
end subroutine intensity

subroutine output(k,soma)
    use var, only : Nq,Ns,qx,qy,rx,ry,mx,my,Intensidade
    implicit none
    integer, intent(in) :: k
    real(8), dimension(Ns), intent(in) :: soma
    integer :: iq,i
    !real(8) :: Imax
    character(60) :: n1,n2

    !Imax = maxval(Intensidade)
    !Intensidade = Intensidade / Imax

    write(n1,"('mssf',I3.3,'.dat')") k
    write(n2,"('mssf',I3.3,'.xyz')") k

    open(20,file=trim(n1))
    open(30,file=trim(n2))

    do iq = 1,Nq*Nq
        write(20,*) qx(iq),qy(iq),Intensidade(iq)
    end do

    write(30,*) Ns
    write(30,*) ' '
    do i = 1,Ns
        write(30,*) rx(i),ry(i),soma(i)*mx(i),soma(i)*my(i)
    end do
    
    close(20)
    close(30)

    return
end subroutine output

program main
    use var, only : Nmssf,Intensidade,Nq,Ns
    implicit none
    integer :: k
    real(8), dimension(:), allocatable :: soma,somaS

    call inicia_var
    allocate(soma(Nq*Nq))
    allocate(somaS(Ns))
    
    somaS = 0.0d0
    Intensidade = 0.0d0
    do k = 1,Nmssf
        call ler_S(somaS)
        call intensity(soma)
        Intensidade = Intensidade + soma
        !call output(k)
    end do

    Intensidade = Intensidade / real(Ns*Nmssf,8)
    somaS = somaS / real(Nmssf,8)
    deallocate(soma)
    call output(1,somaS)
    deallocate(somaS)

end program main
