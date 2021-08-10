module var_inicial
    integer :: rede,contorno
    integer :: nx,ny,Ns,N_viz     !! ny deve ser par !!
    integer :: N_k,N_T,N_svt,N_svk
    integer, dimension(:), allocatable :: Nviz,jviz,cor
    real(8) :: a,Lx,Ly,D,rc,theta
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    real(8), dimension(:), allocatable :: rx,ry,mx,my,Aij
    real(8), dimension(:), allocatable :: xt,yt,xk,yk
    character(200) :: dir1,dir2
end module var_inicial

subroutine parametros_iniciais
    use var_inicial, only : rede,contorno,nx,ny,a,theta,rc
    implicit none

    print*, "Parametros iniciais para a rede:"
    print*, "--------------------------------"
    print*, "Qual a rede a ser simulada?     "
    print*, "Opções: kagome = 1, triangular = 2 ou quadrada = 3"
    read(*,*) rede
    print*, "--------------------------------"
    print*, "Entre com o tamanho da rede: nx e ny"
    read(*,*) nx,ny
    print*, "Entre com o espaçamento de rede:"
    read(*,*) a
    print*, "Entre com o raio de corte (indicado = 6.0d0):"
    read(*,*) rc
    print*, "--------------------------------"
    print*, "A simulação é com ou sem contorno?"
    print*, "Com contorno = 1; Sem contorno = 2"
    read(*,*) contorno
    print*, "Ângulo de rotação do spin:"
    read(*,*) theta
  
    rc = rc*a
    return
end subroutine parametros_iniciais

subroutine unit_cel
    use var_inicial, only : rede,nx,ny,contorno,dir1,theta
    implicit none
    logical :: direx
    character(20) :: tipo

    write(dir1,"('Inicial/')")
    inquire(file=dir1,exist=direx)
    if (direx .eqv. .false.) then
        call system("mkdir " // trim(dir1))
    end if

    select case (rede)
        case (1)
            tipo = 'Kagome'
            if (mod(ny,2) .ne. 0) ny = ny + 1
            call diretorios(tipo,nx,ny,contorno,theta)
            call kagome
        case (2)
            tipo = 'Triang'
            if (mod(ny,2) .ne. 0) ny = ny + 1
            call diretorios(tipo,nx,ny,contorno,theta)
            call triangular
        case (3)
            tipo = 'Square'
            call diretorios(tipo,nx,ny,contorno,theta)
            call quadrada
    end select

    if (contorno == 1) then
        call Aij_com_contorno_x
        print*, "Rede com condições de contorno gerada."
    else if (contorno == 2) then
        call Aij_sem_contorno_x
        print*, "Rede sem condições de contorno gerada."
    end if

end subroutine unit_cel

subroutine diretorios(tipo,nx,ny,cont,theta)
    use var_inicial, only : dir2
    implicit none
    real(8), intent(in) :: theta
    character(8),intent(in) :: tipo
    integer,intent(in) :: nx,ny,cont
    character(10) :: BC
    logical :: direx

    if (cont == 1) then
        BC = "BC"
    else
        BC = "NBC"
    end if

    write(dir2,"('Inicial/',A,'_',I2.2,'x',I2.2,'_',A,'_',I4.4,'/')") trim(tipo),nx,ny,trim(BC),int(100.0d0*theta)

    inquire(file=dir2,exist=direx)
    if (direx .eqv. .false.) then
        call system("mkdir " // trim(dir2))
    end if
    return
end subroutine diretorios

subroutine kagome
    use var_inicial, only : a,pi,nx,ny,Ns,rx,ry,mx,my,Lx,Ly,cor,dir2,theta
    implicit none
    integer :: i,j,k,is
    real(8), dimension(:), allocatable :: x0,y0,ex,ey

    !!Cria a célula unitária de spins na configuração Kagome!!
    Ns = 3*nx*ny
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns),cor(Ns))
    allocate(x0(3),y0(3),ex(3),ey(3))

    x0(1) = 0.0d0
    y0(1) = 0.0d0

    x0(2) = a*cos(60.0d0*pi/180.0d0)
    y0(2) = a*sin(60.0d0*pi/180.0d0)

    x0(3) = a
    y0(3) = 0.0d0

    ex(1) = -cos((30.0d0 + theta)*pi/180.0d0)
    ey(1) = -sin((30.0d0 + theta)*pi/180.0d0)

    ex(2) = -cos((90.0d0 + theta)*pi/180.0d0) 
    ey(2) = -sin((90.0d0 + theta)*pi/180.0d0)

    ex(3) = -cos((150.0d0 + theta)*pi/180.0d0)
    ey(3) = -sin((150.0d0 + theta)*pi/180.0d0)

    is = 0
    do j = 1,ny
        do i = 1,nx
            do k = 1,3
                is = is + 1
                rx(is) = x0(k) + 2.0d0*a*(i-1) + a*mod(j+1,2)
                ry(is) = y0(k) + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
                mx(is) = ex(k)
                my(is) = ey(k)
                cor(is) = k
            end do
        end do
    end do

    Lx = maxval(rx) - minval(rx)
    Ly = maxval(ry) - minval(ry) + a*sin(60.0d0*pi/180.0d0)

    rx = rx - 0.5d0*Lx
    ry = ry - 0.5d0*Ly

    open(50,file=trim(dir2) // 'cor.dat')
    write(50,*) 3
    write(50,*) 30,210
    write(50,*) 90,270
    write(50,*) 150,330
    close(50)

    return
end subroutine kagome

subroutine triangular
    use var_inicial, only : a,pi,nx,ny,Ns,rx,ry,mx,my,Lx,Ly,cor,dir2,theta
    implicit none
    integer :: i,j,k,is
    real(8), dimension(:), allocatable :: x0,y0,ex,ey

    Ns = 3*nx*ny
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns),cor(Ns))
    allocate(x0(3),y0(3),ex(3),ey(3))

    x0(1) = 0.0d0
    y0(1) = 0.0d0

    x0(2) = a*cos(60.0d0*pi/180.0d0)
    y0(2) = a*sin(60.0d0*pi/180.0d0)

    x0(3) = a
    y0(3) = 0.0d0

    ex(1) = -cos((120.0d0 + theta)*pi/180.0d0)
    ey(1) = -sin((120.0d0 + theta)*pi/180.0d0)

    ex(2) = -cos((180.0d0 + theta)*pi/180.0d0) 
    ey(2) = -sin((180.0d0 + theta)*pi/180.0d0)

    ex(3) = -cos((240.0d0 + theta)*pi/180.0d0)
    ey(3) = -sin((240.0d0 + theta)*pi/180.0d0)

    is = 0
    do j = 1,ny
        do i = 1,nx
            do k = 1,3
                is = is + 1
                rx(is) = x0(k) + 2.0d0*a*(i-1) + a*mod(j+1,2)
                ry(is) = y0(k) + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
                mx(is) = ex(k)
                my(is) = ey(k)
                cor(is) = k
            end do
        end do
    end do

    Lx = maxval(rx) - minval(rx)
    Ly = maxval(ry) - minval(ry) + a*sin(60.0d0*pi/180.0d0)

    rx = rx - 0.5d0*Lx
    ry = ry - 0.5d0*Ly

    open(50,file=trim(dir2) // 'cor.dat')
    write(50,*) 3
    write(50,*) 120,300
    write(50,*) 180,0
    write(50,*) 240,60
    close(50)

    return
end subroutine triangular

subroutine quadrada
    use var_inicial, only : a,nx,ny,Ns,rx,ry,mx,my,Lx,Ly,cor,dir2,theta,pi
    implicit none
    integer :: i,j,k,is
    real(8), dimension(:), allocatable :: x0,y0,ex,ey

    Ns = 2*nx*ny
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns),cor(Ns))
    allocate(x0(2),y0(2),ex(2),ey(2))

    x0(1) = 0.5d0*a
    y0(1) = 0.0d0

    x0(2) = 0.0d0
    y0(2) = 0.5d0*a

    ex(1) = cos((0.0d0 + theta)*pi/180.0d0)
    ey(1) = sin((0.0d0 + theta)*pi/180.0d0)

    ex(2) = cos((90.0d0 + theta)*pi/180.0d0) 
    ey(2) = sin((90.0d0 + theta)*pi/180.0d0)

    is = 0
    do j = 1,ny
        do i = 1,nx
            do k = 1,2
                is = is + 1
                rx(is) = x0(k) + a*(i-1)
                ry(is) = y0(k) + a*(j-1)
                mx(is) = ex(k)
                my(is) = ey(k)
                cor(is) = k
            end do
        end do
    end do

    Lx = maxval(rx) - minval(rx) + 0.5d0*a
    Ly = maxval(ry) - minval(ry) + 0.5d0*a

    rx = rx - 0.5d0*Lx
    ry = ry - 0.5d0*Ly

    open(50,file=trim(dir2) // 'cor.dat')
    write(50,*) 2
    write(50,*) 0,180
    write(50,*) 90,270
    close(50)

    return
end subroutine quadrada

subroutine config
    use var_inicial!, only : Ns,rx,ry,mx,my,cor,dir2
    implicit none
    integer :: i

    open(10,file=trim(dir2) // "config0.xyz")
    write(10,*) Ns
    write(10,*) ' '
    do i = 1,Ns
        write(10,*) rx(i),ry(i),mx(i),my(i),cor(i)
    end do
    close(10)

    call flush()
    
    open(10,file=trim(dir2) // "vertice_k.xyz")
    write(10,*) N_k
    write(10,*) ' '
    do i = 1,N_k
        write(10,*) xk(i),yk(i)
    end do
    close(10)

    call flush()
    
    open(10,file=trim(dir2) // "vertice_T.xyz")
    write(10,*) N_T
    write(10,*) ' '
    do i = 1,N_T
        write(10,*) xT(i),yT(i)
    end do
    close(10)

    return
end subroutine config

subroutine Aij_sem_contorno_x
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,rc,dir2
    implicit none
    integer :: i,j
    real(8) :: x,y,dij,D1,D2,Aijmin,Aijmax,soma,Aij

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')
    N_viz = 0
    !rc = 6.0d0*a
    Aijmin = 10000000.0d0
    Aijmax = -10000000.0d0
    do i = 1,Ns
        soma = 0.0d0
        do j = 1,Ns
            if (j.ne.i) then
                x = rx(i) - rx(j)
                y = ry(i) - ry(j)
                if ((sqrt(x*x+y*y) <= rc).and.(sqrt(x*x+y*y)>0.1d0)) then
                    N_viz = N_viz + 1
                    dij = 1.0d0/sqrt(x*x + y*y)
                    x = x*dij
                    y = y*dij
                    D1 = mx(i)*mx(j) + my(i)*my(j)
                    D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                    Aij = (D1 - D2)*dij**3
                    soma = soma + Aij
                    write(11,*) j,Aij
                end if
            end if
        end do
        if (soma > Aijmax) Aijmax = soma
        if (soma < Aijmin) Aijmin = soma
        write(12,*) N_viz
    end do

    print*, "Aij minimo",Aijmin
    print*, "Aij máximo",Aijmax

    close(11)
    close(12)

    return
end subroutine Aij_sem_contorno_x

subroutine Aij_com_contorno_x
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,Lx,Ly,rc,dir2
    implicit none
    integer :: i,j,ni,nj,nc
    real(8) :: x,y,dij,D1,D2,Aijmin,Aijmax,soma,Aij

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    N_viz = 0
    nc = 3
    !rc = 6.0d0*a
    Aijmin = 10000000.0d0
    Aijmax = -10000000.0d0
    do i = 1,Ns
        soma = 0.0d0
        do ni = -nc,nc
            do nj = -nc,nc
                do j = 1,Ns
                    x = rx(i) - rx(j) + real(ni*Lx,8)
                    y = ry(i) - ry(j) + real(nj*Ly,8)
                    if ((sqrt(x*x+y*y) <= rc).and.(sqrt(x*x+y*y)>0.1d0)) then
                        N_viz = N_viz + 1
                        dij = 1.0d0/sqrt(x*x + y*y)
                        x = x*dij
                        y = y*dij
                        D1 = mx(i)*mx(j) + my(i)*my(j)
                        D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                        Aij = (D1 - D2)*dij**3
                        soma = soma + Aij
                        !Aij(N_viz) = (D1 - D2)*dij**3
                        !jviz(N_viz) = j
                        write(11,*) j,Aij
                    end if
                end do
            end do
        end do
        !Nviz(i) = N_viz
        if (soma > Aijmax) Aijmax = soma
        if (soma < Aijmin) Aijmin = soma
        write(12,*) N_viz
    end do

    print*, "Aij minimo",Aijmin
    print*, "Aij máximo",Aijmax

    close(11)
    close(12)
    return
end subroutine Aij_com_contorno_x

subroutine vertice_T
    use var_inicial, only : Ns,nx,ny,N_T,N_svt,xT,yT,a,Lx,Ly,pi,dir2,rx,ry,mx,my
    implicit none
    integer :: i,j,k,ni,nj
    real(8) :: vxu,vyu
    real(8) :: x,y,DD,Lij

    N_T = nx*ny
    N_svt = 6*N_T
    allocate(xT(N_T),yT(N_T))
    
    vxu = 1.5d0*a
    vyu = 0.5d0*a*tan(60.0d0*pi/180.0d0)

    k = 0
    do j = 1,ny
        do i = 1,nx
            k = k + 1
            xT(k) = vxu + 2.0d0*a*(i-1) + a*mod(j+1,2)
            yT(k) = vyu + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
        end do
    end do

    xT = xT - 0.5d0*Lx
    yT = yT - 0.5d0*Ly

    open(1,file=trim(dir2) // "vertice_T.dat")
    do i = 1,N_T
        do ni = -3,3
            do nj = -3,3
                do j = 1,Ns
                    x = xT(i) - rx(j) + real(ni*Lx,8)
                    y = yT(i) - ry(j) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<1.3d0*a .and. DD>0.9d0*a) then
                        Lij = x*mx(j) + y*my(j)
                        write(1,*) j,int(sign(1.d0,Lij))
                    end if
                end do
            end do
        end do
    end do

    do i = 1,Ns
        do ni = -3,3
            do nj = -3,3
                do j = 1,N_T
                    x = xT(j) - rx(i) + real(ni*Lx,8)
                    y = yT(j) - ry(i) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<1.3d0*a .and. DD>0.9d0*a) then
                        write(1,*) j
                    end if
                end do
            end do
        end do
    end do
    close(1)

    return

end subroutine vertice_T

subroutine vertice_K
    use var_inicial, only : Ns,nx,ny,N_k,N_svk,xk,yk,a,Lx,Ly,pi,dir2,rx,ry,mx,my
    implicit none
    integer :: i,j,k,l,ni,nj
    real(8) :: vxu(2),vyu(2)
    real(8) :: x,y,DD,Lij

    N_K = 2*nx*ny
    N_svk = 3*N_k
    allocate(xk(N_k),yk(N_k))
    
    vxu(1) = 0.5d0*a
    vyu(1) = 0.5d0*a*tan(30.0d0*pi/180.0d0)

    vxu(2) = 0.5d0*a
    vyu(2) = a*tan(60.0d0*pi/180.0d0) - vyu(1)

    k = 0
    do j = 1,ny
        do i = 1,nx
            do l = 1,2
                k = k + 1
                xk(k) = vxu(l) + 2.0d0*a*(i-1) + a*mod(j+1,2)
                yk(k) = vyu(l) + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
            end do
        end do
    end do

    xk = xk - 0.5d0*Lx
    yk = yk - 0.5d0*Ly

    open(1,file=trim(dir2) // "vertice_K.dat")
    do i = 1,N_K
        do ni = -3,3
            do nj = -3,3
                do j = 1,Ns
                    x = xk(i) - rx(j) + real(ni*Lx,8)
                    y = yk(i) - ry(j) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<a .and. DD>0.1d0*a) then
                        Lij = x*mx(j) + y*my(j)
                        write(1,*) j,int(sign(1.d0,Lij))
                    end if
                end do
            end do
        end do
    end do

    do i = 1,Ns
        do ni = -3,3
            do nj = -3,3
                do j = 1,N_K
                    x = xk(j) - rx(i) + real(ni*Lx,8)
                    y = yk(j) - ry(i) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<a .and. DD>0.1d0*a) then
                        write(1,*) j
                    end if
                end do
            end do
        end do
    end do
    close(1)

    return

end subroutine vertice_K

subroutine output
    use var_inicial
    implicit none
    real(8) :: theta_h
    !integer :: i

    call config

    print*, "Qual theta para histerese?"
    read(*,*) theta_h

    open(13,file=trim(dir2) // 'input.dat')
    write(13,*) Ns,"   !Ns"
    write(13,*) N_viz,"   !N_viz"
    write(13,*) 11.25d0,"   !Hmax"
    write(13,*) 1.875d0,"   !delH"
    write(13,*) 50.0d0,"   !Campo externo máximo"
    write(13,*) theta_h,"   !Theta para histerese"
    write(13,*) 50000,"    !N_mc"
    write(13,*) 20,"    !N_temp"
    write(13,*) 10.0,"    !Temperatura inicial"
    write(13,*) 0.5,"    !Temperatura final"
    write(13,*) N_svt,N_K,N_T
    close(13)
    if (N_svt .ne. N_svk) then
        print*, 'Erro em N_K e N_T',N_K,N_T
    end if
end subroutine output

program inicial_v2

    call parametros_iniciais
    call unit_cel
    call vertice_k
    call vertice_T
    call output
    
end program inicial_v2