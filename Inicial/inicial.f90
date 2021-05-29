module var_inicial
    integer :: rede,contorno
    integer :: nx,ny,Ns,N_viz     !! ny deve ser par !!
    integer, dimension(:), allocatable :: Nviz,jviz,cor
    real(8) :: a,Lx,Ly,D,rc,theta
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    real(8), dimension(:), allocatable :: rx,ry,mx,my,Aij
    character(200) :: dir1,dir2
end module var_inicial

subroutine parametros_iniciais
    use var_inicial, only : rede,contorno,nx,ny,a,theta
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
    print*, "--------------------------------"
    print*, "A simulação é com ou sem contorno?"
    print*, "Com contorno = 1; Sem contorno = 2; Ising = 3"
    read(*,*) contorno
    print*, "Ângulo de rotação do spin:"
    read(*,*) theta
  
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
    else if (contorno == 3) then
        call Aij_com_contorno_Ising
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

    !if (theta < 10.0d0) then
    !    write(dir2,"('Inicial/',A,'_',I2.2,'x',I2.2,'_',A,'_',F2.0,'/')") trim(tipo),nx,ny,trim(BC),theta
    !else
    !    write(dir2,"('Inicial/',A,'_',I2.2,'x',I2.2,'_',A,'_',F3.0,'/')") trim(tipo),nx,ny,trim(BC),theta
    !end if

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
    ey(2) = sin((90.0d0 + theta)*pi/180.0d0)

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
    use var_inicial, only : Ns,rx,ry,mx,my,cor,dir2
    implicit none
    integer :: i

    open(10,file=trim(dir2) // "config0.xyz")
    write(10,*) Ns
    write(10,*) ' '
    do i = 1,Ns
        write(10,*) rx(i),ry(i),mx(i),my(i),cor(i)
    end do
    close(10)

    return
end subroutine config

subroutine Aij_sem_contorno
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,Aij,Nviz,D,jviz
    implicit none
    integer :: i,j
    real(8) :: x,y,dij,D1,D2
    real(8) :: minA,maxA

    allocate(Aij(Ns*Ns),Nviz(0:Ns),jviz(Ns*Ns))
    Aij(:) = 0.0d0
    Nviz(:) = 0
    jviz(:) = 0
    N_viz = 0
    do i = 1,Ns
        do j = 1,Ns
            if (j.ne.i) then
                N_viz = N_viz + 1
                x = rx(i) - rx(j)
                y = ry(i) - ry(j)
                dij = 1.0d0/sqrt(x*x + y*y)
                x = x*dij
                y = y*dij
                D1 = mx(i)*mx(j) + my(i)*my(j)
                D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                Aij(N_viz) = (D1 - D2)*dij**3
                jviz(N_viz) = j
            end if
        end do
        Nviz(i) = N_viz
    end do

    minA = minval(abs(Aij))
    maxA = maxval(abs(Aij))
    D = 1.0d0/maxA
    Aij = Aij*D

    return
end subroutine Aij_sem_contorno

subroutine Aij_com_contorno
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,Aij,Nviz,D,Lx,Ly,rc,a,jviz
    implicit none
    integer :: i,j,ni,nj,nc
    real(8) :: x,y,dij,D1,D2
    real(8) :: minA,maxA

    allocate(Aij(Ns**3),Nviz(0:Ns),jviz(Ns**3))

    Aij(:) = 0.0d0
    Nviz(:) = 0
    jviz(:) = 0
    N_viz = 0
    nc = 1
    rc = 5.0d0*a
    do i = 1,Ns
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
                        Aij(N_viz) = (D1 - D2)*dij**3
                        jviz(N_viz) = j
                    end if
                end do
            end do
        end do
        Nviz(i) = N_viz
    end do

    minA = 10000.0d0
    do i = 1,N_viz
        if (abs(Aij(i)) .lt. minA) then
            minA = abs(Aij(i))
        end if
    end do
    maxA = maxval(abs(Aij))
    D = 1.0d0/maxA
    Aij = Aij*D

    return
end subroutine Aij_com_contorno

subroutine Aij_sem_contorno_x
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,rc,a,dir2
    implicit none
    integer :: i,j
    real(8) :: x,y,dij,D1,D2,Aijmin,Aijmax,soma,Aij

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')
    N_viz = 0
    rc = 200000.0d0*a
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
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,Lx,Ly,rc,a,dir2
    implicit none
    integer :: i,j,ni,nj,nc
    real(8) :: x,y,dij,D1,D2,Aijmin,Aijmax,soma,Aij

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    N_viz = 0
    nc = 3
    rc = 200.0d0*a
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

subroutine Aij_com_contorno_Ising
    use var_inicial, only : Ns,N_viz,rx,ry,Lx,Ly,rc,dir2
    implicit none
    integer :: i,j,ni,nj,nc
    real(8) :: x,y,Aij

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    N_viz = 0
    nc = 0
    rc = 1.d0
    do i = 1,Ns
        do ni = -nc,nc
            do nj = -nc,nc
                do j = 1,Ns
                    x = rx(i) - rx(j) + real(ni*Lx,8)
                    y = ry(i) - ry(j) + real(nj*Ly,8)
                    if ((sqrt(x*x+y*y) <= rc).and.(sqrt(x*x+y*y)>0.1d0)) then
                        N_viz = N_viz + 1
                        Aij = 1.0d0
                        write(11,*) j,Aij
                    end if
                end do
            end do
        end do
        write(12,*) N_viz
    end do


    close(11)
    close(12)
    return
end subroutine Aij_com_contorno_Ising

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
    close(13)

end subroutine output

subroutine campo
    use var_inicial, only : Ns,rx,ry,mx,my
    implicit none
    integer,parameter :: Np = 400,seed = 83434
    integer :: i,j,k
    real(8) :: dx,dy,x,y,r,Bx,By,xx,yy,r3
    real(8) :: xmin,xmax,ymin,ymax,S(Ns)

    call srand(seed)

    do i = 1,Ns
        if (rand() < 0.3d0) then
            S(i) = -1.0d0
        else
            S(i) = 1.0d0
        end if
    end do

    xmin = minval(rx)
    xmax = maxval(rx)
    ymin = minval(ry)
    ymax = maxval(ry)

    dx = 150.0d0/real(Np-1,8)
    dy = dx

    open(1000,file='campo.dat')

    do i = 1,Np
        x = -75.0d0 + (i-1)*dx
        do j = 1,Np
            y = -75.0d0 + (j-1)*dy
            Bx = 0.0d0
            By = 0.0d0
            if ((x > xmin .and. x < xmax) .and. (y > ymin .and. y < ymax)) then
                Bx = 0.0d0
                By = 0.0d0
            else
                do k = 1,Ns
                    xx = x - rx(k)
                    yy = y - ry(k)
                    r = sqrt(xx**2 + yy**2)
                    xx = xx/r
                    yy = yy/r
                    r3 = r**3
                    Bx = Bx + S(k)*(3.0d0*(mx(k)*xx + my(k)*yy)*xx - mx(k))/r3
                    By = By + S(k)*(3.0d0*(mx(k)*xx + my(k)*yy)*yy - my(k))/r3
                end do
            end if
            write(1000,*) x,y,Bx,By
        end do
    end do
    
    close(1000)

end subroutine campo

program inicial_v2

    call parametros_iniciais
    call unit_cel
    call output
    call campo
    
end program inicial_v2