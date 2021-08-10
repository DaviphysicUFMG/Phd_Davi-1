
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


module var_annealing
   real(8), parameter :: pi = 4.0d0*atan(1.0d0)
   integer :: seed
   integer :: Ns,N_viz,N_mc,N_temp,N_skt,N_K,N_T
   integer, dimension(:), allocatable :: S,Nviz,jviz
   integer, dimension(:), allocatable :: Sk,Vk,Rk,St,Vt,Rt
   real(8), dimension(:), allocatable :: Aij,Bi
   real(8), dimension(:), allocatable :: rx,ry,mx,my
   real(8) :: Ti,Tf,dT
   real(8) :: E_tot,Temp,beta
   character(600) :: dir1,dir2

   integer, parameter :: iunit_conf = 10, iunit_en = 20
end module var_annealing

subroutine inicial
   use var_annealing, only : N_temp,Ti,Tf,dT,iunit_en
   implicit none

   ! Ler arquivos de input !
   
   call ler_input(1)
   call ler_config(2)
   call ler_Aij(3)
   call ler_Nviz(4)
   call ler_Vertices(5)
   
   ! Iniciar Campos !

   call inicia_Bi

   ! Iniciar Diretórios !

   call diretorios

   ! Arquivos de OutPut !

   call En_save(iunit_en,0)

   ! Cálculos Gerais !

   dT = (Tf - Ti)/real(N_temp-1,8)

   return
end subroutine inicial

subroutine ler_input(iunit)
   !use mtmod, only : getseed,sgrnd
   use ranutil, only : initrandom
   use var_annealing, only : Ns,N_viz,N_mc,N_temp,Ti,Tf,N_skt,N_K,N_T!,seed
   implicit none
   integer, intent(in) :: iunit
   real(8) :: a

   open(unit=iunit,file='input.dat',status='old',action='read')
   read(iunit,*) Ns
   read(iunit,*) N_viz
   read(iunit,*) a
   read(iunit,*) a
   read(iunit,*) a
   read(iunit,*) a
   read(iunit,*) N_mc
   read(iunit,*) N_temp
   read(iunit,*) Ti
   read(iunit,*) Tf
   read(iunit,*) N_skt,N_K,N_T
   close(iunit)

   !seed = getseed()
   !call sgrnd(seed)
   call initrandom()

   return
end subroutine ler_input

subroutine ler_config(iunit)
   use var_annealing, only : Ns,S,rx,ry,mx,my
   implicit none
   integer, intent(in) :: iunit
   integer :: i
   real(8) :: a

   allocate(S(Ns),rx(Ns),ry(Ns),mx(Ns),my(Ns))

   open(unit=iunit,file='config0.xyz',status='old',action='read')
   read(iunit,*) i
   read(iunit,*) 
   do i = 1,Ns
      read(iunit,*) rx(i),ry(i),mx(i),my(i),a
   end do
   close(iunit)

   S = 1

   return
end subroutine ler_config

subroutine ler_Aij(iunit)
   use var_annealing, only : N_viz,Aij,jviz
   implicit none
   integer, intent(in) :: iunit
   integer :: i

   allocate(Aij(N_viz),jviz(N_viz))
   open(unit=iunit,file='Aij.dat',status='old',action='read')
   do i = 1,N_viz
      read(iunit,*) jviz(i),Aij(i)
   end do
   close(iunit)

   return
end subroutine ler_Aij

subroutine ler_Nviz(iunit)
   use var_annealing, only : Ns,Nviz
   implicit none
   integer, intent(in) :: iunit
   integer :: i

   allocate(Nviz(0:Ns))
   open(unit=iunit,file='Nviz.dat',status='old',action='read')

   Nviz(0) = 0
   do i = 1,Ns
      read(iunit,*) Nviz(i)
   end do
   close(iunit)

   return
end subroutine ler_Nviz

subroutine ler_Vertices(iunit)
   use var_annealing, only : N_skt,Sk,Vk,Rk,St,Vt,Rt
   implicit none
   integer, intent(in) :: iunit
   integer :: i

   !-----------------------------------------------------------------!

   allocate(St(N_skt),Vt(N_skt),Rt(N_skt))
   open(unit=iunit,file="vertice_T.dat",status='old',action='read')
   do i = 1,N_skt
      read(iunit,*) Vt(i),Rt(i)
   end do
   do i = 1,N_skt
      read(iunit,*) St(i)
   end do
   close(iunit)

   !-----------------------------------------------------------------!
   
   allocate(Sk(N_skt),Vk(N_skt),Rk(N_skt))
   open(unit=iunit,file="vertice_K.dat",status='old',action='read')
   do i = 1,N_skt
      read(iunit,*) Vk(i),Rk(i)
   end do
   do i = 1,N_skt
      read(iunit,*) Sk(i)
   end do
   close(iunit)

   call flush()
   return
end subroutine ler_Vertices

subroutine inicia_Bi
   use var_annealing, only : Ns,Nviz,jviz,Aij,Bi,S,E_tot
   implicit none
   integer :: i,j,k

   allocate(Bi(Ns))

   E_tot = 0.0d0
   Bi(:) = 0.0d0
   do i = 1,Ns
      do k = Nviz(i-1)+1,Nviz(i)
         j = jviz(k)
         Bi(i) = Bi(i) + S(j)*Aij(k)
      end do
      E_tot = E_tot + S(i)*Bi(i)
   end do
   E_tot = 0.5d0*E_tot

   return
end subroutine inicia_Bi

subroutine diretorios
   use var_annealing
   implicit none
   integer :: i
   logical :: direx

   dir1 = 'Resultados/'
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   dir1 = 'Resultados/Annealing/'
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   direx = .true.
   i = 0
   do while(direx .eqv. .true.)
      i = i + 1
      write(dir2,"('Resultados/Annealing/Simula_',I4.4,'/')") i
      inquire(file=dir2,exist=direx)
   end do
   call system("mkdir " // trim(dir2))

   return
end subroutine diretorios

subroutine update(i,dE)
   use var_annealing, only : Nviz,jviz,Aij,Bi,S,E_tot
   implicit none
   integer, intent(in) :: i
   real(8), intent(in) :: dE
   integer :: j,k
   real(8) :: dBi

   do k = Nviz(i-1)+1,Nviz(i)
      j = jviz(k)
      dBi = -2.0d0*S(i)*Aij(k)
      Bi(j) = Bi(j) + dBi
   end do

   S(i) = -S(i)
   E_tot = E_tot + dE
   
   return
end subroutine update

subroutine Monte_Carlo(temps)
   use var_annealing, only : Ns,N_mc,beta,iunit_conf,iunit_en,E_tot
   implicit none
   real(8), intent(in) :: temps
   integer :: imc
   real(8) :: E1,E2

   beta = 1.0d0/temps

   ! Termalização !
   do imc = 1,N_mc
      call metropolis
   end do

   ! Amostragem !
   call config_S(iunit_conf,0)
   call En_save(iunit_en,1)

   E1 = 0.0d0
   E2 = 0.0d0

   do imc = 1,N_mc
      call metropolis
      call samples(temps)
      E1 = E1 + E_tot
      E2 = E2 + E_tot**2

      if (mod(imc,10)==0) then
         call config_S(40,2)
      end if
   end do

   E1 = E1/real(N_mc)
   E2 = E2/real(N_mc)

   write(50,*) temps,E1/real(Ns,8),(E2 - E1**2)/real(Ns,8)*beta**2

   call config_S(iunit_conf,3)
    
   return
end subroutine Monte_Carlo

subroutine metropolis
   use var_annealing, only : Ns,Bi,S,beta
   !use mtmod, only : igrnd,grnd
   use ranutil, only : ranmar,rand_int
   implicit none
   integer :: i,im
   real(8) :: dE

   do im = 1,Ns
      i = rand_int(1,Ns)
      dE = -2.0d0*S(i)*Bi(i)
      if (ranmar() .lt. exp(-beta*dE)) then
         call update(i,dE)
      end if
   end do

   return
end subroutine metropolis

subroutine samples(temps)
   use var_annealing, only : temp,iunit_conf,iunit_en
   implicit none
   real(8), intent(in) :: temps

   temp = temps

   !call config_S(iunit_conf,1)
   call En_save(iunit_en,2)

   return
end subroutine samples

subroutine config_S(iunit,flag)
   use var_annealing, only : Ns,N_mc,rx,ry,mx,my,S,temp,dir2
   implicit none
   integer, intent(in) :: iunit,flag
   integer :: i
   character(60) :: nome

   if (flag == 0) then
      ! Inicia a configuração !
      write(nome,"('config_temp_',f8.2,'.dat')") temp
      open(unit=iunit,file=trim(dir2) // trim(nome))

      write(iunit,*) Ns,N_mc
      do i = 1,Ns
         write(iunit,*) rx(i),ry(i),mx(i),my(i)
      end do
   else if (flag == 1) then
      do i = 1,Ns
         write(iunit,*) (S(i)+1)/2
      end do
   else if (flag == 2) then
      write(iunit,*) Ns
      write(iunit,*) ' '
      do i = 1,Ns
         write(iunit,"(4(2x,f10.5))") rx(i),ry(i),S(i)*mx(i),S(i)*my(i)
      end do
   else if (flag == 3) then
      call flush()
      close(iunit)
   end if

   return
end subroutine config_S

subroutine En_save(iunit,flag)
   use var_annealing, only : Ns,N_mc,N_temp,temp,E_tot,dir2
   implicit none
   integer, intent(in) :: iunit,flag

   if (flag == 0) then
      open(unit=iunit,file=(trim(dir2) // 'energia.dat'))
      write(iunit,*) Ns,N_mc,N_temp
   else if (flag == 1) then
      write(iunit,*) temp
   else if (flag == 2) then
      write(iunit,*) E_tot
   else if (flag == 3) then
      call flush()
      close(iunit)
   end if

   return
end subroutine En_save

subroutine Annealing
   use var_annealing, only : N_temp,Ti,dT,temp
   implicit none
   integer :: i_temp
   real(8) :: temps

   open(40,file='config.xyz')
   open(50,file='en_med.dat')

   do i_temp = 1,N_temp
      temps = Ti + (i_temp-1)*dT
      temp = temps
      call Monte_Carlo(temps)
   end do

   close(40)
   close(50)

   return
end subroutine Annealing

program main
   
   call inicial

   call Annealing

end program main


