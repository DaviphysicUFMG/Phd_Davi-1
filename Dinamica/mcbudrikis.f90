!! Artigo referência: Diversity Enabling Equilibration: Disorder and the Ground State in Artificial Spin Ice
!!                         Zoe Budrikis, Paolo Politi and R. L. Stamps
!!
!!
!!    ListStreamPlot -> data = Partition[#, 2] & /@ data;

!---------------------------------------------------------------------------------!
!           GERADOR DE NÚMEROS ALEATÓRIOS MT19937 (Mersenne Twister Random Numbers)
module mtdefs
   integer,parameter :: rk=selected_real_kind(10,40)
 end module mtdefs
 
 module mtmod
 
   use mtdefs, only : rk_mtmod => rk
   implicit none
   ! integer,parameter :: rk_mtmod=selected_real_kind(10,40) ! This is double precision
   integer,parameter,private :: defaultsd = 4357    ! Default seed  
   integer,parameter,private :: N = 624, N1 = N + 1 ! Period parameters
   integer,dimension(0:N-1),private :: mt      ! Array for the state vector
   integer,private :: mti = N1
   interface mtsave
      module procedure mtsavef
      module procedure mtsaveu
   end interface
   interface mtget
      module procedure mtgetf
      module procedure mtgetu
   end interface
   
 contains
 
   subroutine sgrnd(seed)
     implicit none
     ! Setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
     ! [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]
     integer, intent(in) :: seed    
     mt(0) = iand(seed,-1)
     do mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
     enddo
     return
   end subroutine sgrnd
 
   function grnd()
     implicit none
     real(rk_mtmod) :: grnd    
     ! Period parameters
     integer, parameter :: M = 397, MATA  = -1727483681 ! constant vector a
     integer, parameter :: LMASK =  2147483647          ! least significant r bits
     integer, parameter :: UMASK = -LMASK - 1           ! most significant w-r bits
     ! Tempering parameters
     integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544    
     integer,save :: mag01(0:1)=[0,MATA] ! mag01(x) = x * MATA for x=0,1
 
     integer :: kk,y
 
     if (mti>=N) then               ! generate N words at one time
        if (mti==N+1) then          ! if sgrnd() has not been called,
           call sgrnd( defaultsd ) ! a default initial seed is used
        endif
        do kk=0,N-M-1
           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
           mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        do kk=N-M,N-2
           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
           mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
     endif
     
     y=mt(mti)
     mti = mti + 1 
     y=ieor(y,TSHFTU(y))
     y=ieor(y,iand(TSHFTS(y),TMASKB))
     y=ieor(y,iand(TSHFTT(y),TMASKC))
     y=ieor(y,TSHFTL(y))
     
     if (y<0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
     else
        grnd=dble(y)/(2.0d0**32-1.0d0)
     endif
     
     return
 
   contains
 
     integer function TSHFTU(y)
       integer,intent(in) :: y
       TSHFTU=ishft(y,-11)
       return
     end function TSHFTU
     integer function TSHFTS(y)
       integer,intent(in) :: y
       TSHFTS=ishft(y,7)
       return
     end function TSHFTS
     integer function TSHFTT(y)
       integer,intent(in) :: y
       TSHFTT=ishft(y,15)
       return
     end function TSHFTT
     integer function TSHFTL(y)
       integer,intent(in) :: y
       TSHFTL=ishft(y,-18)
       return
     end function TSHFTL
 
   end function grnd
 
 
   integer function igrnd(l,h)
     implicit none
     integer,intent(in) :: l,h
     real(rk_mtmod) :: u,r
     u=grnd()
     r=(h-l+1)*u+l
     igrnd=int(r)
     return
   end function igrnd
 
   function gaussrnd()
     implicit none
     real(rk_mtmod) :: gaussrnd
     real(rk_mtmod) :: fac,v1,v2,r
     real(rk_mtmod), save :: gset
     integer, save :: iset=0
 
     if (iset==0) then ! Create a new RN
        r=100.0
        do while (r>1.0)
           v1 = 2.0*grnd()-1.0
           v2 = 2.0*grnd()-1.0
           r = v1*v1+v2*v2
        end do
        fac = sqrt(-2.0*log(r)/r)
        gset = v1*fac
        gaussrnd = v2*fac
        iset = 1
     else ! Use the 2nd NR from the previous call
        gaussrnd = gset
        iset = 0
     endif
     return
   end function gaussrnd
   
   
 
   !---------------------------------------------------------------!
   !                                                               !
   ! State saving subroutines.                                     !
   !                                                               !
   !  Usage:  call mtsave( file_name, format_character )           !
   !     or   call mtsave( unit_number, format_character )         !
   !  where   format_character = 'u' or 'U' will save in           !
   !          unformatted form, otherwise state information will   !
   !          be written in formatted form.                        !
   !                                                               !
   !---------------------------------------------------------------!
 
   subroutine mtsavef( fname, forma )
     implicit none
 
     character(*), intent(in) :: fname
     character, intent(in)    :: forma
     
     select case (forma)
     case('u','U')
        open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED',position='APPEND')
        write(10) mti
        write(10) mt       
     case default
        open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED',position='APPEND')
        write(10,*) mti
        write(10,*) mt       
     end select
 
     close(10)    
     return
   end subroutine mtsavef
 
   subroutine mtsaveu(unum,forma)
     implicit none
     integer, intent(in)    :: unum
     character, intent(in)  :: forma
     
     select case (forma)
     case('u','U')
        write(unum) mti
        write(unum) mt       
     case default
        write(unum,*) mti
        write(unum,*) mt       
     end select
     
     return
   end subroutine mtsaveu
 
   subroutine mtgetf(fname,forma)
     implicit none
     character(*), intent(in) :: fname
     character, intent(in)    :: forma
     
     select case (forma)
     case('u','U')
        open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
        read(10) mti
        read(10) mt
     case default
        open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
        read(10,*) mti
        read(10,*) mt
     end select
 
     close(10)
     return
   end subroutine mtgetf
 
   subroutine mtgetu(unum,forma)
     implicit none
     integer, intent(in)    :: unum
     character, intent(in)  :: forma
     
     select case (forma)
     case('u','U')
        read(unum) mti
        read(unum) mt
     case default
        read(unum,*) mti
        read(unum,*) mt
     end select
     
     return
   end subroutine mtgetu  
 
   integer function getseed(info,file)
   
     implicit none
     integer,optional,intent(in) :: info,file
     integer :: t(8),rn,is
     integer,parameter :: LMASK=huge(rn) ! = 0111...111
     integer,parameter :: LUN=676769
     character (len=80) :: rdev0='/dev/urandom',rdev1='/dev/random',rdev
     logical :: openok,readok,printinfo
 
     openok=.true.
     readok=.true.
 
     if (present(file)) then
        if (file==0) then
           rdev=rdev0
        else
           rdev=rdev1
        end if
     else
        rdev=rdev0
     end if
     if (present(info)) then
        printinfo=(info/=0)
     else
        printinfo=.false.
     end if
 
     open(LUN,file=rdev,form='unformatted',access='stream',action='read',iostat=is)
     if (is/=0) then
        openok=.false.
        print *,'open',is
     else
        read(LUN,iostat=is) rn
        if (is/=0) then
           readok=.false.
        end if
     end if
     if (openok) close(LUN)
 
     if (openok.and.readok) then
        rn=iand(rn,LMASK) ! Make it positive, i.e. zero the leftmost bit
        if (printinfo) write(6,'(a,a,a,i0)') 'Seed from ',trim(rdev),': ',rn
     else
        call date_and_time(values=t)
        rn=t(7)+60*(t(6)+60*(t(5)+24*(t(3)-1+31*(t(2)-1+12*t(1)))))+t(8)
        if (printinfo) write(6,'(a,i12)') 'Seed from time:',rn
     end if
 
     getseed=rn
     return
   end function getseed
 
 
end module mtmod
 !-----------------------------------------------------------------!
 !                                                                 !
 ! Usage example:                                                  !
 !                                                                 !
 !    use mtdefs, only : rk                                        !
 !       ! NOTE: 'rk' is the real number kind                      !
 !       ! defined in module mtdefs (see above).                   !
 !       ! Alternatively (see the definitions in the               !
 !       ! beginning of the module) you can define the constant as !
 !       ! integer,parameter :: rk=selected_real_kind(10,40)       !
 !    ...                                                          !
 !    use mtmod                                                    !
 !    ...                                                          !
 !    integer :: seed,i                                            !
 !    real(rk) :: x,g                                              !
 !    ...                                                          !
 !    [ seed=getseed() ]                                           !
 !    call sgrnd(seed)                                             !
 !    x=grnd()                                                     !
 !    g=gaussrnd()                                                 !
 !    i=igrnd(0,100)                                               !
 !    ...                                                          !
 !                                                                 !
 !-----------------------------------------------------------------!

module var_global
   real(8), parameter :: pi = 4.0d0*atan(1.0d0)
   integer :: Ns,N_viz,N_campo,ncor,Nmssf
   integer, dimension(:), allocatable :: S,Nviz,jviz,cor
   real(8), dimension(:), allocatable :: Aij,Bi,Hi,hcut
   real(8), dimension(:), allocatable :: x,y,mx,my
   real(8), dimension(:), allocatable :: Bx,By,theta,ModB
   integer,dimension(:,:),allocatable :: cortab
   real(8) :: BMax,delB,ExtMax,theta_c,D
   character(600) :: dir1,dir2
   integer :: kk
end module var_global

subroutine inicial
   use var_global
   implicit none
   integer :: qual
   logical :: direx

   !!Ler input

   print*, "Histerese = 0; Campo Rotativo = 1"
   read(*,*) qual
   print*, "Número para MSSF:"; read(*,*) Nmssf

   call ler_input(1)
   call ler_config(2)
   call ler_Aij(3)
   call ler_Nviz(4)
   call ler_cortab(5)
   !!Calculo dos campos iniciais
   allocate(Bi(Ns),Hi(Ns))
   call campo_dip
   call campo_flip

   

   if (qual == 0) then
      call campo_zee_H
   else if (qual == 1) then
      call campo_zee
   else
      print*, "Erro"
      stop
   end if
   
   !!Diretórios

   write(dir1,"('Resultados/')")
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   write(dir1,"('Resultados/His/')")
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   write(dir1,"('Resultados/Rot/')")
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   if (qual == 0) then
      kk = 0
      direx = .true.
      do while (direx .eqv. .true.)
         kk = kk + 1
         write(dir2,"('Resultados/His/Simulacao_',I3.3,'/')") kk
         inquire(file=dir2,exist=direx)
      end do
   else
      kk = 0
      direx = .true.
      do while (direx .eqv. .true.)
         kk = kk + 1
         write(dir2,"('Resultados/Rot/Simulacao_',I3.3,'/')") kk
         inquire(file=dir2,exist=direx)
      end do
   end if
   call system ("mkdir " // trim(dir2))

   return

end subroutine inicial

subroutine ler_input(unit_in)
   use mtmod, only : getseed,sgrnd
   use var_global, only : Ns,N_viz,BMax,delB,ExtMax,theta_c,D
   implicit none
   integer, intent(in) :: unit_in
   integer :: seed
   
   open(unit=unit_in,file='input.dat',status='old',action='read')
   read(unit_in,*) D
   read(unit_in,*) Ns
   read(unit_in,*) N_viz
   read(unit_in,*) BMax
   read(unit_in,*) delB
   read(unit_in,*) ExtMax
   read(unit_in,*) theta_c
   close(unit_in)

   seed = getseed()
   print*, seed
   call sgrnd(seed)  

   return
end subroutine ler_input

subroutine ler_config(unit_in)
   use mtmod, only : grnd
   use var_global, only : Ns,S,x,y,mx,my,cor
   implicit none
   integer, intent(in) :: unit_in
   integer :: i
   real(8) :: z

   allocate(S(Ns),x(Ns),y(Ns),mx(Ns),my(Ns),cor(Ns))
   open(unit=unit_in,file='config0.xyz',status='old',action='read')
   read(unit_in,*) i
   read(unit_in,*)
   do i = 1,Ns
      read(unit_in,*) x(i),y(i),mx(i),my(i),cor(i)
   end do
   close(unit_in)

   S(:) = 1
   do i = 1,Ns
      z = grnd()
      if (z .lt. 0.5d0) then
         S(i) = -1
      end if
   end do

   return
end subroutine ler_config

subroutine ler_Aij(unit_in)
   use var_global, only : N_viz,Aij,jviz
   implicit none
   integer, intent(in) :: unit_in
   integer :: i

   allocate(Aij(N_viz),jviz(N_viz))
   open(unit=unit_in,file='Aij.dat',status='old',action='read')
   do i = 1,N_viz
      read(unit_in,*) jviz(i),Aij(i)
   end do
   close(unit_in)

   return
end subroutine ler_Aij

subroutine ler_Nviz(unit_in)
   use var_global, only : Ns,Nviz
   implicit none
   integer, intent(in) :: unit_in
   integer :: i

   allocate(Nviz(0:Ns))
   open(unit=unit_in,file='Nviz.dat',status='old',action='read')

   Nviz(0) = 0
   do i = 1,Ns
      read(unit_in,*) Nviz(i)
   end do

   close(unit_in)
   return
end subroutine ler_Nviz

subroutine ler_cortab(unit_in)
   use var_global, only : ncor,cortab
   implicit none
   integer,intent(in) :: unit_in
   integer :: i

   open(unit=unit_in,file='cor.dat',status='old',action='read')
   read(unit_in,*) ncor
   allocate(cortab(-1:1,ncor))
   cortab = 0
   do i = 1,ncor
      read(unit_in,*) cortab(1,i),cortab(-1,i)
   end do
   close(unit_in) 

end subroutine ler_cortab

subroutine campo_dip
   use var_global, only : Ns,S,Nviz,jviz,Aij,Bi
   implicit none
   integer :: i,j,k

   Bi(:) = 0.0d0
   do i = 1,Ns
      do k = Nviz(i-1)+1,Nviz(i)
         j = jviz(k)
         Bi(i) = Bi(i) + S(j)*Aij(k)
      end do
   end do

   return
end subroutine campo_dip

subroutine campo_zee_H
   use var_global, only : N_campo,ExtMax,Bx,By,pi,theta_c,ModB
   implicit none
   integer :: i
   real(8) :: dB,Bmod

   N_campo = 1000!256

   allocate(Bx(N_campo),By(N_campo),ModB(N_campo))
   dB = 2.0d0*ExtMax/real(N_campo-1,8)

   do i = 1,N_campo/2
      Bmod = -ExtMax + 2.0d0*dB*(i-1)
      ModB(i) = Bmod
      Bx(i) = Bmod*cos(theta_c*pi/180.0d0)
      By(i) = Bmod*sin(theta_c*pi/180.0d0)
   end do

   do i = N_campo/2+1,N_campo
      Bmod = ExtMax - 2.0d0*dB*(i-1-N_campo/2)
      ModB(i) = Bmod
      Bx(i) = Bmod*cos(theta_c*pi/180.0d0)
      By(i) = Bmod*sin(theta_c*pi/180.0d0)
   end do

   return
end subroutine campo_zee_H

subroutine campo_zee
   use var_global, only : N_campo,ExtMax,Bx,By,pi,theta,ModB
   implicit none
   integer :: i

   N_campo = 362!256

   allocate(Bx(1:N_campo),By(1:N_campo),theta(1:N_campo),ModB(1:N_campo))

   do i = 1,N_campo
      theta(i) = i*pi/180.0d0
      Bx(i) = ExtMax*cos(i*pi/180.0d0)
      By(i) = ExtMax*sin(i*pi/180.0d0)
      ModB(i) = sqrt(Bx(i)**2 + By(i)**2)
   end do

   return
end subroutine campo_zee

subroutine campo_flip
   use mtmod, only : gaussrnd
   use var_global, only : Ns,BMax,delB,hcut
   implicit none
   integer :: i

   allocate(hcut(Ns))
   do i = 1,Ns
      hcut(i) = delB*gaussrnd() + Bmax
   end do

   return
end subroutine campo_flip

subroutine abrir_fechar(i,a,b,c)
   use var_global, only : dir2
   implicit none
   integer,intent(in) :: i,a,b,c

   if (i .eq. 1) then
      open(unit = a,file=trim(dir2) // trim("conf.xyz"))
      open(unit = b,file=trim(dir2) // trim("mag.dat"))
      open(unit = c,file=trim(dir2) // trim("campo.xyz"))

      open(100,file=trim(dir2) // trim("Ising.dat"))

   else
      close(a)
      close(b)
      close(c)
   end if

   return
end subroutine abrir_fechar

subroutine budrikis(try)
   use mtmod, only : igrnd
   use var_global, only : Ns,S,Bi,Hi,hcut
   implicit none
   integer, intent(out) :: try
   integer :: SetPos(Ns)
   integer :: i
   real(8) :: E_try

   !!Percorre toda a rede determinando todo o conjunto de spins que são flipáveis!!
   try = 0
   do i = 1,Ns
      E_try = S(i)*(Bi(i) - Hi(i))
      if (E_try .gt. hcut(i)) then !.and. E_try.gt.maxE
         !maxE = E_try
         !j = i
         try = try + 1    !! Como documentado no artigo de Zoe Budrikis em 2012
         SetPos(try) = i
      end if
   end do
   !call update(j)
   !!Flipa um spin aleatório pertencente ao conjunto e atualiza os campos!!
   if (try .ne. 0) then
      i = igrnd(1,try)
      call update(SetPos(i))
   end if

   return
end subroutine budrikis

subroutine update(i)
   use var_global, only : S,Nviz,Aij,Bi,jviz
   implicit none
   integer,intent(in) :: i
   integer :: j,k
   real(8) :: dB

   !! Atualiza o Campo de i nos outros Spin!!
   do k = Nviz(i-1)+1,Nviz(i)
      j = jviz(k)
      dB = -2.0d0*S(i)*Aij(k)
      Bi(j) = Bi(j) + dB
   end do
   !! Atualiza o sentido do Spin i
   S(i) = -S(i)
   
   return
end subroutine update

subroutine config(u_conf)
   use var_global, only : Ns,x,y,mx,my,S,Bi,cor,cortab
   implicit none
   integer, intent(in) :: u_conf
   integer :: i

   write(u_conf,*) Ns
   write(u_conf,*) " "
   do i = 1,Ns
      write(u_conf,*) x(i),y(i),S(i)*mx(i),S(i)*my(i),S(i)*Bi(i),cortab(S(i),cor(i))
   end do

   return
end subroutine config

subroutine config_ising(u_conf)
   use var_global, only: Ns,S
   implicit none
   integer,intent(in) :: u_conf
   integer :: i

   do i = 1,Ns
      write(u_conf,*) (S(i)+1)/2
   end do

   return
end subroutine config_ising

subroutine mag(ibx,iby,MB)
   use var_global
   implicit none
   real(8), intent(in) :: ibx,iby,MB
   real(8) :: Magx,Magy

   Magx = sum(S(:)*mx(:))
   Magy = sum(S(:)*my(:))
   write(30,*) ibx,Magx/real(Ns,8),iby,Magy/real(Ns,8),MB,(Magx*cos(theta_c*pi/180.0d0) + Magy*sin(theta_c*pi/180.0d0))/real(Ns,8)

   return
end subroutine mag

subroutine mssf(iunit,flag)
   use var_global, only: Ns,x,y,mx,my,S,Nmssf,dir2
   implicit none
   integer, intent(in) :: iunit,flag
   integer :: i

   if (flag==1) then
      open(unit = iunit,file=trim(dir2) // trim("mssf.dat"))
      write(iunit,*) Ns,Nmssf
      do i = 1,Ns
         write(iunit,*) x(i),y(i),mx(i),my(i)
      end do
   else if (flag==2) then
      do i = 1,Ns
         write(iunit,*) (S(i)+1)/2
      end do
   else
      close(iunit)
   end if

   return
end subroutine mssf

subroutine one_rev
   use var_global, only : Hi,Bx,By,mx,my,N_campo
   implicit none
   integer :: i,try

   do i = 1,N_campo
      try = 1
      Hi(:) = Bx(i)*mx(:) + By(i)*my(:)
      do while (try>0)
         call budrikis(try)
         !call config(40)
      end do
   end do

end subroutine one_rev

subroutine simu_rev
   use var_global, only : N_campo,Hi,Bx,By,mx,my,ModB,y,Nmssf
   implicit none
   integer :: i,try,N,j
   real(8) :: dM

   dM = 360.0d0/real(N_campo)
   N = N_campo/Nmssf

   j = 0
   do i = 1,N_campo
      try = 1
      Hi(:) = Bx(i)*mx(:) + By(i)*my(:)
      
      do while (try>0)
         call budrikis(try)
         call config_ising(100)
         !call config(40)
      end do

      call config(20)
      call mag(Bx(i),By(i),ModB(i))

      !if (mod(i,N)==0) then
         call mssf(50,2)
      !   j = j + 1
      !end if

      write(40,*) 1
      write(40,*) ModB(i)
      write(40,*) 0.,maxval(y)+1.5d0,Bx(i)/abs(ModB(i)),By(i)/abs(ModB(i))

      call flush()
   end do

   print*, j

end subroutine simu_rev

program main
   implicit none
   integer :: i

   call inicial
   call abrir_fechar(1,20,30,40)
   call mssf(50,1)

   ! Primeira revolução do Campo - Não salva resultados !
   call one_rev

   do i = 1,5
      call simu_rev
      print*, i
   end do

   call abrir_fechar(0,20,30,40)
   call mssf(50,3)

end program main