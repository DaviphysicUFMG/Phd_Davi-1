!---------------------------------------------------------------------------------!
!           GERADOR DE NÚMEROS ALEATÓRIOS MT19937 (Mersenne Twister Random Numbers)
! module mtdefs
!     integer,parameter :: rk=selected_real_kind(10,40)
!   end module mtdefs
  
!   module mtmod
  
!     use mtdefs, only : rk_mtmod => rk
!     implicit none
!     ! integer,parameter :: rk_mtmod=selected_real_kind(10,40) ! This is double precision
!     integer,parameter,private :: defaultsd = 4357    ! Default seed  
!     integer,parameter,private :: N = 624, N1 = N + 1 ! Period parameters
!     integer,dimension(0:N-1),private :: mt      ! Array for the state vector
!     integer,private :: mti = N1
!     interface mtsave
!        module procedure mtsavef
!        module procedure mtsaveu
!     end interface
!     interface mtget
!        module procedure mtgetf
!        module procedure mtgetu
!     end interface
    
!   contains
  
!     subroutine sgrnd(seed)
!       implicit none
!       ! Setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!       ! [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]
!       integer, intent(in) :: seed    
!       mt(0) = iand(seed,-1)
!       do mti=1,N-1
!          mt(mti) = iand(69069 * mt(mti-1),-1)
!       enddo
!       return
!     end subroutine sgrnd
  
!     function grnd()
!       implicit none
!       real(rk_mtmod) :: grnd    
!       ! Period parameters
!       integer, parameter :: M = 397, MATA  = -1727483681 ! constant vector a
!       integer, parameter :: LMASK =  2147483647          ! least significant r bits
!       integer, parameter :: UMASK = -LMASK - 1           ! most significant w-r bits
!       ! Tempering parameters
!       integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544    
!       integer,save :: mag01(0:1)=[0,MATA] ! mag01(x) = x * MATA for x=0,1
  
!       integer :: kk,y
  
!       if (mti>=N) then               ! generate N words at one time
!          if (mti==N+1) then          ! if sgrnd() has not been called,
!             call sgrnd( defaultsd ) ! a default initial seed is used
!          endif
!          do kk=0,N-M-1
!             y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
!             mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
!          enddo
!          do kk=N-M,N-2
!             y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
!             mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
!          enddo
!          y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
!          mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
!          mti = 0
!       endif
      
!       y=mt(mti)
!       mti = mti + 1 
!       y=ieor(y,TSHFTU(y))
!       y=ieor(y,iand(TSHFTS(y),TMASKB))
!       y=ieor(y,iand(TSHFTT(y),TMASKC))
!       y=ieor(y,TSHFTL(y))
      
!       if (y<0) then
!          grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
!       else
!          grnd=dble(y)/(2.0d0**32-1.0d0)
!       endif
      
!       return
  
!     contains
  
!       integer function TSHFTU(y)
!         integer,intent(in) :: y
!         TSHFTU=ishft(y,-11)
!         return
!       end function TSHFTU
!       integer function TSHFTS(y)
!         integer,intent(in) :: y
!         TSHFTS=ishft(y,7)
!         return
!       end function TSHFTS
!       integer function TSHFTT(y)
!         integer,intent(in) :: y
!         TSHFTT=ishft(y,15)
!         return
!       end function TSHFTT
!       integer function TSHFTL(y)
!         integer,intent(in) :: y
!         TSHFTL=ishft(y,-18)
!         return
!       end function TSHFTL
  
!     end function grnd
  
  
!     integer function igrnd(l,h)
!       implicit none
!       integer,intent(in) :: l,h
!       real(rk_mtmod) :: u,r
!       u=grnd()
!       r=(h-l+1)*u+l
!       igrnd=int(r)
!       return
!     end function igrnd
  
!     function gaussrnd()
!       implicit none
!       real(rk_mtmod) :: gaussrnd
!       real(rk_mtmod) :: fac,v1,v2,r
!       real(rk_mtmod), save :: gset
!       integer, save :: iset=0
  
!       if (iset==0) then ! Create a new RN
!          r=100.0
!          do while (r>1.0)
!             v1 = 2.0*grnd()-1.0
!             v2 = 2.0*grnd()-1.0
!             r = v1*v1+v2*v2
!          end do
!          fac = sqrt(-2.0*log(r)/r)
!          gset = v1*fac
!          gaussrnd = v2*fac
!          iset = 1
!       else ! Use the 2nd NR from the previous call
!          gaussrnd = gset
!          iset = 0
!       endif
!       return
!     end function gaussrnd
    
    
  
!     !---------------------------------------------------------------!
!     !                                                               !
!     ! State saving subroutines.                                     !
!     !                                                               !
!     !  Usage:  call mtsave( file_name, format_character )           !
!     !     or   call mtsave( unit_number, format_character )         !
!     !  where   format_character = 'u' or 'U' will save in           !
!     !          unformatted form, otherwise state information will   !
!     !          be written in formatted form.                        !
!     !                                                               !
!     !---------------------------------------------------------------!
  
!     subroutine mtsavef( fname, forma )
!       implicit none
  
!       character(*), intent(in) :: fname
!       character, intent(in)    :: forma
      
!       select case (forma)
!       case('u','U')
!          open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED',position='APPEND')
!          write(10) mti
!          write(10) mt       
!       case default
!          open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED',position='APPEND')
!          write(10,*) mti
!          write(10,*) mt       
!       end select
  
!       close(10)    
!       return
!     end subroutine mtsavef
  
!     subroutine mtsaveu(unum,forma)
!       implicit none
!       integer, intent(in)    :: unum
!       character, intent(in)  :: forma
      
!       select case (forma)
!       case('u','U')
!          write(unum) mti
!          write(unum) mt       
!       case default
!          write(unum,*) mti
!          write(unum,*) mt       
!       end select
      
!       return
!     end subroutine mtsaveu
  
!     subroutine mtgetf(fname,forma)
!       implicit none
!       character(*), intent(in) :: fname
!       character, intent(in)    :: forma
      
!       select case (forma)
!       case('u','U')
!          open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
!          read(10) mti
!          read(10) mt
!       case default
!          open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
!          read(10,*) mti
!          read(10,*) mt
!       end select
  
!       close(10)
!       return
!     end subroutine mtgetf
  
!     subroutine mtgetu(unum,forma)
!       implicit none
!       integer, intent(in)    :: unum
!       character, intent(in)  :: forma
      
!       select case (forma)
!       case('u','U')
!          read(unum) mti
!          read(unum) mt
!       case default
!          read(unum,*) mti
!          read(unum,*) mt
!       end select
      
!       return
!     end subroutine mtgetu  
  
!     integer function getseed(info,file)
    
!       implicit none
!       integer,optional,intent(in) :: info,file
!       integer :: t(8),rn,is
!       integer,parameter :: LMASK=huge(rn) ! = 0111...111
!       integer,parameter :: LUN=676769
!       character (len=80) :: rdev0='/dev/urandom',rdev1='/dev/random',rdev
!       logical :: openok,readok,printinfo
  
!       openok=.true.
!       readok=.true.
  
!       if (present(file)) then
!          if (file==0) then
!             rdev=rdev0
!          else
!             rdev=rdev1
!          end if
!       else
!          rdev=rdev0
!       end if
!       if (present(info)) then
!          printinfo=(info/=0)
!       else
!          printinfo=.false.
!       end if
  
!       open(LUN,file=rdev,form='unformatted',access='stream',action='read',iostat=is)
!       if (is/=0) then
!          openok=.false.
!          print *,'open',is
!       else
!          read(LUN,iostat=is) rn
!          if (is/=0) then
!             readok=.false.
!          end if
!       end if
!       if (openok) close(LUN)
  
!       if (openok.and.readok) then
!          rn=iand(rn,LMASK) ! Make it positive, i.e. zero the leftmost bit
!          if (printinfo) write(6,'(a,a,a,i0)') 'Seed from ',trim(rdev),': ',rn
!       else
!          call date_and_time(values=t)
!          rn=t(7)+60*(t(6)+60*(t(5)+24*(t(3)-1+31*(t(2)-1+12*t(1)))))+t(8)
!          if (printinfo) write(6,'(a,i12)') 'Seed from time:',rn
!       end if
  
!       getseed=rn
!       return
!     end function getseed
  
! end module mtmod
  !-----------------------------------------------------------------!
  !                          MT19937                                !
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
  !    ...                                                          !                                           
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
   integer :: Ns,N_viz,N_mc,N_temp
   integer, dimension(:), allocatable :: S,Nviz,jviz
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
   use var_annealing, only : Ns,N_viz,N_mc,N_temp,Ti,Tf!,seed
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

   call config_S(iunit_conf,1)
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
   use var_annealing, only : Ns,N_mc,temp,E_tot,dir2
   implicit none
   integer, intent(in) :: iunit,flag

   if (flag == 0) then
      open(unit=iunit,file=(trim(dir2) // 'energia.dat'))
      write(iunit,*) Ns,N_mc
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



