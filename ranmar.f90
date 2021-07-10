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
