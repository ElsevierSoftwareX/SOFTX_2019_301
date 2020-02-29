!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: teutils.f90,v 1.0 19-03-2018, IBU
!
!                This source code is part of
!
!   Symbolic Information Flow Measure Code for studying the
!   information flow in dynamical systems
!
!                        VERSION 1.0
!
! Written by Hiqmet Kamberaj.
! Copyright (C) 2018 Hiqmet Kamberaj.
! Check out h.kamberaj@gmail.com for more information.
!
! This program is free software; you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the Free Software Foundation; 
! GPL-3.0
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this program; 
! if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
! Boston, MA 02111-1307 USA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE TEUTILS_CLASS
use sifm_kinds
use TERANDOM_CLASS
implicit none
!

! Default input and output units:
   integer, parameter :: DEFAULT_INPUT_UNIT = 5
   integer, parameter :: DEFAULT_OUTPUT_UNIT = 6
! Number and value of preconnected units
   integer, parameter :: NUMBER_OF_PRECONNECTED_UNITS = 3
   integer, parameter :: PRECONNECTED_UNITS (NUMBER_OF_PRECONNECTED_UNITS) = (/ 0, 5, 6 /)
! Largest allowed unit number (or a large number, if none)
   integer, parameter :: MAX_UNIT_NUMBER = 1000
   INTEGER, PARAMETER :: NPAR_POLY=8
   interface swap 
        module procedure iswap, iswap_vec, rswap, rswap_vec
   end interface
   INTERFACE poly
     MODULE PROCEDURE poly_s,poly_v
   END INTERFACE
   INTERFACE erfcc
     MODULE PROCEDURE erfcc_s,erfcc_v
   END INTERFACE
   interface TransferFunctionInn
       module procedure TransferFunctionInn_2d, TransferFunctionInn_1d
   end interface
   interface TransferFunctionInn_derivative
       module procedure TransferFunctionInn_derivative_2d, TransferFunctionInn_derivative_1d
   end interface
   interface write_csvs
       module procedure writei_csvs, writec_csvs 
   end interface
   interface write_xyzs
       module procedure writei_xyzs, writec_xyzs 
   end interface  
   interface read_csvs
       module procedure readi_csvs, readc_csvs 
   end interface
   interface read_xyzs
       module procedure readi_xyzs, readc_xyzs 
   end interface   
CONTAINS

subroutine mteallocate_r(x, n, var, routine)
real(sifm_real), allocatable, intent(inout) :: x(:)
integer, intent(in) :: n
character(len=*) :: var, routine
integer :: ierror
Allocate( x(n), stat=ierror)
IF (ierror /= 0) then
    write(*,'("Error allocating memory for ", A20, "  in subroutine ", A30)') &
         var, routine
    stop
endif  
return
end subroutine mteallocate_r
!
subroutine mteallocate_i(x, n, var, routine)
integer, allocatable, intent(inout) :: x(:)
integer, intent(in) :: n
character(len=*) :: var, routine
integer :: ierror
Allocate( x(n), stat=ierror)
IF (ierror /= 0) then
    write(*,'("Error allocating memory for ", A20, "  in subroutine ", A30)') &
         var, routine
    stop
endif  
return
end subroutine mteallocate_i
!
subroutine mtedeallocate_r(x, n, var, routine)
integer, intent(in) :: n
real(sifm_real), intent(inout), allocatable :: x(:)
character(len=*) :: var, routine
integer :: ierror
DeAllocate( x, stat=ierror)
IF (ierror /= 0) then
    write(*,'("Error Deallocating memory for ", A20, "  in subroutine ", A30)') &
         var, routine
    stop
endif  
return
end subroutine mtedeallocate_r
!
subroutine mtedeallocate_i(x, n, var, routine)
integer, intent(in) :: n
integer, intent(inout), allocatable :: x(:)
character(len=*) :: var, routine
integer :: ierror
deAllocate( x, stat=ierror)
IF (ierror /= 0) then
    write(*,'("Error deallocating memory for ", A20, "  in subroutine ", A30)') &
         var, routine
    stop
endif  
return
end subroutine mtedeallocate_i
!
function new_unit()  result (result)
      integer :: result
      logical :: exists, opened
      integer :: ios
      do result = 11, max_unit_number
         if (result == DEFAULT_INPUT_UNIT .or. &
             result == DEFAULT_OUTPUT_UNIT) cycle
         if (any (result == PRECONNECTED_UNITS)) cycle
         inquire (unit = result,  &
                  exist = exists,  &
                  opened = opened,  &
                  iostat = ios)
         if (exists .and. .not. opened .and. ios == 0) return
      end do
      result = -1
end function new_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Shuffle the data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine iswap(x,y)
   integer, intent(inout) :: x,y
!   
   integer :: itemp
!
   itemp = x
   x = y
   y = itemp
!  
   return
end subroutine iswap
!
subroutine iswap_vec(D,x,y)
   integer, intent(in) :: D
   integer, dimension(D), intent(inout) :: x,y
!   
   integer, dimension(D) :: itemp
!
   itemp = x
   x = y
   y = itemp
!  
   return
end subroutine iswap_vec
!
subroutine rswap(x,y)
   real(sifm_real), intent(inout) :: x,y
!   
   real(sifm_real) :: itemp
!
   itemp = x
   x = y
   y = itemp
!  
   return
end subroutine rswap
!
subroutine rswap_vec(D,x,y)
   integer, intent(in) :: D
   real(sifm_real), dimension(D), intent(inout) :: x,y
!   
   real(sifm_real), dimension(D) :: itemp
!
   itemp = x
   x = y
   y = itemp
!  
   return
end subroutine rswap_vec
!
subroutine rperm(N0, N, p, IS)
   integer, intent(in) :: N0, N, IS
   integer, dimension(:), intent(out) :: p
!
   integer :: i, N1
   integer :: k, j, ipj, m
   real(sifm_real), dimension(100) :: u
!
   p = (/ (i, i=N0,N) /)
   N1=N-N0+1
! Generate up to 100 U(0,1) numbers at a time.
   do i=1,N1,100
     m = min(N1-i+1, 100)
     call randf_vec(IS,100,u)
     do j=1,m
        ipj = i+j-1
        k = int(u(j)*(N1-ipj+1)) + ipj
        call iswap(p(ipj), p(k))
     end do
   end do
   return
end subroutine rperm
!
subroutine BlockShuffle1D(N,x,m,tau,IS)
integer, intent(in)  :: N,is
integer, intent(in)  :: m,tau
integer, intent(inout) :: x(n)
!
integer, dimension(N) :: indx
integer :: i,k,L
integer, dimension(M) :: I1, I2
integer, dimension(M) :: Temp
real(sifm_real) :: Ur
!
DO i = (M-1)*Tau+1, N
   Ur = randf(IS)
   L = int( ( N - (i + M) )*Ur + i + M)
   I1 = (/(i-(k-1)*Tau, k=1, M)/)
   I2 = (/(L-(k-1)*Tau, k=1, M)/)
   Temp( 1:M )  = X( I1(1:M) )
   X( I1(1:M) ) = X( I2(1:M) )
   X( I2(1:M) ) = Temp( 1:M )
ENDDO
!
return
end subroutine BlockShuffle1D
!
subroutine BlockShuffleND(N,ndim,x,m,tau,IS)
integer, intent(in)  :: N,ndim,is
integer, dimension(ndim), intent(in) :: m,tau
character(len=1), intent(inout) :: x(ndim,n)
!
integer :: i, k, d, L
integer, dimension( imax1d(M) ) :: I1, I2
character(len=1), dimension( imax1d(M) ) :: Temp
real(sifm_real) :: Ur
!
Do d = 1, Ndim
   DO i = (M(d)-1)*Tau(d)+1, N-M(d)
      Ur = randf(IS)
      L = int( ( N - (i + M(d)) )*Ur + i + M(d))
      I1 = (/(i-(k-1)*Tau(d), k=1, M(d))/)
      I2 = (/(L-(k-1)*Tau(d), k=1, M(d))/)
      Temp( 1:M(d) )    = X( d, I1(1:M(d)) )
      X( d,I1(1:M(d)) ) = X( d, I2(1:M(d)) )
      X( d,I2(1:M(d)) ) = Temp( 1:M(d) )
   ENDDO
ENDDO
!
return
end subroutine BlockShuffleND
!
subroutine BlockShuffleND2(N,ndim,x,m,tau,IS)
integer, intent(in)  :: N,ndim,is
integer, dimension(ndim), intent(in) :: m,tau
real(sifm_real), intent(inout) :: x(ndim,n)
!
integer :: i, k, d, L
integer, dimension( imax1d(M) ) :: I1, I2
real(sifm_real), dimension( imax1d(M) ) :: Temp
real(sifm_real) :: Ur
!
Do d = 1, Ndim
   DO i = (M(d)-1)*Tau(d)+1, N-M(d)
      Ur = randf(IS)
      L = int( ( N - (i + M(d)) )*Ur + i + M(d))
      I1 = (/(i-(k-1)*Tau(d), k=1, M(d))/)
      I2 = (/(L-(k-1)*Tau(d), k=1, M(d))/)
      Temp( 1:M(d) )    = X( d, I1(1:M(d)) )
      X( d,I1(1:M(d)) ) = X( d, I2(1:M(d)) )
      X( d,I2(1:M(d)) ) = Temp( 1:M(d) )
   ENDDO
ENDDO
!
return
end subroutine BlockShuffleND2
!
subroutine BlockShuffle1D2(N,x,m,tau,IS)
integer, intent(in)  :: N,is
integer, intent(in) :: m,tau
real(sifm_real), intent(inout) :: x(n)
!
integer :: i, k, L
integer, dimension( M ) :: I1, I2
real(sifm_real), dimension( M ) :: Temp
real(sifm_real) :: Ur
!
DO i = (M-1)*Tau+1, N-M
   Ur = randf(IS)
   L = int( ( N - (i + M) )*Ur + i + M)
   I1 = (/(i-(k-1)*Tau, k=1, M)/)
   I2 = (/(L-(k-1)*Tau, k=1, M)/)
   Temp( 1:M )    = X( I1(1:M) )
   X( I1(1:M) ) = X( I2(1:M) )
   X( I2(1:M) ) = Temp( 1:M )
ENDDO
!
return
end subroutine BlockShuffle1D2
!
subroutine shuffleND(n,ndim,x,is)
integer, intent(in) :: n,ndim,is
character(len=1), intent(inout) :: x(ndim,n)
integer, dimension(n) :: indx
character(len=1) :: temp(ndim,n)
call rperm(1, N, indx, IS)
temp(:,1:n) = x(:,indx(1:n))
x = temp
return
end subroutine shuffleND
!
subroutine shuffleND2(n,ndim,x,is)
integer, intent(in) :: n,ndim,is
real(sifm_real), intent(inout) :: x(ndim,n)

real(sifm_real) :: temp(ndim,n)
integer, dimension(n) :: indx
call rperm(1, N, indx, IS)
temp(:,1:n) = x(:,indx(1:n))
x = temp
return
end subroutine shuffleND2
!
subroutine shuffle1D(n,x,is)
integer, intent(in) :: n,is
integer, intent(inout) :: x(n)
integer, dimension(n) :: indx
integer :: temp(n)
call rperm(1, N, indx, IS)
temp(1:n) = x(indx(1:n))
x = temp
return
end subroutine shuffle1D
!
subroutine shuffle1D2(n,x,is)
integer, intent(in) :: n,is
real(sifm_real), intent(inout) :: x(n)

real(sifm_real) :: temp(n)
integer, dimension(n) :: indx
call rperm(1, N, indx, IS)
temp(1:n) = x(indx(1:n))
x = temp
return
end subroutine shuffle1D2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION psi(xx) RESULT(fn_val)
IMPLICIT NONE
real(sifm_real), INTENT(IN) :: xx
real(sifm_real)             :: fn_val

real(sifm_real) :: dx0 = 1.461632144968362341262659542325721325_sifm_real
real(sifm_real) :: aug, den, piov4 = .785398163397448_sifm_real, sgn, upper,  &
             w, x, xmax1, xmx0, xsmall, z
INTEGER   :: i, m, n, nq
real(sifm_real) :: p1(7) = (/ .895385022981970D-02, .477762828042627D+01,  &
                        .142441585084029D+03, .118645200713425D+04,  &
                        .363351846806499D+04, .413810161269013D+04,  &
                        .130560269827897D+04 /),   &
                  q1(6) = (/ .448452573429826D+02, .520752771467162D+03,  &
                        .221000799247830D+04, .364127349079381D+04,  &
                        .190831076596300D+04, .691091682714533D-05 /)

real(sifm_real) :: p2(4) = (/ -.212940445131011D+01, -.701677227766759D+01,  &
                        -.448616543918019D+01, -.648157123766197D+00 /), &
                  q2(4) = (/  .322703493791143D+02,  .892920700481861D+02,  &
                         .546117738103215D+02,  .777788548522962D+01 /)

xmax1 = ipmpar(3)
xmax1 = MIN(xmax1, one/dpmpar(1))
xsmall = 1.d-9
x = xx
aug = rzero
IF (x >= half) GO TO 200
IF (ABS(x) > xsmall) GO TO 100
IF (x == rzero) GO TO 400
aug = -one / x
GO TO 150
100 w = - x
sgn = piov4
IF (w > rzero) GO TO 120
w = - w
sgn = -sgn
120 IF (w >= xmax1) GO TO 400
nq = INT(w)
w = w - nq
nq = INT(w*4.0_sifm_real)
w = 4.0_sifm_real * (w - nq * .25_sifm_real)
n = nq / 2
IF ((n+n) /= nq) w = one - w
z = piov4 * w
m = n / 2
IF ((m+m) /= n) sgn = - sgn
n = (nq + 1) / 2
m = n / 2
m = m + m
IF (m /= n) GO TO 140
IF (z == rzero) GO TO 400
aug = sgn * ((COS(z) / SIN(z)) * 4.0_sifm_real)
GO TO 150
140 aug = sgn * ((SIN(z) / COS(z)) * 4.0_sifm_real)
150 x = one - x
200 IF (x > 3.0_sifm_real) GO TO 300
den = x
upper = p1(1) * x
DO i = 1, 5
  den = (den + q1(i)) * x
  upper = (upper + p1(i+1)) * x
END DO

den = (upper + p1(7)) / (den + q1(6))
xmx0 = x - dx0
fn_val = den * xmx0 + aug
RETURN
300 IF (x >= xmax1) GO TO 350
w = one / (x * x)
den = w
upper = p2(1) * w

DO i = 1, 3
  den = (den + q2(i)) * w
  upper = (upper + p2(i+1)) * w
END DO

aug = upper / (den + q2(4)) - half / x + aug
350 fn_val = aug + LOG(x)
RETURN
400 fn_val = rzero
RETURN
END FUNCTION psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ipmpar (i) RESULT(fn_val)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER             :: fn_val

SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(one)
  CASE( 9)
    fn_val = MINEXPONENT(one)
  CASE(10)
    fn_val = MAXEXPONENT(one)
  CASE DEFAULT
    RETURN
END SELECT

RETURN
END FUNCTION ipmpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION dpmpar (i) RESULT(fn_val)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
double precision        :: fn_val
! Local variable
double precision :: lone
lone=1.0d0
SELECT CASE (i)
  CASE (1)
    fn_val = 1.D-20  !EPSILON(lone)
  CASE (2)
    fn_val = TINY(lone)
  CASE (3)
    fn_val = HUGE(lone)
END SELECT
RETURN
END FUNCTION dpmpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Function to add two strings  --> Slow routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine addStrings_slow(a,b)
character(len=*), intent(inout) :: a
character(len=*), intent(in)    :: b
!{Local variables}
integer :: n1, n2, n
n1 = len_trim(trim(adjustL(a)))
n2 = len_trim(trim(adjustL(b)))
n  = n1 + n2
a(n1+1:n) = b(1:n2)
return
end subroutine addStrings_slow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Add two Strings -- > Fast routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine addStrings(a,b)
character(len=*), intent(inout) :: a
character(len=*), intent(in)    :: b
a = trim( a ) // b  
return
end subroutine addStrings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE addStrings2(ST,STMAX,STLEN,ADST,ADLEN)
    INTEGER STMAX,STLEN,ADLEN,LIM
    CHARACTER(len=*) ST,ADST
    !
    LIM=ADLEN+STLEN
    IF(LIM > STMAX) THEN
       Stop '<ADDST> TRUNCATION HAS OCCOURRED'
       LIM=STMAX
    ENDIF
    IF(LIM <= STLEN) RETURN
    ST(STLEN+1:LIM)=ADST(1:LIM-STLEN)
    STLEN=LIM
    RETURN
END SUBROUTINE addStrings2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Functions for converting strings to integers and back
! %%DB written by Daniel Barr July, 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function toString(a) result(b)
implicit none
integer, intent(in) :: a
character(len=250) :: b
write(b,*) a
b=adjustl(b)
b=trim(b)
return
end function toString
!
function toInteger(a) result(b)
integer :: b
character(len=*) :: a
read(a,*) b
end function toInteger

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to convert lowercase ASCII characters in a string
! to uppercase. Generalizing this to multinational character sets
! requires a significant effort.
!
! This is just one of many possible implementations. For example, one
! might use a table indexed by the character value, this could be written
! as a function or as a subroutine with an optional output argument, etc.
!
subroutine UPCASE (STRING)
implicit none
character (len=*),INTENT(INOUT) :: STRING
integer I
character(len=1) C

do I = 1,LEN(STRING)
  C = STRING(I:I)
  if ((C >= 'a') .and. (C <= 'z')) then
    STRING(I:I) = CHAR(ICHAR(C)-32)
    end if
  end do
return
end subroutine UPCASE

! Subroutine to convert uppercase ASCII characters in a string
! to lowercase. Generalizing this to multinational character sets
! requires a significant effort.
!
! This is just one of many possible implementations. For example, one
! might use a table indexed by the character value, this could be written
! as a function or as a subroutine with an optional output argument, etc.
!
subroutine DOWNCASE (STRING)
implicit none
character (len=*),INTENT(INOUT) :: STRING
integer I
character(len=1) C

do I = 1,LEN(STRING)
  C = STRING(I:I)
  if ((C >= 'A') .and. (C <= 'Z')) then
    STRING(I:I) = CHAR(ICHAR(C)+32)
    end if
  end do
return
end subroutine DOWNCASE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extract the mean value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function mean(X,n) result(r)
integer, intent(in) :: n
real(sifm_real), intent(in) :: x(n)
real(sifm_real) :: r
!
r = sum( X(1:n) )/real(n,sifm_real)
!
return
end function mean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the minimum and maximum of an array of real and/or integer numbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rmin1d(X)
real(sifm_real), intent(in) :: x(:)
real(sifm_real) :: rmin1d
!{Local variable}
integer :: i
rmin1d=x(1)
do i=2, size(x,1)
   rmin1d = min(rmin1d, x(i))
end do
return
end function rmin1d
!
function rmax1d(X)
real(sifm_real), intent(in) :: x(:)
real(sifm_real) :: rmax1d
!{Local variable}
integer :: i
rmax1d=x(1)
do i=2, size(x,1)
   rmax1d = max(rmax1d, x(i))
end do
return
end function rmax1d
!
pure function imax1d(X) result(r)
integer, intent(in) :: x(:)
integer :: r
!{Local variable}
integer :: i
r=x(1)
do i=2, size(x,1)
   r = max(r, x(i))
end do
return
end function imax1d
!
function imin1d(X) result(r)
integer, intent(in) :: x(:)
integer :: r
!{Local variable}
integer :: i
r=x(1)
do i=2, size(x,1)
   r = min(r, x(i))
end do
return
end function imin1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sort(n,arr,indx)
integer, intent(in) :: n
real(sifm_real), dimension(n), intent(in) :: arr
integer, dimension(n), intent(out) :: indx
!
integer :: i, temp
logical :: swap
real(sifm_real) :: pt
real(sifm_real), dimension(n) :: v
!
FORALL (i=1:N)
   indx(i) = i
END FORALL
v = arr
swap = .true.
do while (swap)
   swap = .false.
   do i = 1, n-1
      if ( v(i) > v(i+1) ) then
           swap = .true.
           temp = v(i); v(i) = v(i+1); v(i+1) = temp
           pt = indx(i); indx(i) = indx(i+1); indx(i+1) = pt
      endif
   enddo
enddo
return
end subroutine sort
!
subroutine sort2(n,r0,indx0,indx)
integer, intent(in) :: n
integer, intent(in) :: indx0(n)
real(sifm_real), intent(in) :: r0(n)
integer, intent(out) :: indx(n)
integer :: i, temp, j
real(sifm_real) :: rt
real(sifm_real), dimension(n) :: r
r=r0
indx=indx0
do i = 2, n
   rt=r(i); temp=indx(i)
   IF (rt >= r(i-1)) Cycle
   r(i)=r(i-1); indx(i)=indx(i-1)
   do j=i-2, 1, -1
      IF ( rt >= r(j) ) Exit
  r(j+1) = r(j); indx(j+1)=indx(j) 
   enddo
   r(j+1)=rt; indx(j+1)=temp
enddo
return
end subroutine sort2
!
function vcopy(Xin,istart,M,istep,Sign) result(Xou)
       integer, intent(in) :: istart, M, istep, Sign
       real(sifm_real), dimension(:), intent(in) :: Xin
       real(sifm_real), dimension(M) :: Xou
       integer :: k
       FORALL (k=1:M)
          Xou(k) = Xin(Istart+Sign*(k-1)*Istep)
       END FORALL
       return
end function vcopy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute covariance matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function getCovariance(Nframes,X, Y) result(Cov)
   integer, intent(in) :: Nframes
   real(sifm_real), intent(in) :: X(Nframes), Y(Nframes)
   real(sifm_real) :: cov
!
   real(sifm_real) :: norm, sxy
!
   norm = one / real(nframes, sifm_real)
   cov = rzero
   sxy = dot_product( x - mean(X, Nframes), y - mean(Y, Nframes) )
   cov = Sxy * norm
   return
end Function getCovariance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function find(i,list) result(r)
integer, intent(in) :: i
integer, dimension(:), intent(in) :: list
logical :: r
integer :: k
r = .FALSE.
DO k=1, size(list,1)
   IF ( list(k) == i) THEN
        r = .TRUE.
        return
   ENDIF
ENDDO
return
end function find
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform a statistical test
! -- Statistical test of significance (Level = 95% is chosen for statProb=0.95)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine statTest(N, r0, val_in, val_out, stat_prob)
implicit none
integer, intent(in) :: N
real(sifm_real), intent(in)  :: r0, val_in
real(sifm_real), intent(out) :: val_out
real(sifm_real), intent(in)  :: stat_prob
!
! Local variables
real(sifm_real) :: ksi, prob, z, u
!
ksi  = half * log( (one + r0) / (one - r0) )
z    = half * log( (one + val_in) / (one - val_in) )
u    = (z - ksi) * sqrt( real(N, sifm_real) )
prob = erfcc( abs(u) / SQRT2 )
IF (prob > stat_prob) THEN
    val_out = rzero
ELSE
    val_out = val_in
ENDIF	
!  
return
end subroutine statTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION erfcc_s(x)
IMPLICIT NONE
REAL(sifm_real), INTENT(IN) :: x
REAL(sifm_real) :: erfcc_s
REAL(sifm_real) :: t,z
REAL(sifm_real), DIMENSION(10) :: coef = (/-1.26551223_sifm_real,1.00002368_sifm_real,&
  0.37409196_sifm_real,0.09678418_sifm_real,-0.18628806_sifm_real,0.27886807_sifm_real,&
  -1.13520398_sifm_real,1.48851587_sifm_real,-0.82215223_sifm_real,0.17087277_sifm_real/)
z=abs(x)
t=one/(one+half*z)
erfcc_s=t*exp(-z*z+poly(t,coef))
if (x < 0.0) erfcc_s=2.0_sifm_real-erfcc_s
END FUNCTION erfcc_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION erfcc_v(x)
IMPLICIT NONE
REAL(sifm_real), DIMENSION(:), INTENT(IN) :: x
REAL(sifm_real), DIMENSION(size(x)) :: erfcc_v,t,z
REAL(sifm_real), DIMENSION(10) :: coef = (/-1.26551223_sifm_real,1.00002368_sifm_real,&
  0.37409196_sifm_real,0.09678418_sifm_real,-0.18628806_sifm_real,0.27886807_sifm_real,&
  -1.13520398_sifm_real,1.48851587_sifm_real,-0.82215223_sifm_real,0.17087277_sifm_real/)
z=abs(x)
t=one/(one+half*z)
erfcc_v=t*exp(-z*z+poly(t,coef))
where (x < 0.0) erfcc_v=2.0_sifm_real-erfcc_v
END FUNCTION erfcc_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION poly_s(x,coeffs)
 REAL(sifm_real), INTENT(IN) :: x
 REAL(sifm_real), DIMENSION(:), INTENT(IN) :: coeffs
 REAL(sifm_real) :: poly_s
!
 REAL(sifm_real) :: pow
 REAL(sifm_real), DIMENSION(:), ALLOCATABLE :: vec
 INTEGER :: i,n,nn
!
 n=size(coeffs)
 if (n <= 0) then
  poly_s=rzero
 else if (n < NPAR_POLY) then
  poly_s=coeffs(n)
  do i=n-1,1,-1
   poly_s=x*poly_s+coeffs(i)
  end do
 else
  allocate(vec(n+1))
  pow=x
  vec(1:n)=coeffs
  do
   vec(n+1)=rzero
   nn=ishft(n+1,-1)
   vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
   if (nn == 1) exit
   pow=pow*pow
   n=nn
  end do
  poly_s=vec(1)
  deallocate(vec)
 end if
 return
END FUNCTION poly_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION poly_v(x,coeffs)
 REAL(sifm_real), DIMENSION(:), INTENT(IN) :: coeffs,x
 REAL(sifm_real), DIMENSION(size(x)) :: poly_v
 INTEGER :: i,n,m
 m=size(coeffs)
 n=size(x)
 if (m <= 0) then
  poly_v=rzero
 else if (m < n .or. m < NPAR_POLY) then
  poly_v=coeffs(m)
  do i=m-1,1,-1
   poly_v=x*poly_v+coeffs(i)
  end do
 else
  do i=1,n
   poly_v(i)=poly_s(x(i),coeffs)
  end do
 end if
 return
END FUNCTION poly_v
!
function square(x) result(r)
    implicit none
    real(sifm_real), intent(in) :: x
    real(sifm_real) :: r
    r = x * x
    return
end function square
!
Subroutine TransferFunctionInn_2d(z, r, ifunc)
    implicit none
    real(sifm_real), intent(in) :: Z(:,:)
    integer, intent(in) :: iFunc
    real(sifm_real), intent(out) :: R(:,:)
!
    integer :: i, j
!
    do i=1, size(Z,1)
       do j=1, size(Z, 2)
          R(i,j) = TransferFunctionInn_0d(Z(i,j), ifunc)
       enddo
    enddo
return
end subroutine TransferFunctionInn_2d
!
subroutine TransferFunctionInn_1d(Z, R, ifunc)
    implicit none
    real(sifm_real), intent(in) :: Z(:)
    integer, intent(in) :: iFunc
    real(sifm_real), intent(out) :: R(:)
!
    integer :: i
!
    do i=1, size(Z)
       R(i) = TransferFunctionInn_0d(Z(i), ifunc)
    enddo
return
end subroutine TransferFunctionInn_1d

function TransferFunctionInn_0d(z, ifunc) result(r)
    implicit none
    real(sifm_real), intent(in) :: Z
    integer, intent(in) :: iFunc
    real(sifm_real) :: R
!
    if (ifunc == 1) THEN
        R = 1.0_sifm_real/(1.0_sifm_real + exp(-Z))
    else if( ifunc == 2 ) THEN
        R = 1.0_sifm_real/(1.0_sifm_real + square(Z))
    else
        R = tanh( Z )
    endif
return
end function TransferFunctionInn_0d
! 
Subroutine TransferFunctionInn_derivative_2d(z, r, ifunc)
    implicit none
    real(sifm_real), intent(in) :: Z(:,:)
    integer, intent(in) :: iFunc
    real(sifm_real), intent(out) :: R(:,:)
!
    integer :: i, j
!
    do i=1, size(Z,1)
       do j=1, size(Z, 2)
          R(i,j) = TransferFunctionInn_derivative_0d(Z(i,j), ifunc)
       enddo
    enddo
return
end subroutine TransferFunctionInn_derivative_2d
!
subroutine TransferFunctionInn_derivative_1d(Z, R, ifunc)
    implicit none
    real(sifm_real), intent(in) :: Z(:)
    integer, intent(in) :: iFunc
    real(sifm_real), intent(out) :: R(:)
!
    integer :: i
!
    do i=1, size(Z)
       R(i) = TransferFunctionInn_derivative_0d(Z(i), ifunc)
    enddo
return
end subroutine TransferFunctionInn_derivative_1d

function TransferFunctionInn_derivative_0d(z, ifunc) result(r)
    implicit none
    real(sifm_real), intent(in) :: Z
    integer, intent(in) :: iFunc
    real(sifm_real) :: R
!
    if (ifunc == 1) THEN
        R = TransferFunctionInn_0d(Z, ifunc) * ( 1.0_sifm_real - TransferFunctionInn_0d(Z, ifunc) )
    else if( ifunc == 2 ) THEN
        R = -2.0_sifm_real * Z * square( TransferFunctionInn_0d(z,ifunc) )
    else
        R = 1.0_sifm_real - square( TransferFunctionInn_0d(z,ifunc) )
    endif
return
end function TransferFunctionInn_derivative_0d
!
subroutine getTimeLag(A_n, Tau)
real(sifm_real), intent(in)  :: A_n(:)
Integer, intent(out) :: tau

!{Local variables}
integer         :: time, t, incr, N
real(sifm_real) :: C, mu_A, sigma2


!{set intial value}
N = size(A_n)

!{compute mean and variances}
mu_A = sum( A_n ) / real(N, sifm_real )
sigma2 = sum( (A_n - mu_A)*(A_n - mu_A) ) / real(N-1,sifm_real)

!{Compute time lag}
t = 1
incr = 1
Loop_TIME: do time = 1, N
   IF (t < N - 1) THEN
       C = sum( A_n(1:(N-t))*A_n((1+t):N) ) / real((N-t), sifm_real)
       C = (C - mu_A*mu_A) / sigma2
       if (C <= rzero) exit
       t = t + incr
       incr = incr + 1
   ENDIF
end do Loop_TIME
Tau = time
if (Tau < 1) Tau = 1
Write(*,*) "Time Lag =   ", Tau
return
end subroutine getTimeLag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the embedded dimension parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_embdparam(Natoms, Ndim, Mopt, Tau)
Implicit None
!
Integer, intent(in) :: Natoms, Ndim
integer, intent(in) :: Mopt(:), Tau(:)
!  
  integer :: i, j, d, nunit
  integer, dimension (Ndim) :: T, M
!
  nunit = new_unit()
  open(unit=nunit, file='embd.txt', status='unknown', action='write')
!  
  j=0
  DO i = 1, Natoms
     DO d = 1, Ndim
        J = J + 1
        M(d) = Mopt(J)
        T(d) = Tau(J)
     ENDDO
     write(nunit,*) (t(d), d=1, Ndim), (m(d), d=1, Ndim)
  ENDDO
  CLOSE(nunit)
!
return
!
end subroutine Print_embdparam
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the embedded dimension parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_embdparam(Natoms, Ndim, Mopt, Topt)
Implicit None
!
Integer, intent(in) :: Natoms, Ndim
integer, intent(out) :: Mopt(:,:), Topt(:,:)
!
  integer :: i, d, nunit
!
  nunit = new_unit()
  open(unit=nunit, file='embd.txt', status='old', action='read')
!  
  DO i = 1, Natoms
     read(nunit,*) (topt(i,d), d=1, Ndim), (mopt(i,d), d=1, Ndim)
  ENDDO
  CLOSE(nunit)
!
return
!
end subroutine read_embdparam
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the trajectories
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_xyz(Nframes, Natoms, Ndim, XYZ)
Implicit None
integer, intent(in) :: Nframes, Natoms, Ndim 
real(sifm_real), intent(in) :: XYZ(:,:,:)
!  
integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='traj.xyz', status='unknown', action='write')
write(nunit, '("SIFM: XYZ Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   DO i = 1, Natoms
      write(nunit, '(A7, 3F12.6)') "FRAME  ", (xyz(i,d,t), d=1, Ndim) 
   ENDDO
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
!
return
!
end subroutine write_xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the symbolic Trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writei_xyzs(Nframes, Natoms, Ndim, XYZs)
Implicit None
integer, intent(in) :: Nframes, Natoms, Ndim 
integer, intent(in) :: XYZs(:,:,:)
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='trajs.xyz', status='unknown', action='write')
write(nunit, '("SIFM: XYZS Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   DO i = 1, Natoms
      write(nunit, '(A7,3I3)') "FRAME  ", (xyzs(i,d,t), d=1, Ndim) 
   ENDDO
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
return
end subroutine writei_xyzs
!
subroutine writec_xyzs(Nframes, Natoms, Ndim, XYZs)
Implicit None
integer, intent(in) :: Nframes, Natoms, Ndim 
character(len=1), intent(in) :: XYZs(:,:,:)
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='trajs.xyz', status='unknown', action='write')
write(nunit, '("SIFM: XYZS Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   DO i = 1, Natoms
      write(nunit, '(A7,3A)') "FRAME  ", (xyzs(i,d,t), d=1, Ndim) 
   ENDDO
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
return
end subroutine writec_xyzs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the trajectories in CSV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_csv(Nframes, Natoms, Ndim, XYZ)
Implicit None
integer, intent(in) :: Nframes, Natoms, Ndim 
real(sifm_real), intent(in) :: XYZ(:,:,:)
!  
integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='traj.csv', status='unknown', action='write')
write(nunit, '("SIFM: CSV Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   write(nunit, *) "FRAME  ", ",", ((xyz(i,d,t),",", d=1, Ndim),i=1,Natoms) 
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
!
return
!
end subroutine write_csv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the symbolic trajectories in CSV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writei_csvs(Nframes, Natoms, Ndim, XYZs)
Implicit None
integer, intent(in) :: Nframes, Natoms, Ndim 
integer, intent(in) :: XYZs(:,:,:)
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='trajs.csv', status='unknown', action='write')
write(nunit, '("SIFM: CSV Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   write(nunit, *) "FRAME  ", ",", ((xyzs(i,d,t),",", d=1, Ndim),i=1,Natoms) 
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
return
end subroutine writei_csvs
!
subroutine writec_csvs(Nframes, Natoms, Ndim, XYZs)
Implicit None
integer, intent(in) :: Nframes, Natoms, Ndim 
character(len=1), intent(in) :: XYZs(:,:,:)
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='trajs.csv', status='unknown', action='write')
write(nunit, '("SIFM: CSV Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   write(nunit, *) "FRAME  ", ",", ((xyzs(i,d,t),",", d=1, Ndim),i=1,Natoms) 
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
return
end subroutine writec_csvs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the trajectories in XYZ format 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xyz(Nframes, Natoms, Ndim, X)
Implicit None
!
integer, intent(in) :: Nframes, Natoms, Ndim 
real(sifm_real), intent(out) :: X(:,:)
!  
  integer :: i, d, j, t, nunit
  character(len=7) :: remark
  character(len=80) :: line
  real(sifm_real), dimension(Ndim) :: xyz
  
!
nunit = new_unit()
open(unit=nunit, file='traj.xyz', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   j = 0
   DO i = 1, Natoms
      read(nunit,'(A80)') line
      IF (line(1:5) == 'FRAME') THEN
          read(line,'(A7, 3F12.6)') remark, ( XYZ(d), d = 1, Ndim ) 
          DO d = 1, Ndim
             j = j + 1
             X(j,t) = xyz(d)
          ENDDO
      ENDIF
   ENDDO
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit)
return
end subroutine read_xyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readc_xyzs(Nframes, Natoms, Ndim, XS)
Implicit None
!
integer, intent(in) :: Nframes, Natoms, Ndim 
character(len=1), intent(out) :: XS(:,:)
!  
integer :: i, d, j, t, nunit
character(len=7) :: remark
character(len=80) :: line
character(len=1), dimension(Ndim) :: xyzs  
!
nunit = new_unit()
open(unit=nunit, file='trajs.xyz', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   j = 0
   DO i = 1, Natoms
      read(nunit,'(A80)') line
      IF (line(1:5) == 'FRAME') THEN
          read(line,'(A7, 3A)') remark, ( XYZS(d), d = 1, Ndim ) 
          DO d = 1, Ndim
             j = j + 1
             XS(j,t) = xyzs(d)
          ENDDO
      ENDIF
   ENDDO
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit)
return
end subroutine readc_xyzs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readi_xyzs(Nframes, Natoms, Ndim, XS)
Implicit None
!
integer, intent(in) :: Nframes, Natoms, Ndim 
integer, intent(out) :: XS(:,:)
!  
integer :: i, d, j, t, nunit
character(len=7) :: remark
character(len=80) :: line
integer, dimension(Ndim) :: xyzs  
!
nunit = new_unit()
open(unit=nunit, file='trajs.xyz', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   j = 0
   DO i = 1, Natoms
      read(nunit,'(A80)') line
      IF (line(1:5) == 'FRAME') THEN
          read(line,'(A7, 3I3)') remark, ( XYZS(d), d = 1, Ndim ) 
          DO d = 1, Ndim
             j = j + 1
             XS(j,t) = xyzs(d)
          ENDDO
      ENDIF
   ENDDO
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit)
return
end subroutine readi_xyzs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readc_csvs(Nframes, Natoms, Ndim, XS)
Implicit None
!
integer, intent(in) :: Nframes, Natoms, Ndim 
character(len=1), intent(out) :: XS(:,:)
!  
  integer           :: i, d, j, t, nunit, ierr
  character(len=7)  :: remark
  character(len=80) :: line
  character(len=1)  :: comma 
  character(len=1), dimension(Natoms,Ndim,Nframes) :: xyzs
  character(len=80), allocatable :: Fields(:)
  integer :: Nfields
  
!
nunit = new_unit()
open(unit=nunit, file='trajs.csv', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   read(nunit,'(A80)') Line
   CALL lenFields(Line, ',', Nfields)
   Allocate( Fields(Nfields), stat = ierr)
   IF (ierr /= 0) Stop 'Stop: Error allocating fields'
   CALL Split(Line, ',', Fields) 
   j = 0
   DO i = 1, Natoms
      DO d = 1, Ndim
         j = j + 1
         read(fields(j+1), *) XS(j,t)
		 write(*,*) XS(j,t)
      ENDDO
   ENDDO
   deallocate( Fields, stat = ierr)
   IF (ierr /= 0) Stop 'Stop: Error deallocating fields'
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit) 
return
end subroutine readc_csvs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readi_csvs(Nframes, Natoms, Ndim, XS)
Implicit None
!
integer, intent(in) :: Nframes, Natoms, Ndim 
integer, intent(out) :: XS(:,:)
!  
  integer           :: i, d, j, t, nunit, ierr
  character(len=7)  :: remark
  character(len=80) :: line
  character(len=1)  :: comma 
  integer, dimension(Natoms,Ndim,Nframes) :: xyzs
  character(len=80), allocatable :: Fields(:)
  integer :: Nfields
  
!
nunit = new_unit()
open(unit=nunit, file='trajs.csv', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   read(nunit,'(A80)') Line
   CALL lenFields(Line, ',', Nfields)
   Allocate( Fields(Nfields), stat = ierr)
   IF (ierr /= 0) Stop 'Stop: Error allocating fields'
   CALL Split(Line, ',', Fields) 
   j = 0
   DO i = 1, Natoms
      DO d = 1, Ndim
         j = j + 1
         read(fields(j+1), *) XS(j,t)
		 write(*,*) XS(j,t)
      ENDDO
   ENDDO
   deallocate( Fields, stat = ierr)
   IF (ierr /= 0) Stop 'Stop: Error deallocating fields'
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit) 
return
end subroutine readi_csvs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the trajectories in CSV format 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_csv(Nframes, Natoms, Ndim, X)
Implicit None
!
integer, intent(in) :: Nframes, Natoms, Ndim 
real(sifm_real), intent(out) :: X(:,:)
!  
  integer           :: i, d, j, t, nunit, ierr
  character(len=7)  :: remark
  character(len=80) :: line
  character(len=1)  :: comma 
  real(sifm_real), dimension(Natoms,Ndim,Nframes) :: xyz
  character(len=80), allocatable :: Fields(:)
  integer :: Nfields
  
!
nunit = new_unit()
open(unit=nunit, file='traj.csv', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   read(nunit,'(A80)') Line
   CALL lenFields(Line, ',', Nfields)
   Allocate( Fields(Nfields), stat = ierr)
   IF (ierr /= 0) Stop 'Stop: Error allocating fields'
   CALL Split(Line, ',', Fields) 
   j = 0
   DO i = 1, Natoms
      DO d = 1, Ndim
         j = j + 1
         read(fields(j+1), *) X(j,t)
		 write(*,*) X(j,t)
      ENDDO
   ENDDO
   deallocate( Fields, stat = ierr)
   IF (ierr /= 0) Stop 'Stop: Error deallocating fields'
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit) 
return
end subroutine read_csv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Split(Line, Symbol, Fields)
implicit none 
character(len=*), intent(inout) :: Line 
character(len=1), intent(in) :: Symbol
character(len=80), intent(out) :: Fields(:)
!
integer :: N, siz, i
character(len=80) :: var
!
Line = adjustL( Line(1:len_trim(Line)) )
siz = len_trim( Line )
var = ''
N = 1
DO i=1, siz 
   IF (Line(i:i) /= symbol) THEN
       CALL addStrings(var, Line(i:i)) 
   ELSE
	   Fields(N) = var
	   var = ''
       N = N + 1
   ENDIF 
ENDDO
return
END SUBROUTINE Split
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
Subroutine lenFields(Line, Symbol,  N)
implicit none 
character(len=*), intent(inout) :: Line 
character(len=1), intent(in) :: Symbol
integer, intent(out) :: N
!
integer :: siz, i
!

Line = adjustL( Line(1:len_trim(Line)) )
siz = len_trim( Line )
N = 0
DO i=1, siz 
   IF (Line(i:i) == symbol) THEN
       N = N + 1
   ENDIF 
ENDDO
return
END SUBROUTINE lenFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WRITE_LTE(Natoms, Nframes, LocalTES, LocalTE, qTEShuffle)
!
Integer, intent(in) :: qTESHUFFLE, Natoms, Nframes
real(sifm_real), intent(in) :: LocalTES(:,:,:), LocalTE(:,:,:)
!
integer             :: t, i, j, nunit
! 
nunit = new_unit()
open(unit=nunit, file='localte.csv', status='unknown', action='write')
write(nunit, *) "Frame",",",( (j,",", j=1, Natoms), i=1, Natoms) 
DO t=1, Nframes
   write(nunit,*) t, ",", ((LocalTE(t,i,j),",", j=1, natoms), i=1, natoms)
ENDDO
!
IF (qteshuffle > 0) THEN
nunit = new_unit()
open(unit=nunit, file='localste.csv', status='unknown', action='write')
write(nunit, *) "Frame",",",( (j,",", j=1, Natoms), i=1, Natoms) 
DO t=1, Nframes
   write(nunit,*) t,",", ((LocalTES(t,i,j),",", j=1, natoms), i=1, natoms)
ENDDO
ENDIF
!
return
!
end subroutine WRITE_LTE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the transfer entropy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine write_TE(Natoms, TES, TE, qTEShuffle)
! Implicit None
! integer, intent(in) :: qTeShuffle, Natoms
! real(sifm_real), intent(in) :: TES(:,:), TE(:,:)
! !  
! integer             :: i, j, nunit
! !  
! nunit = new_unit()
! open(unit=nunit, file='te.csv', status='unknown', action='write')
! write(nunit, *) "Index",",",( j,",", j=1, Natoms )
! DO i=1, Natoms
   ! write(nunit,*) i, ",", (TE(i,j),",", j=1, natoms)
! ENDDO
! !
! IF (qteshuffle > 0) THEN
! nunit = new_unit()
! open(unit=nunit, file='ste.csv', status='unknown', action='write')
! write(nunit, *) "Index",",",( j,",", j=1, Natoms )
! DO i=1, Natoms
   ! write(nunit,*) i, ",",  (TES(i,j), ",", j=1, natoms)
! ENDDO
! ENDIF
! !
! return
! !
! end subroutine write_te
! !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the transfer entropy when running multiple processors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_TE(Natoms, TES, TE, qTeShuffle, Ndim, model)
Implicit None
integer, intent(in) :: Natoms, qTeShuffle
real(sifm_real), intent(in) :: TES(:,:), TE(:,:)
integer, optional :: model, Ndim 
!  
integer             :: i, j, nunit
!  
nunit = new_unit()
open(unit=nunit, file='te.csv', status='unknown', action='write')
IF ( (present(model)) .and. (model == 3) ) THEN 
write(nunit, *) "Index",",",( j,",", j=1, Ndim )
DO i=1, Natoms
   write(nunit,*) i, ",", (TE(i,j), j=1, ndim)
ENDDO
ELSE 
write(nunit, *) "Index",",",( j,",", j=1, Natoms )
DO i=1, Natoms
   write(nunit,*) i, ",", (TE(i,j), j=1, natoms)
ENDDO
ENDIF 
!
IF (qteshuffle > 0) THEN
nunit = new_unit()
open(unit=nunit, file='ste.csv', status='unknown', action='write')
IF ( (present(model)) .and. (model == 3) ) THEN 
write(nunit, *) "Index",",",( j,",", j=1, Ndim )
DO i=1, Natoms
   write(nunit,*) i, ",", (TES(i,j), j=1, ndim)
ENDDO
ELSE 
write(nunit, *) "Index",",",( j,",", j=1, Natoms )
DO i=1, Natoms
   write(nunit,*) i, ",", (TES(i,j), j=1, natoms)
ENDDO
ENDIF 
ENDIF
return
end subroutine write_te
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the Mutual Information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_MI(Natoms, MIS, MI, qMIShuffle)
  implicit none
!
  integer, intent(in) :: qMIShuffle, Natoms 
  real(sifm_real), intent(in) :: MI(:,:), MIS(:,:)
!  
  integer             :: i, j, nunit
!  
  nunit = new_unit()
  open(unit=nunit, file = 'mi.csv', status='unknown', action='write')
  write(nunit, *) "Index",",",( j,",", j=1, Natoms )
  DO i=1, Natoms
     write(nunit,*) i, ",", (MI(i,j),",", j=1, Natoms)
  ENDDO
  close(nunit)
!
  IF (qMIShuffle > 0) THEN
  nunit = new_unit()
  open(unit=nunit, file = 'mis.csv', status='unknown', action='write')
  write(nunit, *) "Index",",",( j,",", j=1, Natoms )
  DO i=1, Natoms
     write(nunit,*) i, ",", ( MIS(i,j), ",", j=1, Natoms )
  ENDDO
  close(nunit)
  ENDIF
!
return
end subroutine write_MI
!
END MODULE TEUTILS_CLASS
