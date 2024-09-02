module matinv_module

use iso_fortran_env, only: nwrite => output_unit

implicit none

integer, parameter :: dp=kind(0.d0)
integer :: i, j, k

contains 

subroutine matinv(n,a,x,ier,idgt)
! This is a lightweight square matrix inversion routine, using the LU-decomposition,
! written specifically for uncert, my parameter uncertainties routine.  
! Error handling is rudimentary: 
! - Checks for zeros on diagonal and exits if there are any 
!    (there shouldn't be any for the covariance matrix)
! - Checks that the product A*X = I to within 10^-(idgt)
!
!  n = dimension of square matrix a(n,n)
!  a(n,n) = input matrix
!  x(n,n) = output matrix (inverse of a)
!  ier = error signal to calling program
!  idgt = requested number of good digits in inverted matrix (i.e., precision is < 1E-idgt)

integer, intent(in)  :: idgt, n
real(dp), intent(in), dimension(n,n) :: a
real(dp), intent(out), dimension(n,n) :: x
integer, intent(out) :: ier 
real(dp), dimension(n,n) :: l, u, c, y

! write(nwrite, '(5x,a)') 'A(i,j) ='
! do i = 1, n
!   write(nwrite, '(5(1x,f9.4))') (a(i,j),j=1,n)
! enddo

! Check for zeros on the diagnals
ier = 0
do i = 1,n
 if (a(i,i) == 0.0) then
  ier = 1
  return
 endif
enddo
    
! Carry out the LU decomposition of square matrix A
do i = 1,n
 do j = 1,n
  if (i <= j) then
   u(i,j) = a(i,j)
   do k = 1,i-1
    u(i,j) = u(i,j) - l(i,k) * u(k,j)
   enddo
  else if (i > j) then
   u(i,j) = 0.0
  end if
  if (i >= j) then
   l(i,j) = a(i,j) / u(j,j)
   do k = 1,j-1
    l(i,j) = l(i,j) - (l(i,k) * u(k,j))/u(j,j)
   enddo
  else if (i < j) then
   l(i,j) = 0.0
  end if
 enddo
enddo

! confirm that L*U = A
!write(nwrite, '(5x,a)') 'L(i,j) ='
!do i = 1, n
!  write(nwrite, '(5(1x,f9.4))') (l(i,j),j=1,n)
!enddo

!write(nwrite, '(5x,a)') 'U(i,j) ='
!do i = 1, n
!  write(nwrite, '(5(1x,f9.4))') (u(i,j),j=1,n)
!enddo

call matrix_mult(n,l,u,c)

!write(nwrite, '(5x,a)') 'L*U ='
!do i = 1, n
!  write(nwrite, '(5(1x,f9.4))') (c(i,j),j=1,n)
!enddo

! calculate the Y matrix which solves L*Y = I
 do i = 1,n
  do j = 1,n
   if (i > j) then
    y(i,j) = 0.0
    do k = 1,i-1
     y(i,j) = y(i,j) - l(i,k)*y(k,j)
    enddo
   else if (i == j) then
    y(i,j) = 1.0
    do k = 1,i-1
     y(i,j) = y(i,j) - l(i,k)*y(k,j)
    enddo
   else if (i < j) then
    y(i,j) = 0.0
   end if
  enddo
 enddo

!write(nwrite, '(5x,a)') 'Y ='
!do i = 1, n
!  write(nwrite, '(5(1x,f9.4))') (y(i,j),j=1,n)
!enddo

! confirm that L*Y = I
call matrix_mult(n,l,y,c)
!write(nwrite, '(5x,a)') 'L*Y ='
!do i = 1, n
!  write(nwrite, '(5(1x,f9.4))') (c(i,j),j=1,n)
!enddo

! calculate the X matrix which solves U*X = Y
 do i = n,1,-1
  do j = 1,n
   x(i,j) = y(i,j)/u(i,i)
   do k=i+1,n
    x(i,j) = x(i,j) - u(i,k)*x(k,j) / u(i,i)
   enddo
  enddo
 enddo

! write(nwrite, '(5x,a)') 'X ='
! do i = 1, n
!   write(nwrite, '(5(1x,f9.4))') (x(i,j),j=1,n)
! enddo

! confirm that A*X = I
 call matrix_mult(n,a,x,c)
!write(nwrite, '(5x,a)') 'A*X ='
do i = 1, n
 do j = 1, n
!  write(nwrite, '(5(1x,f9.4))') (c(i,j),j=1,n)
  if (i == j) then
   if (dabs(c(i,i)-1.0d0) > 10**dfloat(-idgt)) ier = 2   ! ier = 2 for precision error
  else
   if (c(i,j) > 10**dfloat(-idgt)) ier = 2
  endif
 enddo
enddo
end subroutine matinv

subroutine matrix_mult(n,a,b,c)
integer :: i, j, k
integer, intent(in)  :: n
real(8), intent(in), dimension(n,n) :: a
real(8), intent(in), dimension(n,n) :: b
real(8), intent(out), dimension(n,n) :: c

do i = 1,n
 do j = 1,n
  c(i,j) = 0.0
  do k = 1,n
   c(i,j) = c(i,j) + a(i,k) * b(k,j)
  enddo
 enddo
enddo
end subroutine matrix_mult

end module matinv_module

