! Matthew Quenneville
! Phys 395 Assignment 1
!
! September 11, 2014
!
! Approximates a function using a series of Chebyshev polynomials. 
! Creates a matrix indexed by polynomial order and x position.
! Finds coefficients as solution to matrix using Gauss-Jordan Elimination.

program poly
implicit none

! 'n' is number of grid points in x
integer, parameter :: n = 100, m = n-1
integer i, j
real*16 b(0:m,n), c(0:m), y(n), f(n), x(n)
!real*16 test(3,3), test_b(3), test_x(3)
!test(1,:)=(/1,3,2/)
!test(2,:)=(/1,0,2/)
!test(3,:)=(/2,4,0/)
!test_b=(/13,7,10/)
!call gjsolve(test,test_b,test_x)
!write (*,*) test_x
do i = 1,n
   
    ! Approximation grid
    x(i) = 2.0*(i-1.0)/(n-1.0) - 1.0
    f(i) = 1.0/(1.0 + x(i)**2/10.0)
    
    ! Chebyshev polynomial basis
    forall (j=0:m) b(j,i) = cos(j*acos(x(i)))

end do

!Solve for coefficients and calculate approximate function values.
call gjsolve(b,f,c)
y=matmul(b,c)

do i=1,n
   write (*,'(8g20.12)') x(i), f(i), y(i), y(i)-f(i)
end do


! ... figure out how to call gjsolve to find coefficients ...

! ... output the approximation on uniform grid of 1000 points ...

contains
  !Swaps row numbers 'row1' and 'row2' in matrix 'M'
  subroutine swap(M,row1,row2)
    real*16 M(:,:),N(size(M,2))
    integer row1,row2
    N=M(row1,:)
    M(row1,:)=M(row2,:)
    M(row2,:)=N
    !write (*,*) 'Swap: ', row1, row2
  end subroutine swap
  
  !Scales row number 'row' in 'M' by real number 'a'
  subroutine scale(M,row,a)
    real*16 M(:,:),a
    integer row,j
    forall (j=1:size(M,2)) M(row,j)=a*M(row,j)
  end subroutine scale
  
  !Adds row 'row2' scaled by 'fact' to 'row1' in matrix 'M' 
  subroutine add(M,row1,row2,factor)
    real*16 M(:,:),factor
    integer row1,row2,j
    !write (*,*) 'Add: ', row1, row2, factor
    forall (j=1:size(M,2)) M(row1,j) = M(row1,j)+factor*M(row2,j)
  end subroutine add
  
  ! Solve A.x = B by Gauss-Jordan elimination
  subroutine gjsolve(A, B, x)

    real*16 A(:,:), B(:), x(:)
    integer nRows,i,j,maxRow
    real*16 N(size(B),size(x)+1),fact

    nRows=size(B)
    
    !Check if matrix is square
    if (size(B)/=size(x)) then
       write (*,*) "Error, matrix not square"
       return
    end if
    
    !Create augmented matrix
    N(:,1:nRows)=A
    N(:,nRows+1)=B

    !Row reduction
    do j=1,nRows

       !Swap current top row with row with largest leftmost entry to reduce error and avoid zeros.
       maxRow=j+maxloc(N(j:,j),1)-1
       if (maxRow/=j) then
          call swap(N,j,maxRow)
       end if

       !Check for 0 eigenvalues
       if (N(j,j)==0.0) then
          write (*,*) "Error, singular matrix."
          return
       end if
 
       !Reduce all rows below current top row
       do i=j+1,nRows
          if (N(i,j)/=0.0) then
             fact=-N(i,j)/N(j,j)
             call add(N,i,j,fact)
             N(i,j)=0.0
          end if
       end do
    end do

    !forall (j=1:nRows) x(j)=N(j,nRows+1)/N(j,j)

    !Back substitution
    x=0.0
    do i=1,nRows
       do j=nRows+2-i,nRows+1

          ! Adds entries left to right. When last entry is reached, current sum value is negated before final entry is added.
          if (j==nRows+1) then
             x(nRows+1-i)=-x(nRows+1-i)+N(nRows+1-i,j)/N(nRows+1-i,nRows+1-i)
          else
             x(nRows+1-i)=x(nRows+1-i)+N(nRows+1-i,j)*x(j)/N(nRows+1-i,nRows+1-i)
          end if
       end do
    end do
  end subroutine gjsolve

end
