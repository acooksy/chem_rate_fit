	subroutine gaussj(a,n,np,b,m,mp,ier)
c	matrix inversion
	integer*4 m,mp,n,np,i,icol,irow,j,k,l,ll
	integer*4 indxc(np),indxr(np),ipiv(np),ier
	real*8 a(np,np),b(mp),big,dum,pivinv
c	parameter definitions
c	 a = input matrix AND inverted output matrix
c	 np = array dimension of matrix
c	 n = occupied dimension of matrix 
c	 b = solutions array??
c	 mp = array dimension of b
c	 m = occupied dimension of b
	ier=0
	do j=1,n
	 ipiv(j)=0
	enddo
	do i=1,n
	 big=0.d0
	 do j=1,n
	  if(ipiv(j).ne.1)then
	   do k=1,n
	    if (ipiv(k).eq.0) then
	     if (abs(a(j,k)).ge.big)then
	      big=abs(a(j,k))
	      irow=j
	      icol=k
	     endif
	    else if (ipiv(k).gt.1) then
	     ier = 129
	    endif
	   enddo
	  endif
	 enddo
	 ipiv(icol)=ipiv(icol)+1
	 if (irow.ne.icol) then
	  do l=1,n
	   dum=a(irow,l)
	   a(irow,l)=a(icol,l)
	   a(icol,l)=dum
	  enddo
	  dum=b(irow)
	  b(irow)=b(icol)
	  b(icol)=dum
	 endif
	 indxr(i)=irow
	 indxc(i)=icol
	 if (a(icol,icol).eq.0.d0) ier=129
	 pivinv=1.d0/a(icol,icol)
	 a(icol,icol)=1.d0
	 do l=1,n
	  a(icol,l)=a(icol,l)*pivinv
	 enddo
	 b(icol)=b(icol)*pivinv
	 do ll=1,n
	  if(ll.ne.icol)then
	   dum=a(ll,icol)
	   a(ll,icol)=0.d0
	   do l=1,n
	    a(ll,l)=a(ll,l)-a(icol,l)*dum
	   enddo
	   b(ll)=b(ll)-b(icol)*dum
	  endif
	 enddo
	enddo
	do l=n,1,-1
	 if(indxr(l).ne.indxc(l))then
	  do k=1,n
	   dum=a(k,indxr(l))
	   a(k,indxr(l))=a(k,indxc(l))
	   a(k,indxc(l))=dum
	  enddo
	 endif
	enddo
	return
	end
C  (C) Copr. 1986-92 Numerical Recipes Software =$j7,#11+39.d0

