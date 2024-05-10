	subroutine xmrqmin(y,omc,unc,x,dyda,ssq,eps,tol,
     2                     ifit,ix,iy,imx,imy,imax,ier)
c     driver for routine mrqmin
c	[example program EXTENSIVELY modified to manage mrqmin runs]
c  WARNING: dynamic allocation of array dimensions is NOT supported by
c       the original F77, so beware -- some compilers will require
c       that these arrays be explcitly dimensioned, and may not warn
c       you during compilation.  
	implicit double precision(a-h,j-z)
c	integer*4 iy,ix
	real*8 alamda,ssq,fgauss,gasdev,ossq,y(imy),omc(imy),
     2     unc(imy),x(imx),dyda(imy,imx),tol
	integer*4 i,ifit(imx),idum,iter,itst,j,k,mfit,iy,ier
c	 parameter definitions:
c	  y = data array
c	  omc = residuals array
c	  unc = data uncertainty array
c	  iy = # data points
c	  x = *full* parameter array (not just fitted parameters)
c	  ifit = fitting flags array 
c	  ix = # parameters
c	  dyda = jacobian matrix
c	  imx = maximum # parameters (= dimension of covar)
c	  ssq = chi squared
c	  alamda = scaling factor for changes to parameters
c	  res = the residuals subroutine, with format:
c           res(iy,ix,imy,imx,y,x,omc,ssq)
c	      iy = # observations
c	      ix = # parameters
c	      imy = observation array dimension
c	      imx = parameter array dimension
c	      y = data array
c	      x = parameter array
c	      omc = residuals array
c	      ssq = chi squared

	alamda=-1
	call mrqmin(y,omc,unc,iy,x,ifit,ix,dyda,imy,imx,ssq,
     2    alamda,ier)
	k=1
	itst=0
1	k=k+1
	ossq=ssq
	call mrqmin(y,omc,unc,iy,x,ifit,ix,dyda,imy,imx,ssq,
     2    alamda,ier)
c	Convergence complete when ossq-ssq is positive and any convergence 
c	 criteria are satisfied for four consecutive interations.
	if (ssq.gt.ossq) then
	 itst=0
	else if (ssq.le.ossq) then
	 if ((dabs(ossq-ssq).lt.eps).or.(ssq.lt.tol)) then
	  itst=itst+1
	 endif
	endif
	if ((itst.lt.4).and.(k.le.imax)) then
	 goto 1
	endif
	alamda=0.0
	call mrqmin(y,omc,unc,iy,x,ifit,ix,dyda,imy,imx,ssq,
     2    alamda,ier)
	if (k.gt.imax) ier=133
	return
	end
C  (C) Copr. 1986-92 Numerical Recipes Software =$j7,#11+39.


	subroutine mrqmin(y,omc,unc,iy,x,ifit,ix,dyda,imy,imx,ssq,
     2    alamda,ier)
	implicit double precision(a-h,j-z)
c	single step of non-linear least-squares fitting
c	Note: it's up to us to check the fitting criteria, count the
c	 number of iterations, etc.
	integer*4 ix,imx,iy
	real*8 ssq,x(imx),alpha(imx,imx),covar(imx,imx),
     2    unc(imy),y(imy),omc(imy),ossq,
     3    alamda,dyda(imy,imx)
	real*8 atry(imx),beta(imx),da(imx) 
	integer*4 ifit(imx),j,k,l,m,mfit,ier
c	save ossq,atry,beta,da,mfit
c	parameter definitions:
c	 y = data array
c	 omc = residuals array
c	 unc = data uncertainty array
c	 iy = # data points
c	 x = *full* parameter array (not just fitted parameters)
c	 ifit = fitting flags 
c	 ix = # parameters
c	 covar = covariance matrix
c	 alpha = jacobian matrix
c	 imx = maximum # parameters (= dimension of covar)
c	 ssq = chi squared
c	 alamda = scaling factor for changes to parameters

	if(alamda.lt.0.d0)then
	 mfit=0
	 do j=1,ix
	  if (ifit(j).ne.0) mfit=mfit+1
	 enddo
	 alamda=0.001d0
	 call mrqcof(y,omc,unc,iy,x,ifit,ix,alpha,beta,imx,imy,
     2     ssq,dyda,alamda)
	 ossq=ssq
	 do j=1,ix
	  atry(j)=x(j)
	 enddo
	endif
	j=0
	do l=1,ix
	 if(ifit(l).ne.0) then
	  j=j+1
	  k=0
	  do m=1,ix
	   if(ifit(m).ne.0) then
	    k=k+1
	    covar(j,k)=alpha(j,k)
	   endif
	  enddo
	  covar(j,j)=alpha(j,j)*(1.d0+alamda)
	  da(j)=beta(j)
	 endif
	enddo
	call gaussj(covar,mfit,imx,da,1,1,ier)
	if(alamda.eq.0.d0)then
	 call covsrt(covar,imx,ix,ifit,mfit)
	 return
	endif
	j=0
	do l=1,ix
	 if(ifit(l).ne.0) then
	  j=j+1
c       Cooksy kluge programming to prevent bad overshoots...
c	  This opriginally assumed all x's were > 0, but that
c	  doesn't necessarily hold for t0
          if (x(l).ne.0d0) then 
           do while (dabs(da(j)).gt.dabs(x(l))/1.5d0) 
            da(j) = da(j)/1.5d0
           enddo
	  endif
	  atry(l)=x(l)+da(j)
	 endif
	enddo
	call mrqcof(y,omc,unc,iy,atry,ifit,ix,covar,da,imx,imy,
     2    ssq,dyda,alamda)
	if(ssq.lt.ossq)then
	 alamda=0.4d0*alamda
	 ossq=ssq
	 j=0
	 do l=1,ix
	  if(ifit(l).ne.0) then
	   j=j+1
	   k=0
	   do m=1,ix
	    if(ifit(m).ne.0) then
	     k=k+1
	     alpha(j,k)=covar(j,k)
	    endif
	   enddo
	   beta(j)=da(j)
	   x(l)=atry(l)
	  endif
	 enddo
	else
	 alamda=03.0d0*alamda
	 ssq=ossq
	endif
	return
	end
C  (C) Copr. 1986-92 Numerical Recipes Software =$j7,#11+39.d0



	subroutine mrqcof(y,omc,unc,iy,x,ifit,ix,alpha,beta,imx,
     2    imy,ssq,dyda,alamda)
	implicit double precision(a-h,j-z)
c	called by mrqmin
c	[Modified to make "funcs" call compatible with existing subroutine res]
	integer*4 ix,iy,imx
	real*8 ssq,x(imx),alpha(imx,imx),beta(imx),
     2    x2(imx),dy,unc2i,wt,omc(imy),omc2(imy),dyda(imy,imx),
     3    unc(imy),y(imy),alamda,xchk(108)
	integer*4 nalp,mfit,i,j,k,l,m,ifit(imx)
c	integer*4 iy
	mfit=0
	do j=1,ix
	 if (ifit(j).ne.0) mfit=mfit+1
	enddo
	do j=1,mfit
	 do k=1,j
	  alpha(j,k)=0.d0
	 enddo
	 beta(j)=0.d0
	enddo
        do i=1,ix
         xchk(i)=x(i)
        enddo
        call res(iy,ix,imy,imx,y,x,omc,ssq)

c  [this section added to evaluate derivatives numerically]
	do l=1,ix
	 x2(l) = x(l)
	enddo
	do l=1,ix
         if (ifit(l).ne.0) then
	  if (x(l).eq.0d0) delx = 1.0d-4
          if (x(l).ne.0d0) delx = x(l)/1.0d4
          x2(l) = x(l)+delx
          do i=1,ix
           xchk(i)=x2(i)
          enddo
	  call res(iy,ix,imy,imx,y,x2,omc2,ssq)
	  do i=1,iy
	   dyda(i,l) = (omc(i)-omc2(i))/delx
	  enddo
	  x2(l) = x(l)
         else if (ifit(l).eq.0) then
	  do i=1,iy
	   dyda(i,l) = 0d0
	  enddo
         endif
	enddo

	do i=1,iy
	 unc2i=1.d0/(unc(i)*unc(i))
	 domc=omc(i)
	 j=0
	 do l=1,ix
	  if(ifit(l).ne.0) then
	   j=j+1
	   wt=dyda(i,l)*unc2i
	   k=0
	   do m=1,l
	    if(ifit(m).ne.0) then
	     k=k+1
	     alpha(j,k)=alpha(j,k)+wt*dyda(i,m)
	    endif
	   enddo   
	   beta(j)=beta(j)+domc*wt
	  endif
	 enddo
	enddo
	do j=2,mfit
	 do k=1,j-1
	  alpha(k,j)=alpha(j,k)
	 enddo
	enddo
	return
	end
C  (C) Copr. 1986-92 Numerical Recipes Software =$j7,#11+39.d0



	subroutine covsrt(covar,imx,ix,ifit,mfit)
	implicit double precision(a-h,j-z)
c	called by mrqmin
	integer*4 ix,imx
	real*8 covar(imx,imx),swap
	integer*4 mfit,ifit(imx),i,j,k

	do i=mfit+1,ix
	 do j=1,i
	  covar(i,j)=0.d0
	  covar(j,i)=0.d0
	 enddo
	enddo
	k=mfit
	do j=ix,1,-1
	 if(ifit(j).ne.0)then
	  do i=1,ix
	   swap=covar(i,k)
	   covar(i,k)=covar(i,j)
	   covar(i,j)=swap
	  enddo
	  do i=1,ix
	   swap=covar(k,i)
	   covar(k,i)=covar(j,i)
	   covar(j,i)=swap
	  enddo
	  k=k-1
	 endif
	enddo
	return
	end
C  (C) Copr. 1986-92 Numerical Recipes Software =$j7,#11+39.d0



