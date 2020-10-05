	subroutine uncert(stddev,imx,imy,xjac,x,iy,ipar,unc,ist)
c	My major uncertainties routine, adapted from a spectrum analysis
c	 program by Jeff Owrutsky in Rich Saykally's lab ca. 1986.
c	 See older versions for IMSL library calls; this version is written
c	 only for use with Preuss' gaussj.f.    
c	 Andrew Cooksy; San Diego State University 2005
c	stddev = standard deviation, computed in calling program
c	imx  = parameter array dimension
c	imy  = observables array dimension
c	x = parameter array
c	xjac = jacobian matrix
c	iy = # observables
c	ipar = # fitted parameters 
c       unc = experimental uncertainties array

	implicit double precision (a-h,j-z)
	real*8 xjac(800,108),xjact(108,800),var(108,108),mat(108,108),
     2	  matinv(108,108),wkarea(20000),x(108),stddev,unc(800),
     3    s95(108),s99(108),s9995(108),dmat(108)
	integer*4 i1,i2,iy,ipar,iunc,iir,ibad,iidbd1(108),iidbd2(108),
     2    idf(108)
	character*3 npar(108)
	common/ratea/npar

c       evaluate degrees of freedom for uncertainties evaluation by checking
c       jacobian matrix for number of transitions correlated to each parameter
        do i1=1,ipar
          idf(i1)=-1
          do i2=1,iy
            if(dabs(xjac(i2,i1)*x(i1)).ge.unc(i2))idf(i1)=idf(i1)+1
          enddo
	enddo

c	Make the transpose Jacobean matrix xjact
	do i1=1,iy
	  do i2=1,ipar
	    xjact(i2,i1)=xjac(i1,i2)
	  enddo
	enddo

c	Make the inverse variance-covariance matrix matinv
	do i1=1,ipar
         dmat(i1) = 0d0
	  do i2=1,ipar
	    matinv(i1,i2) = 0d0
	    do i3=1,iy
  	      matinv(i1,i2)=xjact(i1,i3)*xjac(i3,i2)+matinv(i1,i2)
	    enddo
	  enddo
	enddo
	iidgt=6

c	Invert matinv to get the variance-covariance matrix
        call gaussj(matinv,ipar,imx,dmat,ipar,imx,ier)
        do i1=1,ipar
         do i2=1,ipar
          mat(i1,i2) = matinv(i1,i2)
         enddo
        enddo

	do i1=1,ipar
	  do i2=1,ipar
	  var(i1,i2)=mat(i1,i2)*stddev**2d0
	  enddo
	enddo

c	computing the parameter uncertainties
	write(7,900)
	write(7,905)iidgt
	if(iir.eq.0)write(7,*)' Successful matrix inversion.'
	if(iir.eq.129)write(7,*)' Warning: singular covariance matrix.'
	if(iir.eq.131)write(7,*)
     2    ' Warning: matrix does not iteratively invert.'
	if(iir.eq.34)write(7,*)
     2    ' Warning: accuracy test failed in inversion.'
	if(iir.ne.0)write(6,*)' Error in covariance matrix inversion.'

	write(7,*)' '
	write(7,907)
	
c	check for uncoupled parameters
	do i1=1,ipar
	  if(var(i1,i1).eq.0d0) write(7,982)npar(i1)
	enddo

c	write uncertainties to output
	ibad=0
	do i1=1,ipar
c     *No* idea why the following line was written the way it was...
c	  if(var(i1,i1).ne.0d0)dvar=5.00*log10(3d0*var(i1,i1))
	  if(var(i1,i1).ne.0d0)dvar=log10(dsqrt(var(i1,i1)))
	  iunc=int(dvar)
	  if((iunc.gt.7).or.(iunc.lt.-8))
     2     write(7,950)npar(i1),x(i1),dsqrt(var(i1,i1)),
     3     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if((iunc.le.7).and.(iunc.gt.1))
     2     write(7,951)npar(i1),x(i1),dsqrt(var(i1,i1)),
     3     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(((iunc.eq.0).or.(iunc.eq.1)).and.(dvar.ge.0d0))
     2     write(7,952)npar(i1),x(i1),
     3     dsqrt(var(i1,i1)),2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(((iunc.eq.0).or.(iunc.eq.1)).and.(dvar.lt.0d0))
     2     write(7,953)npar(i1),x(i1),
     3     dsqrt(var(i1,i1)),2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-1)write(7,954)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-2)write(7,955)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-3)write(7,956)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-4)write(7,957)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-5)write(7,958)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-6)write(7,959)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-7)write(7,960)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
	  if(iunc.eq.-8)write(7,961)npar(i1),x(i1),dsqrt(var(i1,i1)),
     2     2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
c	  flag parameters with uncertainties of at least one third
	  if(dsqrt(var(i1,i1))/dabs(x(i1)).gt.0.33d0)then
	    ibad=ibad+1
	    iidbd1(ibad)=i1
	  endif
	enddo

	if(ist.eq.1)then
	 call studentt(var,idf,i,ipar,s95,s99,s9995)
	 write(7,910)
	 write(7,912)
	 do i1=1,ipar
	  if(s95(i1).gt.0d0)dstt=log10(s95(i1))
	  istt=int(dstt)
	  if(s95(i1).le.0d0)istt=10
	  if((istt.gt.7).or.(istt.lt.-8))
     2     write(7,950)npar(i1),x(i1),s95(i1),s99(i1),s9995(i1)
	  if((istt.le.7).and.(istt.ge.1))
     2     write(7,951)npar(i1),x(i1),s95(i1),s99(i1),s9995(i1)
	  if((istt.eq.0).and.(dstt.ge.0d0))write(7,952)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if((istt.eq.0).and.(dstt.lt.0d0))write(7,953)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-1)write(7,954)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-2)write(7,955)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-3)write(7,956)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-4)write(7,957)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-5)write(7,958)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-6)write(7,959)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-7)write(7,960)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	  if(istt.eq.-8)write(7,961)npar(i1),x(i1),
     2     s95(i1),s99(i1),s9995(i1)
	 enddo

	endif

	if(ibad.eq.1)then
	 write(7,970)
	 write(7,975)(npar(iidbd1(1)))
	endif
	if(ibad.gt.1)then
	 write(7,971)
	 write(7,975)(npar(iidbd1(i1)),i1=1,ibad)
	endif

c	computing covariances
	if(ipar.eq.1) goto 800
	ibad=0
	do i1=1,ipar
	  do i2=1,i1-1
	    if((var(i1,i1).ne.0d0).and.(var(i2,i2).ne.0d0))
     2        var(i1,i2)=var(i1,i2)/dsqrt(var(i1,i1)*var(i2,i2))
	    if((var(i1,i1).eq.0d0).or.(var(i2,i2).eq.0d0))
     2        var(i1,i2)=0d0
	    if(dabs(var(i1,i2)).gt.0.96d0)then
	      ibad=ibad+1
	      iidbd1(ibad)=i1
	      iidbd2(ibad)=i2
	    endif
	  enddo
	enddo
	write(7,980)

	icr1=1
	icr2=8
400	if(icr2.gt.ipar-1)icr2=ipar-1
	write(7,985) (npar(i2),i2=icr1,icr2)
	do i1=icr1+1,ipar
	  i3=i1-1
	  if(i3.gt.icr1+7)i3=icr1+7
	  write(7,987) npar(i1),(var(i1,i2),i2=icr1,i3)
	enddo
	if(icr2.ne.ipar-1)then
	  icr1=icr1+8
	  icr2=icr2+8
	  goto 400
	endif

	if(ibad.gt.0)then
	  write(7,990)
	  do i1=1,ibad
	  write(7,995)npar(iidbd2(i1)),npar(iidbd1(i1)),
     2      var(iidbd1(i1),iidbd2(i1)) 
	  enddo
	endif
800	return

900	format(/,'  UNCERTAINTIES',/)
905	format(' Covariance matrix elements good to roughly ',i2,
     2    ' digits.')
907	format(17x,'value',12x,'one sigma',7x,'two sigma',5x,
     2    'three sigma')
910	format(/,' Student t-distribution confidence levels:',/)
912	format(17x,'value',15x,'95%',13x,'99%',13x,'99.95%')
950	format(1x,a3,11x,e13.7,4x,3(e8.2,8x))
951	format(1x,a3,2x,f11.0,11x,3(f7.0,9x))
952	format(1x,a3,2x,f12.1,13x,3(f5.1,11x))
953	format(1x,a3,2x,f13.2,12x,3(f6.2,10x))
954	format(1x,a3,2x,f14.3,11x,3(f7.3,9x))
955	format(1x,a3,2x,f15.4,10x,3(f8.4,8x))
956	format(1x,a3,2x,f16.5,9x,3(f9.5,7x))
957	format(1x,a3,2x,f17.6,8x,3(f10.6,6x))
958	format(1x,a3,2x,f18.7,7x,3(f11.7,5x))
959	format(1x,a3,2x,f19.8,6x,3(f12.8,4x))
960	format(1x,a3,2x,f20.9,5x,3(f13.9,3x))
961	format(1x,a3,2x,f21.10,4x,3(f14.10,2x))
970	format(/,
     2    ' The following parameter is essentially undetermined:')
971	format(/,
     2    ' The following parameters are essentially undetermined:')
975	format(20(1x,a3),/)
980	format(/,' CORRELATION MATRIX')
982	format(1x,a3,' is completely decoupled from data.',
     2    '  No uncertainty evaluated.')
985 	format(/,6x,8(6x,a3))
987 	format(a3,8(3x,f6.3))
990	format(/,' The following parameters are highly correlated:')
995	format(2x,a3,'  and  ',a3,4x,f5.2)

	end

