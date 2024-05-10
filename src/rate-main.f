c       RATE.F VERSION 2.03  Jul. 30 2018
c       copyright 2004-2024 by Andrew L. Cooksy

c       I.      VERSION NUMBERS
c       2.03 Jul 30 2018
c	     Clean-up the t0 handling; fitting t0 is better now, for
c	     thing negative values are now allowed.
c	     Still some issues: t0 can't exceed t of the first data
c	     point or we get a datum excluded error.
c       2.02 Oct 10 2005
c            Fixed bug (I hope) in integrate subroutine of rate-res.f
c            that overshot contribution to mechanism of A + A -> product
c            reactions vs. A + B -> product.  I now believe that the
c            rate law for these should be d[A}/dt = - k[A]^2 when I
c            had been using d[A}/dt = - 2 k[A]^2.  
c            I just had assumed two A molecules consumed per reaction.  
c            But each A+A collision is per two A molecules
c            (while the A+B collision is per one A molecule), and so
c            the overall consumption from A+A should be divided by two 
c            molecules in each collision *and* multiplied by two molecules 
c            consumed, leaving the factor out front equal to one.
c       2.01 Jul 26 2005
c            Fixed bug in rate-res.f to allow fitting of data using
c            ilog=1; using linear time steps with a converged dt is 
c            generally too demanding for fitting data with a fast
c            induction time.
c       2.00 Adaptation to allow clumsy least squares fitting of rate
c            constants using Preuss' least squares subroutines.
c       1.00 Originally put together for Kinetics and Physical Chemistry
c            courses to predict [X] vs. t plots for known mechanisms,
c            adapted to include both Euler and 4th order Runge-Kutta
c            integration methods, linear and logarithmic time steps,
c            including Arrhenius and pre-exponential T-dependence in 
c            rate constants.

c       II.     BRIEF DESCRIPTION
c	RATE.F predicts chemical concentrations in a general
c	network of elementary reactions, provided all relevant
c	rate laws are known.  The program uses brute-force numerical 
c	integration of the rate laws, and can be tested for accuracy
c	only by decreasing the time increment delta t.
c       Output files:
c        rate.log is the program output summary, containing the 
c               fitted constant values, convergence errors and/or criteria,
c               other program messages.
c        rate.out is the predicted [X] vs. t output, based on the
c               fitted constants annd according to the input parameters
c               ipt, iout, ilog. 
c        rate#.out is the predicted [X] vs. t output for only the #th chemical
c               component selected by ipx, if iout=2.

c       III.    PARTIAL VARIABLE DICTIONARY
c	Variables:
c	xa(i) = previous set of concentrations
c	xt(i)  = current set of concentrations
c	delx(i) = set of incremental corrections to concentrations
c	xp(i)  = set of concentrations to be printed (if iout = 1)
c       irx = # reactions in mechanism (# fwd and rev pairs)
c       icnd = # initial conditions (varying temperature or initial concs)
c	temp = temperature T (K)
c	forward reaction rate(i) = kf(i) * T^(kftp(i)) * exp(-kftx/T)
c	reverse reaction rate(i) = kr(i) * T^(krtp(i)) * exp(-krtx/T)
c	icrf(i) = order of the forward half-reaction(i)
c	icrr(i) = order of the reverse half-reaction(i)
c	icrfi(i,j) = set of reactants for forward half-reaction(i)
c	icrri(i,j) = set of reactants for reverse half-reaction(i)
c	npar(i) = set of parameter names
c	t = time
c	t0 = initial time for reaction (normally zero, but sometimes this is unknown);
c	     this is the time at which the initial concentrations are valid;
c	     all integrations begin from t0, even if it's a fitted parameter
c	tmin = time increment for ilog=1; must be 0 < tmin < tmax
c	tmax = final time for calc'd values (relative to t=0, not necessarily to t0);
c            all observed values must fall between t0 (or tmin for ilog=1) and tmax 
c	     to be counted in fit
c	dt = time increment for calc'd values, = (tmax-t0)/it for ilog=0
c	it = total # samples
c	ipt = # increments to sample between printing concentrations
c	iout = flag selecting output format: 0 more precise; 1 more readable;
c		2 up to six distinct output files, 3 one single CSV file
c	ilog = flag selecting time increments: 0 linear; 1 logarithmic
c	icase = flag setting method used: 0 Euler; 1 4th order Runge-Kutta
c       ifit = array of fitting flags (1 to fit, 0 to fix at input values)
c       itfit = total number of fitted parameters, for estimating std dev
c
c       Remaining development (assuming all this is working:
c         Split to allow either predictive and fitting options
c         One failure mode that surprised me: listing several
c           reactions but only trying to fit the first one can
c           lead to a register overflow problem unless all the
c           later rate constants are set to zero.  The later
c           concentrations can diverge exponentially as the
c           program adjusts the one rate constant independnently
c           of later ones, and the program blindly calculates *all*
c           concentrations.  Some sort of trap for this would be nice.
c         Fitting mechanism to multiple trials with different initial
c           concentrations and temperatures.
c         Fitting t0 is hard under the current algorithm, because the 
c           change in t0 used to measure dyda in xmrqmin is too small
c           to be resolved when the time is digitized.  There should
c           be a better way of treating this term so it can be fit.
c         Lots of traps possible in input.  The observed concentrations
c           need to fall between t0 and tmax, for example, or you
c           waste potentially lots of time.

c*********************** MAIN PROGRAM *************************************

	implicit double precision (a-h,j-z)
c	external res
	real*8 kf(40),kr(40),ktfx(40),ktrx(40),ktfp(40),x0(108),
     2    ktrp(40),xp(28),delx(40),delx2(40),delx3(40),delx4(40),
     3    x(108),y(2000),yt(2000),omc(2000),unc(2000),xjac(2000,108)
	integer*4 icrf(40),icrr(40),icrfi(40,5),icrri(40,5),ipx(28),
     2    ifit(108),ya(2000),yc(2000)
	character*3 npar(108)
	character*80 garbage,header

	common/rates/t,dt,ktfx,ktrx,ktfp,ktrp,tmin,tmax,
     2    iconc,irx,icnd,icrf,icrr,icrfi,icrri,it,ipt,iout,ilog,
     3    icase,ipx,itfit
        common/rate/yc,ya,yt,y
        common/ratea/npar
        common/max/imrx,imconc,imy,imx

c        data imrx,imconc,imy,imx/10,10,100,10/
       data imrx,imconc,imy,imx/40,27,2000,108/

c	initialization
	do i=1,imx
	  x(i) = 0d0
	  ifit(i)=0
	enddo
        do i=1,imy
          y(i) = 0d0
          yt(i) = 0d0
          ya(i) = 0
          yc(i) = 0
          unc(i) = 0d0
          omc(i) = 0d0
        enddo

	open(unit=4,file='rate.dat')
	read(4,900)header
	read(4,900)garbage
	read(4,*)iconc,irx,icnd
	read(4,900)garbage
	read(4,*)it,tmin,tmax,iy,ipt
	read(4,900)garbage
	read(4,*)iout,ilog,icase,ipred
        read(4,900)garbage
        read(4,*)eps,tol,imax

c       read the initial concentrations
	read(4,*)(npar(i),i=1,iconc+2)
        do i1=1,icnd
	 read(4,*)(x(i),i=(i1-1)*(iconc+2)+1,i1*(iconc+2))
	 read(4,*)(ifit(i),i=(i1-1)*(iconc+2)+1,i1*(iconc+2))
         if (i1.gt.1) then
          do i2=1,iconc+2
           npar((i1-1)*(iconc+2)+i2)=npar(i2)
          enddo
         endif
        enddo
	read(4,*)(ipx(i),i=1,iconc+1)
        iflg = 0
        do i1=1,icnd 
         if (ifit((i1-1)*(iconc+2)+2).eq.1) iflg = 1
        enddo
        if (iflg.eq.1) then
         write(6,*)"Warning: fitting t0 usually won't work because time"
         write(6,*)"is a discrete (rather than continuous variable)."
         write(6,*)"Discretization of the curves is not fine enough to"
         write(6,*)"detect changes w.r.t. t0 except for dt/t0 < 10E-7."
         iflg = 0
        endif

c       ixcm is just the x array index of the last concentration
        ixcm = icnd*(iconc+2)

c	read in reaction mechanism
	read(4,900)garbage
	do i=1,irx
	  read(4,*)npar(ixcm+2*i-1),kf(i),ifit(ixcm+2*i-1),ktfp(i),
     2             ktfx(i),icrf(i),(icrfi(i,i1),i1=1,icrf(i))
	  read(4,*)npar(ixcm+2*i),kr(i),ifit(ixcm+2*i),ktrp(i),ktrx(i),
     2             icrr(i),(icrri(i,i1),i1=1,icrr(i))
	  x(ixcm+2*i-1) = kf(i)
	  x(ixcm+2*i) = kr(i)
	enddo
	do i=1,imx
	  x0(i) = x(i)
	enddo

c	read in experimental concentration data: species, time, value, unc
	read(4,900)garbage
	do i=1,iy
	  read(4,*)yc(i),ya(i),yt(i),y(i),unc(i)
	enddo

    	close(unit=4)
c       end reading from input file

c       ix is the total number of parameters for this run
c	ix = 2*irx + ixcm    = fwd & rev rate constants per reaction + icnd x
c	                           (temp + t0 + initial concentrations) 
	ix = 2*irx + ixcm     

        if (ipred.eq.1) then
         call res(iy,ix,imy,imx,y,x,omc,ssq)
         stop
        endif

c       find number of fitted parameters
        itfit = 0
        do i=1,ix
          itfit = itfit + ifit(i)
        enddo

        write(6,*)' Least squares fitting.'
                write(6,*)' '

        ier = 0
        call xmrqmin(y,omc,unc,x,xjac,ssq,eps,tol,
     2               ifit,ix,iy,imx,imy,imax,ier)

        open(unit=7,file='rate.log')
        write(7,900)header
        stddev = 0d0
        if (iy.gt.itfit+1) stddev = dsqrt(ssq/dfloat(iy-itfit-1))
        write(7,920)ssq,ier,stddev
c Error message statements
        if((ier.eq.129).or.(ier.eq.132))write(7,*)
     2   ' Singularity in Jacobean.'
        if((ier.eq.133).or.(ier.eq.4))write(7,*)' imax exceeded.'
        if(ier.eq.39)write(7,*)' Convergence criteria too stringent.'
        if(ier.eq.38)write(7,*)' Zero Jacobean.'
        if(ier.eq.7)write(7,*)' Convergence too slow.'
        if(ier.eq.6)write(7,*)' Five steps taken with max step length.'
        if(ier.eq.2)
     2   write(7,*)' Convergence to non-critical point suspected.'
        if(ier.eq.0)write(7,*)' Successful convergence.'
        if(ier.ne.0)write(7,919)
        write(7,949)
        do i=1,irx
         write(7,*)" "
         if ((icrf(i).eq.1).and.(icrr(i).eq.1)) then
          write(7,951)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.1).and.(icrr(i).eq.2)) then
          write(7,952)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.1).and.(icrr(i).eq.3)) then
          write(7,953)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.2).and.(icrr(i).eq.1)) then
          write(7,954)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.2).and.(icrr(i).eq.2)) then
          write(7,955)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.2).and.(icrr(i).eq.3)) then
          write(7,956)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.3).and.(icrr(i).eq.1)) then
          write(7,957)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.3).and.(icrr(i).eq.2)) then
          write(7,958)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.3).and.(icrr(i).eq.3)) then
          write(7,959)i,npar(ixcm+2*i-1),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).gt.3).or.(icrr(i).gt.3)) then
          write(7,950)i,npar(ixcm+2*i-1),npar(ixcm+2*i),
     2     (npar(icrfi(i,i1)+2),i1=1,icrf(i)), 
     3     (npar(icrri(i,i1)+2),i1=1,icrr(i))
         endif
        enddo
        write(7,*)' '
        write(7,*)' Fit criteria'
        write(7,*)' '
        write(7,921)iy,ix
        dt = (tmax-t0)/dfloat(it)
        write(7,922)it,dt
        write(7,923)eps,tol,imax
        write(7,960)
c       have to call res again to refresh the omc array, which otherwise
c        may not correspond to the final ssq value
        call res(iy,ix,imy,imx,y,x,omc,ssq)
        do i=1,iy
         ycalc = y(i)-omc(i)
         write(7,961)yc(i),npar(ya(i)+2),yt(i),y(i),ycalc,omc(i),unc(i)
        enddo
        write(7,962)
        do i=1,ix
         if (i.le.ixcm) write(7,963)npar(mod(i-1,iconc+2)+1),x0(i),x(i)
         if (i.gt.ixcm) write(7,963)npar(i),x0(i),x(i)
        enddo

        ipar=0
        do i1=1,ix
         if (ifit(i1).eq.1) then
          ipar=ipar+1
          x(ipar)=x(i1)
          npar(ipar)=npar(i1)
          do i2=1,iy
           xjac(i2,ipar)=xjac(i2,i1)
          enddo
         endif
        enddo
        ist = 0
        call uncert(stddev,imx,imy,xjac,x,iy,ipar,unc,ist)

	close(unit=7)
900     format(a80)
919     format(' WARNING: Failure to converge with Preuss least ',/,
     2   'squares routine does NOT necessarily report best values.')
920     format('  ssq = ',g9.3,5x,'  ier = ',i4,/,'  std dev = ',g10.4)
921     format('  # observables =',i4,', # variables =',i3)
922     format('  # time steps =',i9,', delta t =',g12.6,' (linear t)')
923     format('  Convergence criteria: eps =',g8.2,3x,'tol =',
     2   g8.2,3x,'imax =',i5,'.')
949     format(/,' reaction mechanism')
950     format('  reaction #',i2,'  rate constants ',a3,', ',a3,/,
     2   ' reactants: ',5(a3,1x),/,'  products: ',5(a3,1x))
951     format('  reaction #',i2,6x,a3,/,3x,
     2                 12x,a3,' <---> ',a3,/,20x,a3)
952     format('  reaction #',i2,6x,a3,/,3x,
     2                 12x,a3,' <---> ',a3,' + ',a3,/,20x,a3)
953     format('  reaction #',i2,6x,a3,/,3x,
     2                 12x,a3,' <---> ',a3,' + ',a3,' + ',a3,/,20x,a3)
954     format('  reaction #',i2,6x,a3,/,3x,
     2         6x,a3,' + ',a3,' <---> ',a3,/,20x,a3)
955     format('  reaction #',i2,6x,a3,/,3x,
     2         6x,a3,' + ',a3,' <---> ',a3,' + ',a3,/,20x,a3)
956     format('  reaction #',i2,6x,a3,/,3x,
     2         6x,a3,' + ',a3,' <---> ',a3,' + ',a3,' + ',a3,/,20x,a3)
957     format('  reaction #',i2,6x,a3,/,3x,
     2   a3,' + ',a3,' + ',a3,' <---> ',a3,/,20x,a3)
958     format('  reaction #',i2,6x,a3,/,3x,
     2   a3,' + ',a3,' + ',a3,' <---> ',a3,' + ',a3,/,20x,a3)
959     format('  reaction #',i2,6x,a3,/,3x,
     2   a3,' + ',a3,' + ',a3,' <---> ',a3,' + ',a3,' + ',a3,/,20x,a3)
960     format('   species',/,'conds',6x,'t',12x,'y',12x,'ycalc',7x,
     2   'omc',10x,'unc')
961     format(1x,i2,1x,a3,3(1x,g12.4),2(1x,g12.4))
962     format(/,' initial and final parameter values')
963     format(1x,a3,2(4x,g12.5))

	end

