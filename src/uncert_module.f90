       module uncert_module
!       use common_ratea
       use common_module
       use matinv_module, only : matinv

       implicit none

       contains

       subroutine uncert(stddev,imx,imy,xjac,xfit,iy,ipar,ist,ifile_log)
!       My major uncertainties routine, adapted from a spectrum analysis
!        program by Jeff Owrutsky in Rich Saykally's lab ca. 1986.
!        See older versions for IMSL library calls or use with Preuss' gaussj.f.    
!        Andrew Cooksy; San Diego State University 2005
!       stddev = weighted, unnormalized standard deviation, computed in calling program
!       imx  = parameter array dimension
!       imy  = observables array dimension
!       x = parameter array
!       xjac = jacobian matrix
!       iy = # observables
!       ipar = # fitted parameters 
!       unc = experimental uncertainties array

       integer, intent(in) :: imx,imy,iy,ipar,ist
       real(dp), intent(in) :: stddev,xfit(108)
       real(dp), intent(in) :: xjac(iy,ipar)
       integer :: i,i1,i2,i3,iunc,ibad,iidbd1(108),iidbd2(108), &
         idf(108),icr1,icr2,istt,idgt,ier,ifile_log
       real(dp) :: a(ipar,ipar),ainv(ipar,ipar),dmat(ipar),var(ipar,ipar), &
         xjact(ipar,iy)
       real(dp) :: dstt,dvar,s95(108),s99(108),s9995(108)

!       For the benefit of the Student t uncertainties only:
!       evaluate degrees of freedom for uncertainties evaluation by checking
!       jacobian matrix for number of transitions correlated to each parameter
        do i1=1,ipar
          idf(i1)=-1
          do i2=1,iy
            if(dabs(xjac(i2,i1)*xfit(i1)).ge.unc(i2))idf(i1)=idf(i1)+1
          enddo
       enddo

!       Make the transpose Jacobean matrix xjact
       do i1=1,iy
         do i2=1,ipar
           xjact(i2,i1)=xjac(i1,i2)
         enddo
       enddo

!       Make the inverse variance-covariance matrix ainv
       do i1=1,ipar
         dmat(i1) = 0d0
         do i2=1,ipar
           ainv(i1,i2) = 0d0
           do i3=1,iy
               ainv(i1,i2)=xjact(i1,i3)*xjac(i3,i2)+ainv(i1,i2)
           enddo
         enddo
       enddo

!       Invert ainv to get the variance-covariance matrix a
!        call gaussj(matinv,ipar,imx,dmat,ipar,imx,ier)
       idgt = 6
        call matinv(ipar,ainv,a,ier,idgt)
!      error handling for matinv
       if (ier /= 0) then
        write(ifile_log,*)' Error in matrix inversion:'
        if (ier == 1) then
         write(ifile_log,*)'  at least one diagonal element in covariance matrix is zero'
         return
        endif
        if (ier == 2) write(ifile_log,*)'  precision in matrix inversion is below spec idgt =',idgt
       endif
       do i1=1,ipar
         do i2=1,ipar
         var(i1,i2)=a(i1,i2)*stddev**2d0
         enddo
       enddo

!       computing the parameter uncertainties
       write(ifile_log,900)

       write(ifile_log,*)' '
       write(ifile_log,907)
       
!       check for uncoupled parameters
       do i1=1,ipar
         if(var(i1,i1).eq.0d0) write(ifile_log,982)npar(i1)
       enddo

!       write uncertainties to output
       ibad=0
       do i1=1,ipar
!     *No* idea why the following line was written the way it was...
!         if(var(i1,i1).ne.0d0)dvar=5.00*log10(3d0*var(i1,i1))
         if(var(i1,i1).ne.0d0)dvar=log10(dsqrt(var(i1,i1)))
         iunc=int(dvar)
         if((iunc.gt.7).or.(iunc.lt.-8)) &
          write(ifile_log,950)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if((iunc.le.7).and.(iunc.gt.1)) &
          write(ifile_log,951)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(((iunc.eq.0).or.(iunc.eq.1)).and.(dvar.ge.0d0)) &
          write(ifile_log,952)npar(i1),xfit(i1), &
          dsqrt(var(i1,i1)),2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(((iunc.eq.0).or.(iunc.eq.1)).and.(dvar.lt.0d0)) &
          write(ifile_log,953)npar(i1),xfit(i1), &
          dsqrt(var(i1,i1)),2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-1)write(ifile_log,954)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-2)write(ifile_log,955)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-3)write(ifile_log,956)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-4)write(ifile_log,957)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-5)write(ifile_log,958)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-6)write(ifile_log,959)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-7)write(ifile_log,960)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
         if(iunc.eq.-8)write(ifile_log,961)npar(i1),xfit(i1),dsqrt(var(i1,i1)), &
          2d0*dsqrt(var(i1,i1)),3d0*dsqrt(var(i1,i1))
!         flag parameters with uncertainties of at least one third
         if(dsqrt(var(i1,i1))/dabs(xfit(i1)).gt.0.33d0)then
           ibad=ibad+1
           iidbd1(ibad)=i1
         endif
       enddo

       if(ist.eq.1)then
        call studentt(var,idf,i,ipar,s95,s99,s9995)
        write(ifile_log,910)
        write(ifile_log,912)
        do i1=1,ipar
         if(s95(i1).gt.0d0)dstt=log10(s95(i1))
         istt=int(dstt)
         if(s95(i1).le.0d0)istt=10
         if((istt.gt.7).or.(istt.lt.-8)) &
          write(ifile_log,950)npar(i1),xfit(i1),s95(i1),s99(i1),s9995(i1)
         if((istt.le.7).and.(istt.ge.1)) &
          write(ifile_log,951)npar(i1),xfit(i1),s95(i1),s99(i1),s9995(i1)
         if((istt.eq.0).and.(dstt.ge.0d0))write(ifile_log,952)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if((istt.eq.0).and.(dstt.lt.0d0))write(ifile_log,953)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-1)write(ifile_log,954)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-2)write(ifile_log,955)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-3)write(ifile_log,956)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-4)write(ifile_log,957)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-5)write(ifile_log,958)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-6)write(ifile_log,959)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-7)write(ifile_log,960)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
         if(istt.eq.-8)write(ifile_log,961)npar(i1),xfit(i1), &
          s95(i1),s99(i1),s9995(i1)
        enddo

       endif

       if(ibad.eq.1)then
        write(ifile_log,970)
        write(ifile_log,975)(npar(iidbd1(1)))
       endif
       if(ibad.gt.1)then
        write(ifile_log,971)
        write(ifile_log,975)(npar(iidbd1(i1)),i1=1,ibad)
       endif

!       computing covariances
       if(ipar /= 1) then
       ibad=0
       do i1=1,ipar
         do i2=1,i1-1
           if((var(i1,i1).ne.0d0).and.(var(i2,i2).ne.0d0)) &
             var(i1,i2)=var(i1,i2)/dsqrt(var(i1,i1)*var(i2,i2))
           if((var(i1,i1).eq.0d0).or.(var(i2,i2).eq.0d0)) &
             var(i1,i2)=0d0
           if(dabs(var(i1,i2)).gt.0.96d0)then
             ibad=ibad+1
             iidbd1(ibad)=i1
             iidbd2(ibad)=i2
           endif
         enddo
       enddo
       write(ifile_log,980)

! Writing the covariance matrix, up to 8 elements on one line
       icr1=1
       icr2=8
       if(icr2.gt.ipar-1)icr2=ipar-1
       do while (icr2 <= ipar-1)
        write(ifile_log,985) (npar(i2),i2=icr1,icr2)
        do i1=icr1+1,ipar
         i3=i1-1
         if(i3.gt.icr1+7)i3=icr1+7
         write(ifile_log,987) npar(i1),(var(i1,i2),i2=icr1,i3)
        enddo
        if (icr2 == ipar-1) then
         exit
        else if (icr2 /= ipar-1) then
         icr1 = icr1+8
         icr2 = icr2+8
         if (icr2 > ipar-1) icr2 = ipar-1
        endif
       enddo  ! while (icr2 /= ipar-1)

! Rat out any very high correlations
       if(ibad.gt.0)then
         write(ifile_log,990)
         do i1=1,ibad
         write(ifile_log,995)npar(iidbd2(i1)),npar(iidbd1(i1)), &
           var(iidbd1(i1),iidbd2(i1)) 
         enddo
       endif
       
       endif    !  if(ipar /= 1) 
       return

900       format(/,'  UNCERTAINTIES',/)
905       format(' Covariance matrix elements good to roughly ',i2, &
         ' digits.')
907       format(17x,'value',12x,'one sigma',7x,'two sigma',5x, &
         'three sigma')
910       format(/,' Student t-distribution confidence levels:',/)
912       format(17x,'value',15x,'95%',13x,'99%',13x,'99.95%')
950       format(1x,a3,11x,e13.7,4x,3(e8.2,8x))
951       format(1x,a3,2x,f11.0,11x,3(f7.0,9x))
952       format(1x,a3,2x,f12.1,13x,3(f5.1,11x))
953       format(1x,a3,2x,f13.2,12x,3(f6.2,10x))
954       format(1x,a3,2x,f14.3,11x,3(f7.3,9x))
955       format(1x,a3,2x,f15.4,10x,3(f8.4,8x))
956       format(1x,a3,2x,f16.5,9x,3(f9.5,7x))
957       format(1x,a3,2x,f17.6,8x,3(f10.6,6x))
958       format(1x,a3,2x,f18.7,7x,3(f11.7,5x))
959       format(1x,a3,2x,f19.8,6x,3(f12.8,4x))
960       format(1x,a3,2x,f20.9,5x,3(f13.9,3x))
961       format(1x,a3,2x,f21.10,4x,3(f14.10,2x))
970       format(/, &
         ' The following parameter is essentially undetermined:')
971       format(/, &
         ' The following parameters are essentially undetermined:')
975       format(20(1x,a3),/)
980       format(/,' CORRELATION MATRIX')
982       format(1x,a3,' is completely decoupled from data.', &
         '  No uncertainty evaluated.')
985        format(/,6x,8(6x,a3))
987        format(a3,8(3x,f6.3))
990       format(/,' The following parameters are highly correlated:')
995       format(2x,a3,'  and  ',a3,4x,f5.2)

       end subroutine uncert


       subroutine studentt(var,idf,ilin,ipar,s95,s99,s9995)
!       Calculates Student's t distribution, using data table from CRC.
       real(dp) :: t600(34),t750(34),t900(34),t950(34),t975(34),t990(34)
       real(dp) :: t995(34),t999(34),t9995(34)
       real(dp), intent(in) :: var(ipar,ipar)
       real(dp), intent(out) :: s95(108),s99(108),s9995(108)
       integer :: i1,idf2(108)
       integer, intent(in) :: ilin,idf(108),ipar
       data t600/.325,.289,.277,.271,.267,.265,.263,.262,.261,2*.260, &
         2*.259,3*.258,5*.257,9*.256,.255,2*.254,.253/
       data t750/1.000,.816,.765,.741,.727,.718,.711,.706,.703,.700, &
         .697,.695,.694,.692,.691,.690,.689,2*.688,.687,2*.686, &
         2*.685,.684,.684,.684,.683,.683,.683,.681,.679,.677,.674/
       data t900/3.078,1.886,1.638,1.533,1.476,1.440,1.415,1.397,1.383, &
         1.372,1.363,1.356,1.30,1.345,1.341,1.337,1.333,1.330,1.328, &
         1.325,1.323,1.321,1.319,1.318,1.316,1.315,1.314,1.313,1.311, &
         1.310,1.303,1.296,1.289,1.282/
       data t950/6.314,2.920,2.353,2.132,2.015,1.943,1.895,1.860,1.833, &
         1.812,1.796,1.782,1.771,1.761,1.753,1.746,1.740,1.734,1.729, &
         1.725,1.721,1.717,1.714,1.711,1.708,1.706,1.703,1.701,1.699, &
         1.697,1.684,1.671,1.658,1.645/
       data t975/12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262, &
         2.228,2.201,2.179,2.160,2.145,2.131,2.120,2.110,2.101,2.093, &
         2.086,2.080,2.074,2.069,2.064,2.060,2.056,2.052,2.048,2.045, &
         2.042,2.021,2.000,1.980,1.960/
       data t990/31.821,6.965,4.541,3.747,3.365,3.143,2.998,2.896,2.821, &
         2.764,2.718,2.681,2.650,2.624,2.602,2.583,2.567,2.552,2.539, &
         2.528,2.518,2.508,2.500,2.492,2.485,2.479,2.473,2.467,2.462, &
         2.457,2.423,2.390,2.358,2.326/
       data t995/63.657,9.925,5.841,4.604,4.032,3.707,3.499,3.355,3.250, &
         3.169,3.106,3.055,3.012,2.977,2.947,2.921,2.898,2.878,2.861, &
         2.845,2.831,2.819,2.807,2.797,2.787,2.779,2.771,2.763,2.756, &
         2.750,2.704,2.660,2.617,2.576/
       data t999/318.0,22.30,10.20,7.173,5.893,5.208,4.785,4.501,4.297, &
         4.144,4.025,3.930,3.852,3.787,3.733,3.686,3.646,3.610,3.579, &
         3.552,3.527,3.505,3.485,3.467,3.450,3.435,3.421,3.408,3.396, &
         3.385,3.307,3.232,3.160,3.090/
       data t9995/636.619,31.598,12.924,8.610,6.869,5.959,5.408,5.041, &
         4.781,4.587,4.437,4.318,4.221,4.140,4.073,4.015,3.965,3.922, &
         3.883,3.850,3.819,3.792,3.767,3.745,3.725,3.707,3.690,3.674, &
         3.659,3.646,3.551,3.460,3.373,3.291/

!       initialize
       do i1=1,ipar
         s95(i1)=0d0
         s99(i1)=0d0
         s9995(i1)=0d0
       enddo

       do i1=1,ipar
        idf2(i1)=idf(i1)
        if(idf2(i1).le.0) then 
         cycle 
        endif
        if((idf(i1).ge.30).and.(idf(i1).le.39))idf2(i1)=30
        if((idf(i1).ge.40).and.(idf(i1).le.59))idf2(i1)=31
        if((idf(i1).ge.60).and.(idf(i1).le.119))idf2(i1)=32
        if((idf(i1).ge.120).and.(idf(i1).le.499))idf2(i1)=33
        if(idf(i1).ge.500)idf2(i1)=34
        s95(i1)=dsqrt(var(i1,i1))*t950(idf2(i1))/dsqrt(dfloat(idf(i1)))
        s99(i1)=dsqrt(var(i1,i1))*t990(idf2(i1))/dsqrt(dfloat(idf(i1)))
        s9995(i1)=dsqrt(var(i1,i1))*t9995(idf2(i1))/ &
           dsqrt(dfloat(idf(i1)))
       enddo
       return
       end subroutine studentt

       end module uncert_module
