!*********************** MAIN PROGRAM *************************************
       program main
       use residuals_module
       use common_module
       use uncert_module
       use minpack_module

       implicit none

       real(dp) :: eps,tol,ssq,stddev,stddvn,t0,ycalc = 0.0_dp
       real(dp) :: x0(108),xp(28),xfit(108),omc(2100),wa(82200),fvec(2100) = 0.0_dp
       integer :: i,i1,i2,iwa(40),lwa,ier,iflg,imax,info,ipar, &
         ipred,ist,ixcm,iy = 0
!       integer ::  imrx,imconc,imy,imx = 0
       integer :: ifile_dat,ifile_log,ifile_out 
       character(len=80) :: garbage,header
       real(dp),allocatable :: xjac(:,:) 

!       initialization
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

       open(newunit=ifile_dat,file='rate.dat',status="old")
       read(ifile_dat,*)header
       read(ifile_dat,*)garbage
       read(ifile_dat,*)iconc,irx,icnd
       read(ifile_dat,*)garbage
       read(ifile_dat,*)it,tmin,tmax,iy,ipt
       read(ifile_dat,*)garbage
       read(ifile_dat,*)iout,ilog,icase,ipred
       read(ifile_dat,*)garbage
       read(ifile_dat,*)eps,tol,imax

!       read the initial concentrations
       read(ifile_dat,*)(npar(i),i=1,iconc+2)
        do i1=1,icnd
        read(ifile_dat,*)(x(i),i=(i1-1)*(iconc+2)+1,i1*(iconc+2))
        read(ifile_dat,*)(ifit(i),i=(i1-1)*(iconc+2)+1,i1*(iconc+2))
         if (i1.gt.1) then
          do i2=1,iconc+2
           npar((i1-1)*(iconc+2)+i2)=npar(i2)
          enddo
         endif
        enddo
       read(ifile_dat,*)(ipx(i),i=1,iconc+1)
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

!       ixcm is the x array index of the last concentration
        ixcm = icnd*(iconc+2)

!       read in reaction mechanism
       read(ifile_dat,*)garbage
       do i=1,irx
         read(ifile_dat,*)npar(ixcm+2*i-1),kf(i),ifit(ixcm+2*i-1),ktfp(i), &
                  ktfx(i),icrf(i),(icrfi(i,i1),i1=1,icrf(i))
         read(ifile_dat,*)npar(ixcm+2*i),kr(i),ifit(ixcm+2*i),ktrp(i),ktrx(i), &
                  icrr(i),(icrri(i,i1),i1=1,icrr(i))
         x(ixcm+2*i-1) = kf(i)
         x(ixcm+2*i) = kr(i)
       enddo
       do i=1,imx
         x0(i) = x(i)
       enddo

!       read in experimental concentration data: species, time, value, unc
       read(ifile_dat,*)garbage
       sumunc = 0d0
       do i=1,iy
         read(ifile_dat,*)yc(i),ya(i),yt(i),y(i),unc(i)
         sumunc=sumunc+1d0/(unc(i)*unc(i))
       enddo
!      the weights wt(i) to be used in determining the sum of squared deviations
       do i=1,iy
         wt(i) = dfloat(iy)/(sumunc*unc(i)*unc(i))
       enddo
       
           close(ifile_dat)
!       end reading from input file

!       ix is the total number of parameters for this run
!       ix = 2*irx + ixcm    = fwd & rev rate constants per reaction + icnd x
!                                  (temp + t0 + initial concentrations) 
       ix = 2*irx + ixcm     

        if (ipred.eq.1) then
!         call res(iy,ix,imy,imx,y,x,omc,ssq)
          call fcn(iy,ix,x,y,ier)
         stop
        endif

!       find number of fitted parameters
! also, lmdif1 fits all the values of array x, so separate the fitted from fixed values
        itfit = 0
        do i=1,ix
          if (ifit(i)==1) then
            itfit = itfit + 1
            xfit(itfit) = x(i)
          endif
        enddo

        write(6,*)' Least squares fitting.'
                write(6,*)' '

        allocate(xjac(iy,itfit))

        ier = 0
! The original Press Numerical Recipes call
!        call xmrqmin(y,omc,unc,x,xjac,ssq,eps,tol,
!     2               ifit,ix,iy,imx,imy,imax,ier)
!  lwa is the length of the work array wa in lmdif1
       lwa = itfit*iy + 5*itfit + iy
       call lmdif1(fcn, iy, itfit, xfit, fvec, tol, info, iwa, wa, lwa)
! fvec returned by lmdif1 is actually the weighted omc array
       do i=1,iy
         omc(i) = fvec(i)/dsqrt(wt(i))
       enddo

        open(newunit=ifile_log,file='rate.log',status="replace")
        write(ifile_log,900)header
        stddev = 0d0
        if (iy > itfit+1) stddev = dsqrt(ssq/dfloat(iy-itfit-1))
        write(ifile_log,920)ssq,ier,stddev
! Error message statements
        if((ier.eq.129).or.(ier.eq.132))write(ifile_log,*) &
        ' Singularity in Jacobean.'
        if((ier.eq.133).or.(ier.eq.4))write(ifile_log,*)' imax exceeded.'
        if(ier.eq.39)write(ifile_log,*)' Convergence criteria too stringent.'
        if(ier.eq.38)write(ifile_log,*)' Zero Jacobean.'
        if(ier.eq.7)write(ifile_log,*)' Convergence too slow.'
        if(ier.eq.6)write(ifile_log,*)' Five steps taken with max step length.'
        if(ier.eq.2) &
        write(ifile_log,*)' Convergence to non-critical point suspected.'
        if(ier.eq.0)write(ifile_log,*)' Successful convergence.'
        if(ier.ne.0)write(ifile_log,919)
        write(ifile_log,949)
        do i=1,irx
         write(ifile_log,*)" "
         if ((icrf(i).eq.1).and.(icrr(i).eq.1)) then
          write(ifile_log,951)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.1).and.(icrr(i).eq.2)) then
          write(ifile_log,952)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.1).and.(icrr(i).eq.3)) then
          write(ifile_log,953)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.2).and.(icrr(i).eq.1)) then
          write(ifile_log,954)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.2).and.(icrr(i).eq.2)) then
          write(ifile_log,955)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.2).and.(icrr(i).eq.3)) then
          write(ifile_log,956)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.3).and.(icrr(i).eq.1)) then
          write(ifile_log,957)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.3).and.(icrr(i).eq.2)) then
          write(ifile_log,958)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).eq.3).and.(icrr(i).eq.3)) then
          write(ifile_log,959)i,npar(ixcm+2*i-1), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i)),npar(ixcm+2*i)
         else if ((icrf(i).gt.3).or.(icrr(i).gt.3)) then
          write(ifile_log,950)i,npar(ixcm+2*i-1),npar(ixcm+2*i), &
          (npar(icrfi(i,i1)+2),i1=1,icrf(i)),  &
          (npar(icrri(i,i1)+2),i1=1,icrr(i))
         endif
        enddo
        write(ifile_log,*)' '
        write(ifile_log,*)' Fit criteria'
        write(ifile_log,*)' '
        write(ifile_log,921)iy,ix,itfit
        dt = (tmax-t0)/dfloat(it)
        write(ifile_log,922)it,dt
        write(ifile_log,923)eps,tol,imax
        write(ifile_log,960)
!       have to call res again to refresh the omc array, which otherwise
!        may not correspond to the final ssq value
!       MINPACK's LMDIF1 requires fcn to be called in the format
!        call fcn(m, n, x, fvec, iflag)
!       where m = # observables, n = # parameters to be fitted, 
!         x = array of the parameters to be fitted, fvec = array of calculated function values
!        call res(iy,ix,imy,imx,y,x,omc,ssq)
        call fcn(iy,itfit,xfit,fvec,ier)

       ssq = 0d0
        do i=1,iy
         omc(i) = fvec(i)/dsqrt(wt(i))
         ycalc = y(i)-omc(i)
         write(ifile_log,961)yc(i),npar(ya(i)+2),yt(i),y(i),ycalc,omc(i),unc(i)
! try weighting ssq by the experimental unceratinties
         ssq = ssq + fvec(i)*fvec(i) 
        enddo
       stddev = 0d0
       if (iy > itfit+1) then
         stddev = dsqrt(ssq/dfloat(iy-itfit-1))
         stddvn = dsqrt(ssq*sumunc/(dfloat(iy-itfit-1)*dfloat(iy)))
       endif

       write(6,924)ssq,stddev,stddvn

        write(ifile_log,962)
        do i=1,ix
         if (i.le.ixcm) write(ifile_log,963)npar(mod(i-1,iconc+2)+1),x0(i),x(i)
         if (i.gt.ixcm) write(ifile_log,963)npar(i),x0(i),x(i)
        enddo

        ier = 0
        call fdjac2(fcn, iy, itfit, xfit, fvec, xjac, iy, ier, eps, wa)

! From this point forward, we only care about the fitted values of x and corresponing npar
        i2=1
        do i1=1,ix
         if (ifit(i1).eq.1) then
          npar(i2)=npar(i1)
          i2=i2+1
         endif
        enddo

        ist = 0

        call uncert(stddev,imx,imy,xjac,xfit,iy,itfit,ist,ifile_log)

        close(ifile_log)

900     format(a80)
919     format(' WARNING: Failure to converge with Press least ',/, &
        'squares routine does NOT necessarily report best values.')
920     format('  ssq = ',g9.3,5x,'  ier = ',i4,/,'  std dev = ',g10.4)
921     format('  # observables =',i4,', # variables =',i3, &
        ', # fitted variables = ',i3)
922     format('  # time steps =',i9,', delta t =',g12.6,' (linear t)')
923     format('  Convergence criteria: eps =',g8.2,3x,'tol =', &
        g8.2,3x,'imax =',i5,'.')
924     format('  ssq = ',g10.4,5x,'std dev = ',g10.4,5x,'norm std dev = ',g10.4)
949     format(/,' reaction mechanism')
950     format('  reaction #',i2,'  rate constants ',a3,', ',a3,/, &
        ' reactants: ',5(a3,1x),/,'  products: ',5(a3,1x))
951     format('  reaction #',i2,6x,a3,/,3x, &
                      12x,a3,' <---> ',a3,/,20x,a3)
952     format('  reaction #',i2,6x,a3,/,3x, &
                      12x,a3,' <---> ',a3,' + ',a3,/,20x,a3)
953     format('  reaction #',i2,6x,a3,/,3x, &
                      12x,a3,' <---> ',a3,' + ',a3,' + ',a3,/,20x,a3)
954     format('  reaction #',i2,6x,a3,/,3x, &
              6x,a3,' + ',a3,' <---> ',a3,/,20x,a3)
955     format('  reaction #',i2,6x,a3,/,3x, &
              6x,a3,' + ',a3,' <---> ',a3,' + ',a3,/,20x,a3)
956     format('  reaction #',i2,6x,a3,/,3x, &
              6x,a3,' + ',a3,' <---> ',a3,' + ',a3,' + ',a3,/,20x,a3)
957     format('  reaction #',i2,6x,a3,/,3x, &
        a3,' + ',a3,' + ',a3,' <---> ',a3,/,20x,a3)
958     format('  reaction #',i2,6x,a3,/,3x, &
        a3,' + ',a3,' + ',a3,' <---> ',a3,' + ',a3,/,20x,a3)
959     format('  reaction #',i2,6x,a3,/,3x, &
        a3,' + ',a3,' + ',a3,' <---> ',a3,' + ',a3,' + ',a3,/,20x,a3)
960     format('   species',/,'conds',6x,'t',12x,'y',12x,'ycalc',7x, &
        'omc',10x,'unc')
961     format(1x,i2,1x,a3,3(1x,g12.4),2(1x,g12.4))
962     format(/,' initial and final parameter values')
963     format(1x,a3,2(4x,g12.5))

       end program main

