       module residuals_module
!       use common_ratea
!       use common_rate
!       use common_rates
!       use common_temp
       use common_module

       implicit none

       real(dp) :: v(40)

       contains

!       subroutine res(iy,ix,imy,imx,y,x,omc,ssq)
       subroutine fcn(iy,itfit,xfit,fvec,ier)
!       ix is the total number of parameters for this run
!       itfit is the total number of *fitted* parameters for this run
!       iy is the total number of data points for this run
!       imx, imy are the array dimensions for x and y

       real(dp) :: x2(28),x3(28),x4(28),xa(28),xb(28),xc(28),xd(28),xp(28)=0.0_dp
       real(dp) :: delx(40)=0.0_dp,delx2(40)=0.0_dp,delx3(40)=0.0_dp,delx4(40)=0.0_dp,dlogt=0.0_dp
       real(dp) :: ycalc(2100),stddev,stddvn,ssq,t0,tnew,tlmin,xt(28)=0.0_dp
       integer :: ichk(2100),itime(2100)=0
       integer :: i,i0,i1,i2,i3,imx,imy,iptot,ixcm,ikf,ikr,index=0
       integer :: ifile_out,ifile_out1,ifile_out2,ifile_out3,ifile_out4, &
         ifile_out5,ifile_out6
        character(len=3) :: xpar(108),cpar(108)
       integer, intent(in) :: iy,itfit
       real(dp), intent(in) :: xfit(ix)
       real(dp), intent(out) :: fvec(iy)
       integer, intent(inout) :: ier

!       ichk is a flag to ensure that all data are involved in fit
        do i=1,iy
         ichk(i) = 0
        enddo

!       transfer values of xfit into the x array
        i1=0
        do i=1,ix
          if (ifit(i) == 1) then
            i1 = i1 + 1
            x(i) = xfit(i1)
          endif
        enddo

!       set up the output table 
       if(iout == 0) then
         open(newunit=ifile_out,file='rate.out',status="replace")
       else if(iout >= 1) then
         iptot=0
         do i=1,iconc+1
           if(ipx(i) == 1)then
             iptot=iptot+1
           endif
         enddo
         if(iout == 1) then
           open(newunit=ifile_out,file='rate.out',status="replace")
         else if(iout == 3) then
           open(newunit=ifile_out,file='rate.out',status="replace")
         else if(iout == 2) then
           if(iptot > 5) open(newunit=ifile_out6,file='rate6.out',status="replace")
           if(iptot > 4) open(newunit=ifile_out5,file='rate5.out',status="replace")
           if(iptot > 3) open(newunit=ifile_out4,file='rate4.out',status="replace")
           if(iptot > 2) open(newunit=ifile_out3,file='rate3.out',status="replace")
           if(iptot > 1) open(newunit=ifile_out2,file='rate2.out',status="replace")
           open(newunit=ifile_out1,file='rate1.out',status="replace")
         endif
       endif

!       convert x array to k's.  
        ixcm = icnd*(iconc+2)
        do i=1,irx
         ikr = ixcm+2*i
         ikf = ikr-1
         kf(i) = x(ikf)
         kr(i) = x(ikr)
        enddo

!       Begin loop over initial conditions
        do i0 = 1,icnd

!       convert x array to current temp, t0, initial concs xa.
!        xa, temp, and t0 are all local to the current set of conditions.
        temp = x((i0-1)*(iconc+2)+1)
        t0 = x((i0-1)*(iconc+2)+2)
        do i=1,iconc
         xa(i) = x((i0-1)*(iconc+2)+2+i) 
        enddo

!       establish time increments, which may vary with condition index i0 if t0 varies
       if(ilog == 0) dt = (tmax-t0)/dfloat(it)
       if(ilog == 1) dlogt = (dlog10(tmax-t0)-dlog10(tlmin))/dfloat(it)

!       write initial concentrations
       t=t0
       xa(iconc+1) = 0d0
       do i=1,iconc
         xa(iconc+1) = xa(iconc+1)+xa(i)
         cpar(i) = npar(i+2)
       enddo
!       note that xa(iconc+1) is the total concentration of all species
       cpar(iconc+1) = 'tot'

!       get the output table started
       if(iout == 0) then
         write(ifile_out,900)t
         do i=1,iconc
           write(ifile_out,905)cpar(i),xa(i)
         enddo
         write(ifile_out,905)cpar(iconc+1),xa(iconc+1)
       else if(iout >= 1) then
         index=0
         do i=1,iconc+1
           if(ipx(i) == 1)then
             index=index+1
             xp(index)=xa(i)
             xpar(index)=cpar(i)
           endif
         enddo
         if(iout == 1) then
           if(iptot > 5) then
             write(ifile_out,950)(xpar(i),i=1,iptot)
             write(ifile_out,955)t,(xp(i),i=1,iptot)
           else if(iptot <= 5) then
             write(ifile_out,951)(xpar(i),i=1,iptot)
             write(ifile_out,956)t,(xp(i),i=1,iptot)
           endif
         else if(iout == 3) then
           if(iptot > 5) then
             write(ifile_out,952)(xpar(i),i=1,iptot)
             write(ifile_out,957)t,(xp(i),i=1,iptot)
           else if(iptot <= 5) then
             write(ifile_out,953)(xpar(i),i=1,iptot)
             write(ifile_out,958)t,(xp(i),i=1,iptot)
           endif
         else if(iout == 2) then
           if(iptot > 5) then
             write(ifile_out6,951)xpar(6)
             write(ifile_out6,956)t,xp(6)
           endif
           if(iptot > 4) then
             write(ifile_out5,951)xpar(5)
             write(ifile_out5,956)t,xp(5)
           endif
           if(iptot > 3) then
             write(ifile_out4,951)xpar(4)
             write(ifile_out4,956)t,xp(4)
           endif
           if(iptot > 2) then
             write(ifile_out3,951)xpar(3)
             write(ifile_out3,956)t,xp(3)
           endif
           if(iptot > 1) then
             write(ifile_out2,951)xpar(2)
             write(ifile_out2,956)t,xp(2)
           endif
           write(ifile_out1,951)xpar(1)
           write(ifile_out1,956)t,xp(1)
         endif
       endif

!       determine the values of i that correspond to observed values
       do i=1,iy
        if (yc(i) == i0) then
          if (ilog == 0) itime(i) = int((yt(i)-t0)/dt) + 1
          if (ilog == 1) itime(i) =  &
            int( (dlog10(yt(i)-t0) - dlog10(tlmin))/dlogt) + 1
          if (itime(i) > it) then
            write(6,*)" Reset tmax in input file to be greater than latest time in data set.  Exiting"
            stop
          endif
         endif
       enddo

!       Begin loop over time
!       establish time t, which is measured from time t0, not necessarily the same as yt
       t = t0
       do i1 = 1,it
          do i=1,27
           xt(i)=0d0
          enddo
         if(ilog == 0) t = t + dt
         if(ilog == 1) then
           tnew = 10**(dlog10(tlmin)+dfloat(i1)*dlogt)
           dt = tnew-t
           t=tnew
         endif

!         Integrate kinetic equations
         call integrate(xa,xt,delx)

!         Add fourth order Runge-Kutta method corrections if called for
         if(icase == 1)then
           do i=1,iconc
             xb(i) = xa(i) + delx(i)/2d0
           enddo
           call integrate(xb,x2,delx2)
           do i=1,iconc
             xc(i) = xa(i) + delx2(i)/2d0
           enddo
           call integrate(xc,x3,delx3)
           do i=1,iconc
             xd(i) = xa(i) + delx3(i)
           enddo
           call integrate(xd,x4,delx4)
           do i=1,iconc
             xt(i) = xa(i) + (delx(i) + 2d0*delx2(i) + 2d0*delx3(i) &
                          + delx4(i))/6d0
           enddo
         endif
           
!         Move new concentrations into old array
         xa(iconc+1) = 0d0
         do i=1,iconc
           if (xt(i) >= 0d0) xa(i)=xt(i)
           if (xt(i) < 0d0) xa(i)=0d0
           xa(iconc+1) = xa(iconc+1)+xa(i)
         enddo

!         Report the new concentrations if at the right interval
         if(mod(i1,ipt) == 0)then
           if(iout == 0) then
             write(ifile_out,900)t
             do i=1,iconc
               if((log(xa(i)) < 9).and.(log(xa(i)) > -10)) &
                 write(ifile_out,905)cpar(i),xa(i)
               if((log(xa(i)) >= 9).or.(log(xa(i)) <= -10)) &
                 write(ifile_out,906)cpar(i),xa(i)
             enddo
             if((log(xa(iconc+1)) < 9).and.(log(xa(iconc+1)) > -10)) &
               write(ifile_out,905)cpar(iconc+1),xa(iconc+1)
             if((log(xa(iconc+1)) >= 9).or.(log(xa(iconc+1)) <= -10)) &
               write(ifile_out,906)cpar(iconc+1),xa(iconc+1)
           else if(iout >= 1) then
             index=0
             do i=1,iconc+1
               if(ipx(i) == 1)then
                 index=index+1
                 xp(index)=xa(i)
               endif
             enddo
             if(iout == 1) then
               if(mod(i1,ipt*20) == 0)then
                 if(iptot > 5) write(ifile_out,950)(xpar(i),i=1,iptot)
                 if(iptot <= 5) write(ifile_out,951)(xpar(i),i=1,iptot)
               endif
               if(iptot > 5)write(ifile_out,955)t,(xp(i),i=1,iptot)
               if(iptot <= 5)write(ifile_out,956)t,(xp(i),i=1,iptot)
             else if(iout == 2) then
               if(iptot > 5) write(ifile_out6,956)t,xp(6)
               if(iptot > 4) write(ifile_out5,956)t,xp(5)
               if(iptot > 3) write(ifile_out4,956)t,xp(4)
               if(iptot > 2) write(ifile_out3,956)t,xp(3)
               if(iptot > 1) write(ifile_out2,956)t,xp(2)
               write(ifile_out1,956)t,xp(1)
             else if(iout == 3) then
               if(iptot > 5)write(ifile_out,957)t,(xp(i),i=1,iptot)
               if(iptot <= 5)write(ifile_out,958)t,(xp(i),i=1,iptot)
             endif
           endif
         endif

!       update ycalc values if any are at this time and set of conditions
         do i2=1,iy
           if ((yc(i2) == i0).and.(i1 == itime(i2))) then
              ycalc(i2) = xa(ya(i2))
              ichk(i2) = 1
            endif
         enddo
       enddo
!       End loop over time

        enddo
!       End loop over initial conditions

!       verify that all data were used
        do i2=1,iy
         if (ichk(i2).ne.1) then
          write(6,*) "datum ",i2," skipped somehow; exiting"
          stop
         endif
        enddo

! weight ssq by the experimental unceratinties
       ssq = 0d0
       do i2=1,iy
         fvec(i2) = dsqrt(wt(i2)) * (y(i2) - ycalc(i2))
         ssq = ssq + fvec(i2)*fvec(i2) 
       enddo
       stddev = 0d0
       if (iy > itfit+1) then
         stddev = dsqrt(ssq/dfloat(iy-itfit-1))
         stddvn = dsqrt(ssq*sumunc/(dfloat(iy-itfit-1)*dfloat(iy)))
       endif
       write(6,940)ssq,stddev,stddvn

900       format('t = ',f20.10)
905       format(2x,a7,2x,f20.10)
906       format(2x,a7,2x,e20.10)
940       format("ssq = ",g10.4,5x,"std dev = ",g10.4,5x,"norm std dev = ",g10.4)
950       format(8x,8(2x,a7))
951       format(12x,8(6x,a7))
952       format('t',7x,8(',',1x,a7))
953       format('t',11x,8(',',5x,a7))
955       format(e8.2,8(1x,e8.2))
956       format(e12.6,5(1x,e12.6))
957       format(e8.2,8(',',e8.2))
958       format(e12.6,5(',',e12.6))
990       format('Negative concentration detected for ',a7,'. Exiting.')

1000       close(unit=4)
       close(unit=ifile_out)
       if(iout == 2) then
        if(iptot > 5) close(ifile_out6)
        if(iptot > 4) close(ifile_out5)
        if(iptot > 3) close(ifile_out4)
        if(iptot > 2) close(ifile_out3)
        if(iptot > 1) close(ifile_out2)
        close(ifile_out1)
        endif

        return
       end subroutine fcn
!       end subroutine res



       subroutine integrate(xa,xt,delx)
!       Calculates the change in concentration for each component
!        as Delta [X] = (delta t) * v
       real(dp) :: dx=0.0_dp
       real(dp),intent(in) :: xa(28)
       real(dp),intent(out) :: xt(28),delx(40)
       integer :: i1,i2,i3

        do i1=1,irx
         v(i1)=0d0
        enddo

!       get the current reaction velocity v = k[A][B]...
!        for each reaction based on current concentrations xa
       call velocity(xa,v)

!       loop over all chemical components
       do i1=1,iconc
         dx = 0d0
!         loop over all reactions
         do i2=1,irx
!           loop over all reactants for forward, then reverse
!            but make sure not to count any reactant twice
            i3=1
            do while (i3 <= icrf(i2))
              if (icrfi(i2,i3)==i1) then
                dx = dx - v(i2)
                exit
              endif
              i3 = i3 + 1
            enddo
            i3=1
            do while (i3 <= icrr(i2))
              if (icrri(i2,i3)==i1) then
                dx = dx + v(i2)
                exit
              endif
	      i3 = i3 + 1
            enddo
         enddo
  
         delx(i1) = dx*dt
         xt(i1) = xa(i1) + delx(i1)
       enddo

       return
       end subroutine integrate

       subroutine velocity(xa,v)
       real(dp) :: for,rev,tp,tx
       integer :: ipx(28),i1,i2,index
       real(dp),intent(in) :: xa(28)
       real(dp),intent(out) :: v(40)

       do i1=1,irx
!         Forward half-reaction
         tx=1d0
         if(ktfx(i1).ne.0d0)tx=exp(-ktfx(i1)/temp)
         tp=1d0
         if(ktfp(i1).ne.0d0)tp=temp**ktfp(i1)
         for = kf(i1)*tp*tx
         do i2=1,icrf(i1)
           for = for*xa(icrfi(i1,i2))
         enddo

!         Reverse half-reaction
         tx=1d0
         if(ktrx(i1).ne.0d0)tx=exp(-ktrx(i1)/temp)
         tp=1d0
         if(ktrp(i1).ne.0d0)tp=temp**ktrp(i1)
         rev = kr(i1)*tp*tx
         do i2=1,icrr(i1)
           rev = rev*xa(icrri(i1,i2))
         enddo
         v(i1) = for - rev
       enddo

       return
       end subroutine velocity

       end module residuals_module
