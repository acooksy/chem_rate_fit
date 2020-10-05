	subroutine res(iy,ix,imy,imx,y,x,omc,ssq)
c       ix is the total number of parameters for this run
c       iy is the total number of data points for this run
c       imx, imy are the array dimensions for x and y

        implicit double precision (a-h,j-z)
	real*8 x2(28),x3(28),x4(28),xb(28),xc(28),xd(28),xp(28),xt(28),
     2    xa(28),delx(40),delx2(40),delx3(40),delx4(40),
     2    kf(40),kr(40),ktfx(40),ktrp(40),ktrx(40),ktfp(40),
     3    x(108),omc(800),y(800),ycalc(800),yt(800)
	integer*4 ipx(28),icrf(40),icrr(40),icrfi(40,5),icrri(40,5),
     2    ichk(800),itime(800),ya(800),yc(800)
        character*3 xpar(108),npar(108),cpar(108)
        common/rates/t,dt,ktfx,ktrx,ktfp,ktrp,tlmin,tmax,
     2    iconc,irx,icnd,icrf,icrr,icrfi,icrri,it,ipt,iout,ilog,
     3    icase,ipx,itfit
	common/rate/yc,ya,yt
	common/ratea/npar
        common/temp/temp,kf,kr

c       ichk is a flag to ensure that all data are involved in fit
        do i=1,iy
         ichk(i) = 0
        enddo

c	Get the output table set up
	if(iout.eq.0) then
	  open(unit=8,file='rate.out')
	else if(iout.ge.1) then
	  iptot=0
	  do i=1,iconc+1
	    if(ipx(i).eq.1)then
	      iptot=iptot+1
	    endif
	  enddo
	  if(iout.eq.1) then
	    open(unit=8,file='rate.out')
	  else if(iout.eq.3) then
	    open(unit=8,file='rate.out')
	  else if(iout.eq.2) then
	    if(iptot.gt.5) open(unit=1,file='rate6.out')
	    if(iptot.gt.4) open(unit=2,file='rate5.out')
	    if(iptot.gt.3) open(unit=3,file='rate4.out')
	    if(iptot.gt.2) open(unit=8,file='rate3.out')
	    if(iptot.gt.1) open(unit=9,file='rate2.out')
	    open(unit=10,file='rate1.out')
	  endif
	endif

c       convert x array to k's.  
        ixcm = icnd*(iconc+2)
        do i=1,irx
         ikr = ixcm+2*i
         ikf = ikr-1
         kf(i) = x(ikf)
         kr(i) = x(ikr)
        enddo

c       Begin loop over initial conditions
        do i0 = 1,icnd

c       convert x array to current temp, t0, initial concs xa.
c        xa, temp, and t0 are all local to the current set of conditions.
        temp = x((i0-1)*(iconc+2)+1)
        t0 = x((i0-1)*(iconc+2)+2)
        do i=1,iconc
         xa(i) = x((i0-1)*(iconc+2)+2+i) 
        enddo

c       establish time increments, which may vary with condition index i0 if t0 varies
	if(ilog.eq.0) dt = (tmax-t0)/dfloat(it)
	if(ilog.eq.1) dlogt = (dlog10(tmax-t0)-dlog10(tlmin))/dfloat(it)

c	write initial concentrations
	t=t0
	xa(iconc+1) = 0d0
	do i=1,iconc
	  xa(iconc+1) = xa(iconc+1)+xa(i)
	  cpar(i) = npar(i+2)
	enddo
c	note that xa(iconc+1) is the total concentration of all species
	cpar(iconc+1) = 'tot'

c	get the output table started
	if(iout.eq.0) then
	  write(8,900)t
	  do i=1,iconc
	    write(8,905)cpar(i),xa(i)
	  enddo
	  write(8,905)cpar(iconc+1),xa(iconc+1)
	else if(iout.ge.1) then
	  index=0
	  do i=1,iconc+1
	    if(ipx(i).eq.1)then
	      index=index+1
	      xp(index)=xa(i)
	      xpar(index)=cpar(i)
	    endif
	  enddo
	  if(iout.eq.1) then
	    if(iptot.gt.5) then
	      write(8,950)(xpar(i),i=1,iptot)
	      write(8,955)t,(xp(i),i=1,iptot)
	    else if(iptot.le.5) then
	      write(8,951)(xpar(i),i=1,iptot)
	      write(8,956)t,(xp(i),i=1,iptot)
	    endif
	  else if(iout.eq.3) then
	    if(iptot.gt.5) then
	      write(8,952)(xpar(i),i=1,iptot)
	      write(8,957)t,(xp(i),i=1,iptot)
	    else if(iptot.le.5) then
	      write(8,953)(xpar(i),i=1,iptot)
	      write(8,958)t,(xp(i),i=1,iptot)
	    endif
	  else if(iout.eq.2) then
	    if(iptot.gt.5) then
	      write(1,951)xpar(6)
	      write(1,956)t,xp(6)
	    endif
	    if(iptot.gt.4) then
	      write(2,951)xpar(5)
	      write(2,956)t,xp(5)
	    endif
	    if(iptot.gt.3) then
	      write(3,951)xpar(4)
	      write(3,956)t,xp(4)
	    endif
	    if(iptot.gt.2) then
	      write(8,951)xpar(3)
	      write(8,956)t,xp(3)
	    endif
	    if(iptot.gt.1) then
	      write(9,951)xpar(2)
	      write(9,956)t,xp(2)
	    endif
	    write(10,951)xpar(1)
	    write(10,956)t,xp(1)
	  endif
	endif

c	determine the values of i that correspond to observed values
	do i=1,iy
	 if (yc(i).eq.i0) then
          if (ilog.eq.0) itime(i) = int((yt(i)-t0)/dt) + 1
          if (ilog.eq.1) itime(i) = 
     2       int( (dlog10(yt(i)-t0) - dlog10(tlmin))/dlogt) + 1
         endif
	enddo

c	Begin loop over time
c	establish time t, which is measured from time t0, not necessarily the same as yt
	t = t0
	do i1 = 1,it
          do i=1,27
           xt(i)=0d0
          enddo
	  if(ilog.eq.0) t = t + dt
	  if(ilog.eq.1) then
	    tnew = 10**(dlog10(tlmin)+dfloat(i1)*dlogt)
	    dt = tnew-t
	    t=tnew
	  endif

c	  Integrate kinetic equations
	  call integrate(xa,xt,delx)

c	  Add fourth order Runge-Kutta method corrections if called for
	  if(icase.eq.1)then
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
	      xt(i) = xa(i) + (delx(i) + 2d0*delx2(i) + 2d0*delx3(i)
     2                     + delx4(i))/6d0
	    enddo
	  endif
	    
c	  Move new concentrations into old array
	  xa(iconc+1) = 0d0
	  do i=1,iconc
	    if (xt(i).ge.0d0) xa(i)=xt(i)
	    if (xt(i).lt.0d0) xa(i)=0d0
	    xa(iconc+1) = xa(iconc+1)+xa(i)
	  enddo

c	  Report the new concentrations if at the right interval
	  if(mod(i1,ipt).eq.0)then
	    if(iout.eq.0) then
	      write(8,900)t
	      do i=1,iconc
	        if((log(xa(i)).lt.9).and.(log(xa(i)).gt.-10))
     2            write(8,905)cpar(i),xa(i)
	        if((log(xa(i)).ge.9).or.(log(xa(i)).le.-10))
     2            write(8,906)cpar(i),xa(i)
	      enddo
	      if((log(xa(iconc+1)).lt.9).and.(log(xa(iconc+1)).gt.-10))
     2          write(8,905)cpar(iconc+1),xa(iconc+1)
	      if((log(xa(iconc+1)).ge.9).or.(log(xa(iconc+1)).le.-10))
     2          write(8,906)cpar(iconc+1),xa(iconc+1)
	    else if(iout.ge.1) then
	      index=0
	      do i=1,iconc+1
	        if(ipx(i).eq.1)then
	          index=index+1
	          xp(index)=xa(i)
	        endif
	      enddo
	      if(iout.eq.1) then
	        if(mod(i1,ipt*20).eq.0)then
	          if(iptot.gt.5) write(8,950)(xpar(i),i=1,iptot)
	          if(iptot.le.5) write(8,951)(xpar(i),i=1,iptot)
	        endif
	        if(iptot.gt.5)write(8,955)t,(xp(i),i=1,iptot)
	        if(iptot.le.5)write(8,956)t,(xp(i),i=1,iptot)
	      else if(iout.eq.2) then
	        if(iptot.gt.5) write(1,956)t,xp(6)
	        if(iptot.gt.4) write(2,956)t,xp(5)
	        if(iptot.gt.3) write(3,956)t,xp(4)
	        if(iptot.gt.2) write(8,956)t,xp(3)
	        if(iptot.gt.1) write(9,956)t,xp(2)
	        write(10,956)t,xp(1)
	      else if(iout.eq.3) then
	        if(iptot.gt.5)write(8,957)t,(xp(i),i=1,iptot)
	        if(iptot.le.5)write(8,958)t,(xp(i),i=1,iptot)
	      endif
	    endif
	  endif

c	update ycalc values if any are at this time
	  do i2=1,iy
	    if ((yc(i2).eq.i0).and.(i1.eq.itime(i2))) then
              ycalc(i2) = xa(ya(i2))
              ichk(i2) = 1
            endif
	  enddo
	enddo
c       End loop over time

        enddo
c       End loop over initial conditions

c       verify that all data were used
        do i2=1,iy
         if (ichk(i2).ne.1) then
          write(6,*) "datum ",i2," skipped somehow; exiting"
          stop
         endif
        enddo

	ssq = 0d0
	do i2=1,iy
	  omc(i2) = y(i2) - ycalc(i2)
	  ssq = ssq + omc(i2)*omc(i2)
	enddo
        stddev = 0d0
        if (iy.gt.itfit+1) stddev = dsqrt(ssq/dfloat(iy-itfit-1))
        write(6,*)"ssq = ",ssq,"    std dev = ",stddev

900	format('t = ',f20.10)
905	format(2x,a7,2x,f20.10)
906	format(2x,a7,2x,e20.10)
950	format(8x,8(2x,a7))
951	format(12x,8(6x,a7))
952	format('t',7x,8(',',1x,a7))
953	format('t',11x,8(',',5x,a7))
955	format(e8.2,8(1x,e8.2))
956	format(e12.6,5(1x,e12.6))
957	format(e8.2,8(',',e8.2))
958	format(e12.6,5(',',e12.6))
990	format('Negative concentration detected for ',a7,'.  Exiting.')

1000	close(unit=4)
	close(unit=8)
	if(iout.eq.2) then
	 if(iptot.gt.5) close(unit=1)
	 if(iptot.gt.4) close(unit=2)
	 if(iptot.gt.3) close(unit=3)
	 if(iptot.gt.2) close(unit=8)
	 if(iptot.gt.1) close(unit=9)
	 close(unit=10)
        endif

        return
	end



	subroutine integrate(xa,xt,delx)
c       Calculates the change in concentration for each component
c        as Delta [X] = (delta t) * v
	implicit double precision(a-h,j-z)
	real*8 xt(28),xa(28),
     2    kf(40),kr(40),ktfx(40),ktrx(40),ktfp(40),ktrp(40),
     3    delx(40),v(40)
	integer*4 ipx(28),icrf(40),icrr(40),icrfi(40,5),icrri(40,5)
	common/rates/t,dt,ktfx,ktrx,ktfp,ktrp,tlmin,tmax,
     2    iconc,irx,icnd,icrf,icrr,icrfi,icrri,it,ipt,iout,ilog,
     3    icase,ipx,itfit

        do i1=1,irx
         v(i1)=0d0
        enddo

c       get the current reaction velocity v = k[A][B]...
c        for each reaction based on current concentrations xa
	call velocity(xa,v)

c       loop over all chemical components
	do i1=1,iconc
	  dx = 0d0
c         loop over all reactions
	  do 250 i2=1,irx
c           loop over all reactants for forward, then reverse
c            but make sure not to count any reactant twice
            do 100 i3=1,icrf(i2)
              if(icrfi(i2,i3).eq.i1) then
                dx = dx - v(i2)
                goto150
              endif
100         continue
150         do 200 i3=1,icrr(i2)
              if(icrri(i2,i3).eq.i1) then
                dx = dx + v(i2)
                goto250
              endif
200         continue
250       continue
  
	  delx(i1) = dx*dt
	  xt(i1) = xa(i1) + delx(i1)
	enddo

	return
	end


	subroutine velocity(xa,v)
	implicit double precision(a-h,j-z)
	real*8 xa(28),
     2    kf(40),kr(40),ktfx(40),ktrx(40),ktfp(40),ktrp(40),v(40)
	integer*4 ipx(28),icrf(40),icrr(40),icrfi(40,5),icrri(40,5)
	common/rates/t,dt,ktfx,ktrx,ktfp,ktrp,tlmin,tmax,
     2    iconc,irx,icnd,icrf,icrr,icrfi,icrri,it,ipt,iout,ilog,
     3    icase,ipx,itfit
        common/temp/temp,kf,kr

	do i1=1,irx
c	  Forward half-reaction
	  tx=1d0
	  if(ktfx(i1).ne.0d0)tx=exp(-ktfx(i1)/temp)
	  tp=1d0
	  if(ktfp(i1).ne.0d0)tp=temp**ktfp(i1)
	  for = kf(i1)*tp*tx
	  do i2=1,icrf(i1)
	    for = for*xa(icrfi(i1,i2))
	  enddo

c	  Reverse half-reaction
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
	end

