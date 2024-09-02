       module common_module
       implicit none

       integer, parameter :: dp=kind(0.d0)

       character(len=3) :: npar(108)

       integer :: imrx,imconc,imy,imx
       integer :: yc(2100),ya(2100)
       integer :: icrf(40),icrr(40),icrfi(40,5),icrri(40,5),ipx(28), &
           iconc,irx,icnd,ix,it,ipt,iout,ilog,icase,itfit,ifit(108)
       real(dp) :: y(2100),yt(2100),wt(2100),unc(2100),sumunc
       real(dp) :: t,dt,ktfx(40),ktrx(40),ktfp(40),ktrp(40),tmin,tmax
       real(dp) :: x(108)
       real(dp) :: temp,kf(40),kr(40)

       data imrx,imconc,imy,imx/40,27,2100,108/

       end module common_module

