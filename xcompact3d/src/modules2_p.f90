  module tools_p
!********************************

! Module for various tools for processing data etc
implicit none
contains

subroutine conv_p(itime,dt,cutOff,ux1,uy1,uz1,uxR1_p,uyR1_p,uzR1_p,ta1)
!#########################

use decomp_2d
use decomp_2d_io
use param_p, only: foldName_p
use var_p, only: dummyString
use MPI

implicit none

real(mytype), intent(in) :: cutOff
real(mytype),dimension(:,:,:) :: ux1,uy1,uz1,uxR1_p,uyR1_p,uzR1_p,ta1
integer :: mpiErr,itime
real(mytype) :: maxErr,totErr,mETemp,tETemp,dt

if(cutOff .ge. 0d0) then   
  ta1 = sqrt((ux1-uxR1_p)**2+(uy1-uyR1_p)**2 + (uzR1_p-uz1)**2) ! uxR1_p no need any longer
  
  mETemp = maxval(ta1); tETemp = sum(ta1)
  
  call MPI_REDUCE(mETemp,maxErr,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,mpiErr)
  call MPI_REDUCE(tETemp,totErr,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)
 
  uxR1_p = ux1; uyR1_p = uy1 ; uzR1_p = uz1 ! Replace references
 

  if(nrank == 0) then
      open(163,file=trim(foldName_p)//'maxErr.txt',status='unknown',position='append')
      write(163,*) itime*dt,' ',maxErr
      close(163)
      open(163,file=trim(foldName_p)//'totErr.txt',status='unknown',position='append')
      write(163,*) itime*dt,' ',totErr
      close(163)
      if(maxErr .le. cutOff) then
         write(*,*) 'Code converged with max diff betwn consecutive viz of: ', maxErr
         call decomp_2d_abort(163,'Code converged')
      endif
  endif
endif

!#########################
end subroutine conv_p

!@@10A start
!#########################
 subroutine streakline(ux1,uy1,uz1,itime,isave,halt)
!#########################

use decomp_2d
use param_p
use var_p
use param, only: dx,dy,dz,dt,xlx,yly,zlz
use variables
use mpi

implicit none


integer :: i,j,ierr,freezeIterEnd
integer :: ixMin, ixMax, iyMin, iyMax, izMin, izMax
integer :: iyMinR,iyMaxR,izMinR,izMaxR !Relative positions in pencil
integer, intent(in) :: itime,isave
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1,uy1,uz1
integer, intent(out) :: halt
real(mytype), dimension(2,2,2) :: uxAdj,uyAdj,uzAdj,uxAdj1,uyAdj1,uzAdj1
real(mytype) :: xd,yd,zd ! Check wiki page of trilinear interpoln for convention
real(mytype) :: uxInt,uyInt,uzInt ! Interpolated velocity at current position of particle

1234 format(ES15.8,a,ES15.8,a,ES15.8,a,ES15.8)
1235 format(I0.6,a,I0.6,a,I0.6,a,ES15.8,a,ES15.8,a,ES15.8)

freezeIterEnd = min((itime-1)/skipTime_p+1,noOfReleases_p)
!if skip = 1, for t=1,fI = 1, t=2, fI =2 etc till noOfRel exceeded  
!if skip = 2, for t=1,2:fI=1, t=3,4 fI=2 and so on

! Add the positions of frozen particles just before release in 1 file visualization
if(nrank == 0) then
   i = mod(itime,skipTime_p) !if skip=2,itime=2,4,.. for i=0
   if(i == 0) then  
      j = itime/skipTime_p+1   !That which is to be unfrozen in next iter
      if(j <= noOfReleases_p) then
         do i= 1,strkPosns_p
            if(strkSwitch_p .ne. 2) then
               write(666,1235) itime,'  ',i,'  ',j,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
            endif
            if(strkSwitch_p .ne. 1) then
               write(dummystring,'(a,i3.3,a,i0.3)') 'streaklines/path',i,'At',j
               open(10,file=trim(dummystring),status='replace')
               write(10,*) 'Time        X        Y        Z'
               write(10,1234) (j-1)*skipTime_p*dt,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
               close(10)
            endif
         enddo
      endif
   endif   
endif


! Check if for all positions, atleast one particle has exited, if so start bufferTime
halt = 0 ! Default
if(sum(strkExit) == strkPosns_p) then
   bufferTime_p = bufferTime_p +dt
   if(bufferTime_p >= (isave-2)*dt) then
      halt = 1
      return
   endif
endif


do i = 1,strkPosns_p
   do j = 1,freezeIterEnd

   ! Check if within domain, else take note for given position and skip rest
     if(strkX_p(i,j) .ge. xlx) then
        strkExit(i) = 1
        cycle
     endif
     if(strkY_p(i,j) .le. -yly/2d0) then
        strkExit(i) = 1
        cycle
     endif
     if(strkY_p(i,j) .ge.  yly/2d0) then
        strkExit(i) = 1
        cycle
     endif
     if(strkZ_p(i,j) .le. -zlz/2d0) then
        strkExit(i) = 1
        cycle
     endif
     if(strkZ_p(i,j) .ge.  zlz/2d0) then
        strkExit(i) = 1
        cycle
     endif
   
   ! Find adjacent points
      ixMin = floor(strkX_p(i,j)/dx)+1;       ixMax = ixMin+1
      iyMin = floor(strkY_p(i,j)/dy)+nym/2+1; iyMax = iyMin+1
      izMin = floor(strkZ_p(i,j)/dz)+nzm/2+1; izMax = izMin+1

      xd = (strkX_p(i,j)-(ixMin-1)*dx)/dx        ! (x-x0)/(x1-x0)
      yd = (strkY_p(i,j)-(iyMin-nym/2-1)*dy)/dy
      zd = (strkZ_p(i,j)-(izMin-nzm/2-1)*dz)/dz

      ! note: Modification for STRETCHING not implemented yet
       ! Relative position in pencil
      iyMinR = iyMin-xstart(2)+1
      iyMaxR = iyMax-xstart(2)+1 
      izMinR = izMin-xstart(3)+1 
      izMaxR = izMax-xstart(3)+1 

      uxAdj = 0d0; uyAdj = 0d0; uzAdj = 0d0 ! Make sure all are zero b4 summing

      if(iyMin .ge. xstart(2) .and. iyMin .le. xend(2)) then
         if(izMin .ge. xstart(3) .and. izMin .le. xend(3)) then
            uxAdj(1,1,1) = ux1(ixMin,iyMinR,izMinR)
            uyAdj(1,1,1) = uy1(ixMin,iyMinR,izMinR)
            uzAdj(1,1,1) = uz1(ixMin,iyMinR,izMinR)
            uxAdj(2,1,1) = ux1(ixMax,iyMinR,izMinR)
            uyAdj(2,1,1) = uy1(ixMax,iyMinR,izMinR)
            uzAdj(2,1,1) = uz1(ixMax,iyMinR,izMinR)
         endif 
         if(izMax .ge. xstart(3) .and. izMax .le. xend(3)) then
            uxAdj(1,1,2) = ux1(ixMin,iyMinR,izMaxR)
            uyAdj(1,1,2) = uy1(ixMin,iyMinR,izMaxR)
            uzAdj(1,1,2) = uz1(ixMin,iyMinR,izMaxR)
            uxAdj(2,1,2) = ux1(ixMax,iyMinR,izMaxR)
            uyAdj(2,1,2) = uy1(ixMax,iyMinR,izMaxR)
            uzAdj(2,1,2) = uz1(ixMax,iyMinR,izMaxR)
         endif
      endif
      if(iyMax .ge. xstart(2) .and. iyMax .le. xend(2)) then
        if(izMin .ge. xstart(3) .and. izMin .le. xend(3)) then
           uxAdj(1,2,1) = ux1(ixMin,iyMaxR,izMinR)
           uyAdj(1,2,1) = uy1(ixMin,iyMaxR,izMinR)
           uzAdj(1,2,1) = uz1(ixMin,iyMaxR,izMinR)
           uxAdj(2,2,1) = ux1(ixMax,iyMaxR,izMinR)
           uyAdj(2,2,1) = uy1(ixMax,iyMaxR,izMinR)
           uzAdj(2,2,1) = uz1(ixMax,iyMaxR,izMinR)
        endif
        if(izMax .ge. xstart(3) .and. izMax .le. xend(3)) then
           uxAdj(1,2,2) = ux1(ixMin,iyMaxR,izMaxR)
           uyAdj(1,2,2) = uy1(ixMin,iyMaxR,izMaxR)
           uzAdj(1,2,2) = uz1(ixMin,iyMaxR,izMaxR)
           uxAdj(2,2,2) = ux1(ixMax,iyMaxR,izMaxR)
           uyAdj(2,2,2) = uy1(ixMax,iyMaxR,izMaxR)
           uzAdj(2,2,2) = uz1(ixMax,iyMaxR,izMaxR)
         endif
      endif

      ! Sum up all and send to all
      call mpi_allreduce(uxAdj(1,1,1), uxAdj1(1,1,1), 8, real_type, mpi_sum,mpi_comm_world, ierr)
      call mpi_allreduce(uyAdj(1,1,1), uyAdj1(1,1,1), 8, real_type, mpi_sum,mpi_comm_world, ierr)
      call mpi_allreduce(uzAdj(1,1,1), uzAdj1(1,1,1), 8, real_type, mpi_sum,mpi_comm_world, ierr)

      uxAdj = uxAdj1
      uyAdj = uyAdj1
      uzAdj = uzAdj1

  
      ! Compute u,v,w at given point
      call triLinIterpolate(xd,yd,zd,uxAdj,uxInt)   
      call triLinIterpolate(xd,yd,zd,uyAdj,uyInt)   
      call triLinIterpolate(xd,yd,zd,uzAdj,uzInt)   

      strkX_p(i,j) = strkX_p(i,j)+ dt*uxInt
      strkY_p(i,j) = strkY_p(i,j)+ dt*uyInt
      strkZ_p(i,j) = strkZ_p(i,j)+ dt*uzInt
      
      ! Write data into file
      if(nrank == 0 .and. mod(itime,strkVizFreq_p) == 0) then
         if(strkSwitch_p .ne. 1) then
             write(dummystring,'(a,i3.3,a,i0.3)') 'streaklines/path',i,'At',j
             open(113,file=trim(dummystring),status='unknown',position='append')
             write(113,1234) itime*dt,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
             close(113)
         endif
         if(strkSwitch_p .ne. 2) then
            write(666,1235) itime,'  ',i,'  ',j,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
         endif
      endif  
      
      enddo
enddo 

if(nrank == 0) then
! Write restart positions
      if(mod(itime,2*isave) == 0) then
        open(10,file='streaklines/RestartPositions.dat',status='replace')
        write(10,*) strkPosns_p, '# Total particle positions'
        write(10,*) noOfReleases_p, '# Total particle releases'
        write(10,*) 'Position ReleaseTime X Y Z'
        do j = 1,noOfReleases_p
            do i = 1,strkPosns_p
               write(10,'(i0.3,a,i0.3,a,es15.8,a,es15.8,a,es15.8)')  i,'  ',j,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
            enddo
         enddo
         close(10)
      elseif(mod(itime,isave) == 0) then
        open(10,file='streaklines/RestartPositionsA.dat',status='replace')
        write(10,*) strkPosns_p, '# Total particle positions'
        write(10,*) noOfReleases_p, '# Total particle releases'
        write(10,*) 'Position ReleaseTime X Y Z'
        do j = 1,noOfReleases_p
            do i = 1,strkPosns_p
               write(10,'(i0.3,a,i0.3,a,es15.8,a,es15.8,a,es15.8)')  i,'  ',j,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
            enddo
         enddo
         close(10)
      endif
endif

!#########################
 end subroutine streakline
!#########################

!#########################
 subroutine visuNHalt(rho1, ux1, uy1, uz1, pp3, phi1, ep1)
!#########################
use decomp_2d, only : decomp_2d_abort
use decomp_2d, only : nrank, mytype, xsize, ph1
use visu, only  : write_snapshot, end_snapshot
use var, only : itime
use var, only : numscalar, nrhotime
use var, only : nzmsize, npress
use case, only : visu_case

use print_module

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

integer :: num

call write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime, num)
call visu_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1, num)
call end_snapshot(itime, num)

!call VISU_INSTA (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
!      ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
!      ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
!      ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)

if (nrank == 0) call printInfo_p('All positions have at least one particle which has &
                        & exited domain, so simulation done!')
call decomp_2d_abort(164,'All exited')

!#########################
 end subroutine visuNHalt
!#########################

!#########################
 subroutine triLinIterpolate(xd,yd,zd,uAdj,uInt)   
!#########################
 use decomp_2d, only: mytype
 implicit none
 
 real(mytype), intent(in) :: xd,yd,zd,uAdj(2,2,2)
 real(mytype), intent(out):: uInt
 real(mytype) :: c00,c01,c10,c11,c0,c1 
 
 c00 = uAdj(1,1,1)*(1d0-xd)+uAdj(2,1,1)*xd
 c01 = uAdj(1,1,2)*(1d0-xd)+uAdj(2,1,2)*xd
 c10 = uAdj(1,2,1)*(1d0-xd)+uAdj(2,2,1)*xd
 c11 = uAdj(1,2,2)*(1d0-xd)+uAdj(2,2,2)*xd

 c0 = c00*(1-yd)+c10*yd
 c1 = c01*(1-yd)+c11*yd

 uInt = c0*(1-zd)+c1*zd

!#########################
 end subroutine triLinIterpolate
!#########################

!@@10A End

!********************************
  end module tools_p
!********************************
