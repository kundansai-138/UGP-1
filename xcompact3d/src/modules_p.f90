module param_converging
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_con,del_con,del_cf_con,S_con,con_angle_degree,len_con
	
end module param_converging

module param_diverging
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_div,del_div,del_cf_div,S_div,div_angle_degree,len_div
	
end module param_diverging

module param_swirling
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_p,del_p,del_cf,S_p, len_nz	! Tolstoy Debbarma jet with nozzle
	
end module param_swirling

module param_swirling_nonozzle
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_nonz,del_nonz,del_cf_nonz,S_nonz,rao_nonz
	
end module param_swirling_nonozzle

module param_swirling_PM
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_pm, del_pm, S_pm
	
end module param_swirling_PM	
	
  module param_p
!********************************

! Module to set common parmeters I will be using
  use decomp_2d, only : mytype
  implicit none
  
  character, parameter :: foldName_p*4='out/'    ! Output Folder name
  
  real(mytype), parameter :: blowUpVel_p=1d1                       ! Velocity above which code is assumed to have blown
  
  integer, save :: i_retain_checkpoint, bet_start, bet_end, i_del_stat
  
  real(mytype), save :: convT
  
  integer, save :: convModulo_p  !@@ time iteration difference between two adjacent convergence velocities
  
  real(mytype), save :: convVelDiff_p	!Velocity difference between two successive visualizations, freq of viz
  
  ! Hysteresis
  integer, save :: ihys_td, istreakline, iSFD, iSFD_first_run
  real(mytype), save :: STarg_td, TTarg_td
  
  ! Sub-domain
  integer, save :: isubdom, subdom_together
  integer, save :: nx_first, nx_last, ny_first, ny_last, nz_first, nz_last
  
  integer, save :: iplanez0
  integer, save :: nx_first_z0, nx_last_z0, ny_first_z0, ny_last_z0, nz_first_z0, nz_last_z0
  
  integer, save :: iplaney0
  integer, save :: nx_first_y0, nx_last_y0, ny_first_y0, ny_last_y0, nz_first_y0, nz_last_y0
  
  !@@10A Streaklines
  integer,save :: strkSwitch_p,skipTime_p,strkPosns_p,noOfReleases_p,strkVizFreq_p
  real(mytype) :: bufferTime_p
  
  ! SFD parameters
  real(mytype), save :: wc_p, xi_p ! Filter cutoff freq and control parameter
  
  end module param_p
  

  module var_p 
!********************************

! Module to initialize common variables I will be using

 use decomp_2d, only : mytype
 implicit none

 character,save :: dummyString*200                         ! String used for multiple purposes

! Probes
  integer,save :: totPb_p						  	  !@@7A total no of probes	
  integer, allocatable, save, dimension(:) :: pbXPos_p,pbYPos_p,pbZPos_p  !@@7A Probe positions
  real(mytype),save,allocatable,dimension(:,:,:) :: pbU_p ! Probe save till isave 

!@@10A Streaklines
  real(mytype),allocatable,save, dimension(:,:) :: strkX_p,strkY_p,strkZ_p
  integer, allocatable,save, dimension(:) :: strkExit
  
  end module var_p
  
!#########################
  subroutine readProbes_p
!#########################

! Subroutine to read and initialize probe location

use decomp_2d, only: mytype
use var_p
use variables, only: nxm,nym,nzm
use param, only: xlx,yly,zlz,irestart,icheckpoint

implicit none

real(mytype) :: xPb,yPb,zPb	! Probe positions
integer :: i

!Probe file: probes.prm read
open(891,file = 'probes_points.prm',status='old')
read(891,*) dummyString
read(891,*) dummyString
read(891,*) dummyString

! No of probes
read(891,*) totPb_p
read(891,*) dummyString

! Allocate that many probes
allocate(pbXPos_p(totPb_p))
allocate(pbYPos_p(totPb_p))
allocate(pbZPos_p(totPb_p))
allocate(pbU_p(totPb_p,3,icheckpoint))

! Create files if starting afresh
if(irestart.eq.0) then
do i=1,totPb_p
   write(dummyString,'(a,i3.3)') 'probes/Pb',i
   open(892,file=trim(dummyString),status='replace')
   close(892)
enddo
endif

! Read positions & calculate indices
do i=1,totPb_p
   read(891,*) xPb,yPb,zPb
! Probe Global indices
   pbxPos_p(i) = nint(xPb/xlx*nxm)+1
   pbyPos_p(i) = nint(yPb/yly*nym)+1+nym/2
   pbzPos_p(i) = nint(zPb/zlz*nzm)+1+nzm/2
enddo

close(891)

!#########################
  end subroutine readProbes_p
!#########################

!#########################
subroutine probes_p(ux,uy,uz,t)
!#########################

! Subroutine to write velocity at given probe position
use var_p
use decomp_2d, only: mytype, xsize,xstart,xend
use param, only: icheckpoint
implicit none

real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
integer, intent(in) :: t
integer :: i,j,k,pbNo

! Expensive way of doing

do pbNo = 1,totPb_p

  if(pbyPos_p(pbNo) .ge. xstart(2) .and. pbyPos_p(pbNo) .le. xend(2) ) then
  if(pbzPos_p(pbNo) .ge. xstart(3) .and. pbzPos_p(pbNo) .le. xend(3) ) then

    i = pbxPos_p(pbNo)
    j = pbyPos_p(pbNo)-xstart(2)+1
    k = pbzPos_p(pbNo)-xstart(3)+1
    pbU_p(pbNo,1,mod(t-1,icheckpoint)+1) = ux(i,j,k)
    pbU_p(pbNo,2,mod(t-1,icheckpoint)+1) = uy(i,j,k)
    pbU_p(pbNo,3,mod(t-1,icheckpoint)+1) = uz(i,j,k)

    if(mod(t,icheckpoint) == 0) then
       write(dummyString,'(a,i3.3)') 'probes/Pb',pbNo
       open(892,file=trim(dummyString),status='unknown',action='write', &
       form='formatted',position="append")
       do i = 1,icheckpoint
           write(892,*)  t-icheckpoint+i, pbU_p(pbNo,1,i), pbU_p(pbNo,2,i), pbU_p(pbNo,3,i)
       enddo
       close(892)
    endif

  endif
  endif

enddo

!#########################
end subroutine probes_p
!#########################  

module print_module
  use var_p, only: dummyString
  implicit none
contains
subroutine printInfo_p(str1)

! subroutine to print a string if passed, else print 
! dummyString

character(len=*), optional, intent(in) :: str1

if(present(str1)) write(dummyString,*) str1

write(*,*) trim(dummyString)

end subroutine printInfo_p
end module print_module


!@@10A start
!#########################
  subroutine readStreakVars_p
!#########################
 use decomp_2d
 use var_p
 use param_p
 use param
 use print_module
 implicit none

 integer :: i,j,nr,noOfRings,totTimeIters
 real(mytype) :: rLow,rHigh,rStep,rCurr,flowThruTime,flowThruFctr,strkType

 
 1234 format(ES15.8,a,ES15.8,a,ES15.8,a,ES15.8)
 1235 format(I0.6,a,I0.6,a,I0.6,a,ES15.8,a,ES15.8,a,ES15.8)
 ! Read prm file
 open(10,file='streaklines.prm',status='old')
  do i=1,3
    read(10,*)  dummystring
 enddo
 
 read(10,*) strkSwitch_p 
 if(strkSwitch_p == 0) then  ! Streaklines Off
    close(10)
    return
 endif
 
 if(nrank == 0) bufferTime_p = 0d0 ! How long code has run after all particles have reached exit
 if(nrank == 0) call printInfo_p('Streakline vizualization:')
 read(10,*) strkType
 read(10,*) strkVizFreq_p
 read(10,*) skipTime_p
 read(10,*) noOfReleases_p
 read(10,*) flowThruFctr
 if(strkType == 1 .and. nrank ==0 ) call printInfo_p('Radial positioning')

 if(strkType == 3) then
   call restartStreaks ! Restarting from prev simulation
   return
 endif
  
 ! Set no of times particles are released
 if(noOfReleases_p <= 0) then ! Automatic calculation based on flowThru Time
    flowThruTime = xlx/u1*flowThruFctr
    noOfReleases_p = nint(flowThruTime/dt/skipTime_p)
 endif
   

 ! Read method type and characteristics
 do i=1,4
    read(10,*)  dummystring
 enddo
 if(strkType == 1) then ! Radial allocation
    read(10,*) rLow
    read(10,*) noOfRings
    read(10,*) rHigh
    read(10,*) nr
    strkPosns_p = noOfRings*nr+1 ! +1 to include (0,0,0) 
    rStep = (rHigh-rLow)/(noOfRings-1)
 elseif(strkType == 2) then ! Positional allocation
    do i=1,7
       read(10,*)  dummystring
    enddo
    read(10,*) strkPosns_p
 endif

 ! Allocate position arrays
 allocate(strkX_p(strkPosns_p,noOfReleases_p))
 allocate(strkY_p(strkPosns_p,noOfReleases_p))
 allocate(strkZ_p(strkPosns_p,noOfReleases_p))
 allocate(strkExit(strkPosns_p))
 strkExit = 0 
 ! Initialize positions for first set to be released & copy it to rest
 if(strkType == 1) then ! Radial allocn
    strkX_p(:,1) = 0d0
    do i = 1,noOfRings
       if(noOfRings == 1) then
           rCurr = rLow
       else
           rCurr = rLow+(i-1)*rStep
       endif
       do j = 1,nr
          strkY_p((i-1)*nr+j,1) = rCurr*cos(2d0*pi*(j-1d0)/nr)
          strkZ_p((i-1)*nr+j,1) = rCurr*sin(2d0*pi*(j-1d0)/nr)
       enddo
    enddo
    strkY_p(strkPosns_p,1) = 0d0 ! Zero point pushed to last
    strkZ_p(strkPosns_p,1) = 0d0
 else                ! Positional allocn   
    read(10,*) dummystring
    do i=1,strkPosns_p
       read(10,*) strkX_p(i,1), strkY_p(i,1), strkZ_p(i,1)  
    enddo
 endif

 do i = 2,noOfReleases_p
    strkX_p(:,i) = strkX_p(:,1)
    strkY_p(:,i) = strkY_p(:,1)
    strkZ_p(:,i) = strkZ_p(:,1)
 enddo
 close(10)

! Print info and store initial data in file
  if(nrank == 0) then
     write(dummyString,'(a,i0.3)')'No of positions = ', strkPosns_p
     call printInfo_p
     write(dummyString,'(a,i0.3)')'No of releases = ', noOfReleases_p
     call printInfo_p
     write(dummyString,'(a,i0.3)')'Total = ',noOfReleases_p*strkPosns_p
     call printInfo_p
     write(dummyString,'(a,i0.3)')'with successive itime interval = ', skipTime_p
     call printInfo_p
     if(strkSwitch_p .ne. 1) then
     j=1
        do i = 1,strkPosns_p
           write(dummystring,'(a,i3.3,a,i0.3)') 'streaklines/path',i,'At',j
           open(10,file=trim(dummystring),status='replace')
           write(10,*) 'Time        X        Y        Z'
           write(10,1234) 0d0,'  ',strkX_p(i,j),'  ',strkY_p(i,j),'  ',strkZ_p(i,j)
           close(10)
         enddo
     endif
  endif

! all in 1data file
  if(nrank == 0) then
    if(strkSwitch_p .ne. 2) then
       write(dummystring,'(a,i4.4,a)')'streaklines/allIn1',0,'.dat'
       open(666,file=trim(dummystring),status='replace')
!       write(10,*) 'Time    PositionNo    ReleaseNo    X        Y        Z'
       do i = 1,strkPosns_p
          write(666,1235) 0,'  ',i,'  ',1,'  ',strkX_p(i,1),'  ',strkY_p(i,1),'  ',strkZ_p(i,1)
       enddo
     endif
  endif


!#########################
  end subroutine readStreakVars_p
!#########################

!#########################
  subroutine restartStreaks
!#########################
! Read streak number and positions based on 
 use var_p
 use param_p
 use param, only: ifirst, icheckpoint
 use decomp_2d, only: nrank
 use print_module
 
 implicit none
 integer :: i,j,i1,j1 
 open(10,file='streaklines/RestartPositions.dat')
 read(10,*) strkPosns_p
 read(10,*) noOfReleases_p
 read(10,*) dummystring
 allocate(strkX_p(strkPosns_p,noOfReleases_p))
 allocate(strkY_p(strkPosns_p,noOfReleases_p))
 allocate(strkZ_p(strkPosns_p,noOfReleases_p))
 allocate(strkExit(strkPosns_p))
 strkExit = 0
 do j = 1,noOfReleases_p
     do i = 1,strkPosns_p
        read(10,*)  i1,j1,strkX_p(i,j),strkY_p(i,j),strkZ_p(i,j)
     enddo
 enddo
 close(10)
 if(nrank ==0) then
   write(dummyString,*) 'Restarting from file = ', ifirst/icheckpoint
   call printInfo_p
 endif
 ! Replace previous runs half written allIn1 File: ie starting new
 if(nrank == 0) then
    if(strkSwitch_p .ne. 2) then
       write(dummystring,'(a,i4.4,a)')'streaklines/allIn1',ifirst/icheckpoint,'.dat'
       open(666,file=trim(dummystring),status='replace')
     endif
  endif

!#########################
  end subroutine restartStreaks
!#########################

!@@10A end