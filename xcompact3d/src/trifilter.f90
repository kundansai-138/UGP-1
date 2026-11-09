module trifilter

USE variables, only: nx,ny,nz
USE decomp_2d
USE param

implicit none

! Rohan's Filter Parameters for LES modelling
real(mytype),save :: f_alphaX,f_alphaY,f_alphaZ
integer,save :: f_oneSideState,f_stateX,f_stateY,f_stateZ
integer,save :: f_order
integer,save :: f_leftSurfState,f_rightSurfState,f_topSurfState, &
                     f_bottomSurfState,f_frontSurfState,f_backSurfState

integer,save :: bstiX,beniX,bstjX,benjX,bstkX,benkX
integer,save :: bstiY,beniY,bstjY,benjY,bstkY,benkY
integer,save :: bstiZ,beniZ,bstjZ,benjZ,bstkZ,benkZ

integer,save :: f_stiX,f_stjX,f_stkX,f_eniX,f_enjX,f_enkX
integer,save :: f_stiY,f_stjY,f_stkY,f_eniY,f_enjY,f_enkY
integer,save :: f_stiZ,f_stjZ,f_stkZ,f_eniZ,f_enjZ,f_enkZ

! Rohan's filter used for LES modelling
real(mytype),save,allocatable,dimension(:,:) :: f_rhs2,f_rhs4,f_rhs6,f_rhs8,f_rhs10 
real(mytype),allocatable,save, dimension(:) :: f_triAx,f_triAy,f_triAz
real(mytype),allocatable,save, dimension(:) :: f_triBx,f_triBy,f_triBz
real(mytype),allocatable,save, dimension(:) :: f_triCx,f_triCy,f_triCz
real(mytype),allocatable,save, dimension(:) :: f_triA_Bpx,f_tri_Bpx
real(mytype),allocatable,save, dimension(:) :: f_triA_Bpy,f_tri_Bpy
real(mytype),allocatable,save, dimension(:) :: f_triA_Bpz,f_tri_Bpz
!Stores one sided filter coefficients
real(mytype),save,allocatable,dimension(:,:,:) :: f_rhs_one !second indice indicates the point number away from the


contains              
 !debug: see if you can reduce cache misses in this file

!#########################
 subroutine allocVars_p


    if(f_order == 10) then
      allocate(f_rhs2(2,3),f_rhs4(3,3),f_rhs6(4,3),f_rhs8(5,3),f_rhs10(6,3))
    else if(f_order == 8) then
      allocate(f_rhs2(2,3),f_rhs4(3,3),f_rhs6(4,3),f_rhs8(5,3))
    else if(f_order == 6) then
      allocate(f_rhs2(2,3),f_rhs4(3,3),f_rhs6(4,3))
    else if(f_order == 4) then
      allocate(f_rhs2(2,3),f_rhs4(3,3))
    else if(f_order == 2) then
      allocate(f_rhs2(2,3))
    endif

    if(f_oneSideState == 1) then
      if(f_order == 10) then
        allocate(f_rhs_one(11,4,3))
      else if(f_order == 8) then
        allocate(f_rhs_one(11,3,3))
      else if(f_order == 6) then
        allocate(f_rhs_one(11,2,3))
      else if(f_order == 4) then
        allocate(f_rhs_one(11,1,3))
      endif
    endif

    allocate(f_triAx(nx),f_triBx(nx),f_triCx(nx))
    allocate(f_triAy(ny),f_triBy(ny),f_triCy(ny))
    allocate(f_triAz(nz),f_triBz(nz),f_triCz(nz))
    allocate(f_triA_Bpx(nx),f_tri_Bpx(nx))
    allocate(f_triA_Bpy(ny),f_tri_Bpy(ny))
    allocate(f_triA_Bpz(nz),f_tri_Bpz(nz))

!#########################
  end subroutine allocVars_p
!#########################
                                                  

subroutine filterRHS(rhs,fX,nx,ny,nz,sti,eni,stj,enj,stk,enk,isPeriodic,dir)
 !Sets the RHS for the thomas algorithm for matrix inversion of the filters
  integer,intent(IN) :: nx,ny,nz,dir,sti,eni,stj,enj,stk,enk
  real(mytype),dimension(nx,ny,nz) :: rhs,fX
  logical,intent(IN) :: isPeriodic
  integer :: i,j,k,l,fOneStX,fOneEnX
  integer :: ipl,iml,actualOrder
  real(mytype) :: actualCoeff
  
  rhs=fX !Needed because we do not want to change the boundary values.
	 !Thus all loops in this file use the skipboundary parameters
    
  if(isPeriodic) then
    do k=stk,enk
      do j=stj,enj
	do i=sti,eni
	  rhs(i,j,k) = 0.
	  do l=0,f_order/2
	    
	    call calcFilterRHSCoeff(f_order,l+1,dir,actualCoeff)
	    
	    if(dir == 1) then!X-direction
	      ipl = i+l
	      iml = i-l
	      call findCorrectedPeriodicIndices(ipl,nx) 
	      call findCorrectedPeriodicIndices(iml,nx)
	      rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(ipl,j,k)+fX(iml,j,k))
	    else if(dir == 2) then!Y-direction
	      ipl = j+l
	      iml = j-l
	      call findCorrectedPeriodicIndices(ipl,ny) 
	      call findCorrectedPeriodicIndices(iml,ny)
	      rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,ipl,k)+fX(i,iml,k))
	    else if(dir == 3) then!Z-direction
	      ipl = k+l
	      iml = k-l
	      call findCorrectedPeriodicIndices(ipl,nz) 
	      call findCorrectedPeriodicIndices(iml,nz)
	      rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,j,ipl)+fX(i,j,iml))
	    endif
	  
	  enddo
	enddo
      enddo
    enddo
  else
    do k=stk,enk
      do j=stj,enj
	do i=sti,eni
	  rhs(i,j,k) = 0.
	  if(f_oneSideState == 0) then !order reduces near the boundary
	    
	    call calcActualOrder(actualOrder,i,j,k,sti,eni,stj,enj,stk,enk,dir)
	    
	    do l=0,actualOrder/2
		
		call calcFilterRHSCoeff(actualOrder,l+1,dir,actualCoeff)
		
		if(dir == 1) then!X-direction
		  ipl = i+l
		  iml = i-l
		  rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(ipl,j,k)+fX(iml,j,k))
		  if((i==sti).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaX*fX(i-1,j,k)
		  if((i==eni).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaX*fX(i+1,j,k)
		else if(dir == 2) then!Y-direction
		  ipl = j+l
		  iml = j-l
		  rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,ipl,k)+fX(i,iml,k))
		  if((j==stj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaY*fX(i,j-1,k)
		  if((j==enj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaY*fX(i,j+1,k)
		else if(dir == 3) then!Z-direction
		  ipl = k+l
		  iml = k-l
		  rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,j,ipl)+fX(i,j,iml))
		  if((k==stk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaZ*fX(i,j,k-1)
		  if((k==enk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaZ*fX(i,j,k+1)
		endif
	      enddo  
	    else !One sided schemes used near the boundary

	      do l=0,f_order/2
	
		call calcFilterRHSCoeff(f_order,l+1,dir,actualCoeff)
		
		if(dir == 1) then!X-direction
		  ipl = i+l
		  iml = i-l
		  
		  if(l==0) call calcOneSidedFilterIndices(sti,eni,f_order,fOneStX,fOneEnX)
		  
		  if((i.lt.fOneStX).and.(l==0)) then
		    rhs(i,j,k) = f_rhs_one(1,i-1,dir)*fX(1,j,k) + f_rhs_one(2,i-1,dir)*fX(2,j,k) + f_rhs_one(3,i-1,dir)*fX(3,j,k) + &
				 f_rhs_one(4,i-1,dir)*fX(4,j,k) + f_rhs_one(5,i-1,dir)*fX(5,j,k) + f_rhs_one(6,i-1,dir)*fX(6,j,k) + &
				 f_rhs_one(7,i-1,dir)*fX(7,j,k) + f_rhs_one(8,i-1,dir)*fX(8,j,k) + f_rhs_one(9,i-1,dir)*fX(9,j,k) + &
				 f_rhs_one(10,i-1,dir)*fX(10,j,k) + f_rhs_one(11,i-1,dir)*fX(11,j,k)
		    if(i==sti) rhs(i,j,k) = rhs(i,j,k) - f_alphaX*fX(i-1,j,k) 
		  else if((i.gt.fOneEnX).and.(l==0)) then
		    rhs(i,j,k) = f_rhs_one(1,eni-i+1,dir)*fX(eni+1,j,k) + f_rhs_one(2,eni-i+1,dir)*fX(eni,j,k) + &
				 f_rhs_one(3,eni-i+1,dir)*fX(eni-1,j,k) + f_rhs_one(4,eni-i+1,dir)*fX(eni-2,j,k) + &
				 f_rhs_one(5,eni-i+1,dir)*fX(eni-3,j,k) + f_rhs_one(6,eni-i+1,dir)*fX(eni-4,j,k) + &
				 f_rhs_one(7,eni-i+1,dir)*fX(eni-5,j,k) + f_rhs_one(8,eni-i+1,dir)*fX(eni-6,j,k) + &
				 f_rhs_one(9,eni-i+1,dir)*fX(eni-7,j,k) + f_rhs_one(10,eni-i+1,dir)*fX(eni-8,j,k) + &
				 f_rhs_one(11,eni-i+1,dir)*fX(eni-9,j,k)
		    if(i==eni) rhs(i,j,k) = rhs(i,j,k) - f_alphaX*fX(i+1,j,k)
		  else if((i.ge.fOneStX).and.(i.le.fOneEnX)) then
		    rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(ipl,j,k)+fX(iml,j,k))
		    if((i==eni).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaX*fX(i+1,j,k)!needed for f_order == 2
		    if((i==sti).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaX*fX(i-1,j,k)!needed for f_order == 2
		  endif
		else if(dir == 2) then!Y-direction
		  ipl = j+l
		  iml = j-l
		  
		  if(l==0) call calcOneSidedFilterIndices(stj,enj,f_order,fOneStX,fOneEnX)
		  
		  if((j.lt.fOneStX).and.(l==0)) then
		    rhs(i,j,k) = f_rhs_one(1,j-1,dir)*fX(i,1,k) + f_rhs_one(2,j-1,dir)*fX(i,2,k) + f_rhs_one(3,j-1,dir)*fX(i,3,k) + &
				 f_rhs_one(4,j-1,dir)*fX(i,4,k) + f_rhs_one(5,j-1,dir)*fX(i,5,k) + f_rhs_one(6,j-1,dir)*fX(i,6,k) + &
				 f_rhs_one(7,j-1,dir)*fX(i,7,k) + f_rhs_one(8,j-1,dir)*fX(i,8,k) + f_rhs_one(9,j-1,dir)*fX(i,9,k) + &
				 f_rhs_one(10,j-1,dir)*fX(i,10,k) + f_rhs_one(11,j-1,dir)*fX(i,11,k)
		  if(j==stj) rhs(i,j,k) = rhs(i,j,k) - f_alphaY*fX(i,j-1,k) 
		  else if((j.gt.fOneEnX).and.(l==0)) then
		    rhs(i,j,k) = f_rhs_one(1,enj-j+1,dir)*fX(i,enj+1,k) + f_rhs_one(2,enj-j+1,dir)*fX(i,enj,k) + &
				 f_rhs_one(3,enj-j+1,dir)*fX(i,enj-1,k) + f_rhs_one(4,enj-j+1,dir)*fX(i,enj-2,k) + &
				 f_rhs_one(5,enj-j+1,dir)*fX(i,enj-3,k) + f_rhs_one(6,enj-j+1,dir)*fX(i,enj-4,k) + &
				 f_rhs_one(7,enj-j+1,dir)*fX(i,enj-5,k) + f_rhs_one(8,enj-j+1,dir)*fX(i,enj-6,k) + &
				 f_rhs_one(9,enj-j+1,dir)*fX(i,enj-7,k) + f_rhs_one(10,enj-j+1,dir)*fX(i,enj-8,k) + &
				 f_rhs_one(11,enj-j+1,dir)*fX(i,enj-9,k)
		    if(j==enj) rhs(i,j,k) = rhs(i,j,k) - f_alphaY*fX(i,j+1,k)
		  else if((j.ge.fOneStX).and.(j.le.fOneEnX)) then
		    rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,ipl,k)+fX(i,iml,k))
		    if((j==enj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaY*fX(i,j+1,k)!needed for f_order == 2
		    if((j==stj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaY*fX(i,j-1,k)!needed for f_order == 2 
		  endif
		else if(dir == 3) then!Z-direction
		  ipl = k+l
		  iml = k-l
		  
		  if(l==0) call calcOneSidedFilterIndices(stk,enk,f_order,fOneStX,fOneEnX)
		  
		  if((k.lt.fOneStX).and.(l==0)) then
		    rhs(i,j,k) = f_rhs_one(1,k-1,dir)*fX(i,j,1) + f_rhs_one(2,k-1,dir)*fX(i,j,2) + f_rhs_one(3,k-1,dir)*fX(i,j,3) + &
				 f_rhs_one(4,k-1,dir)*fX(i,j,4) + f_rhs_one(5,k-1,dir)*fX(i,j,5) + f_rhs_one(6,k-1,dir)*fX(i,j,6) + &
				 f_rhs_one(7,k-1,dir)*fX(i,j,7) + f_rhs_one(8,k-1,dir)*fX(i,j,8) + f_rhs_one(9,k-1,dir)*fX(i,j,9) + &
				 f_rhs_one(10,k-1,dir)*fX(i,j,10) + f_rhs_one(11,k-1,dir)*fX(i,j,11)
		    if(k==stk) rhs(i,j,k) = rhs(i,j,k) - f_alphaZ*fX(i,j,k-1) 
		  else if((k.gt.fOneEnX).and.(l==0)) then
		    rhs(i,j,k) = f_rhs_one(1,enk-k+1,dir)*fX(i,j,enk+1) + f_rhs_one(2,enk-k+1,dir)*fX(i,j,enk) + &
				 f_rhs_one(3,enk-k+1,dir)*fX(i,j,enk-1) + f_rhs_one(4,enk-k+1,dir)*fX(i,j,enk-2) + &
				 f_rhs_one(5,enk-k+1,dir)*fX(i,j,enk-3) + f_rhs_one(6,enk-k+1,dir)*fX(i,j,enk-4) + &
				 f_rhs_one(7,enk-k+1,dir)*fX(i,j,enk-5) + f_rhs_one(8,enk-k+1,dir)*fX(i,j,enk-6) + &
				 f_rhs_one(9,enk-k+1,dir)*fX(i,j,enk-7) + f_rhs_one(10,enk-k+1,dir)*fX(i,j,enk-8) + &
				 f_rhs_one(11,enk-k+1,dir)*fX(i,j,enk-9)
		    if(k==enk) rhs(i,j,k) = rhs(i,j,k) - f_alphaZ*fX(i,j,k+1)
		  else if((k.ge.fOneStX).and.(k.le.fOneEnX)) then
		    rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,j,ipl)+fX(i,j,iml))
		    if((k==enk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaZ*fX(i,j,k+1)!needed for f_order == 2
		    if((k==stk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - f_alphaZ*fX(i,j,k-1)!needed for f_order == 2 
		  endif
		endif
	      
	      enddo
	    
	    endif
	  
	enddo
      enddo
    enddo
  endif
  
 end subroutine filterRHS
 
 subroutine filterPerturbRHS(pert,nx,ny,nz,sti,eni,stj,enj,stk,enk,dir)
 !This is needed for periodic domains where the inversion is obtained using the Sherman Morrison Formula
  integer,intent(IN) :: nx,ny,nz,dir,sti,eni,stj,enj,stk,enk
  real(mytype),dimension(nx,ny,nz),intent(INOUT) :: pert
  integer :: i,j,k
  
  if(dir == 1) then
    do k=stk,enk 
      do j=stj,enj
	pert(1,j,k) = -1.
	do i=2,nx-1
	  pert(i,j,k) = 0.
	enddo
	pert(nx,j,k) = f_alphaX
      enddo
    enddo  
  else if(dir == 2) then
    do k=stk,enk 
      do i=sti,eni
	pert(i,1,k) = -1.
	do j=2,ny-1
	  pert(i,j,k) = 0.
	enddo
	pert(i,ny,k) = f_alphaY
      enddo
    enddo
  else if(dir == 3) then
    do j=stj,enj 
      do i=sti,eni
	pert(i,j,1) = -1.
	do k=2,nz-1
	  pert(i,j,k) = 0.
	enddo
	pert(i,j,nz) = f_alphaZ
      enddo
    enddo
  endif
 end subroutine filterPerturbRHS
 
 subroutine calcOneSidedFilterIndices(stX,enX,actualOrder,f_stX,f_enX)
  !stX and enX are boundary corrected parameters i.e. bstiX,beniX,bstjY,benjY etc
  integer,intent(IN) :: stX,enX,actualOrder
  integer,intent(OUT) :: f_stX,f_enX
  
  if(actualOrder == 2) then
    f_stX = stX
    f_enX = enX
  else if(actualOrder == 4) then
    f_stX = stX+1
    f_enX = enX-1
  else if(actualOrder == 6) then
    f_stX = stX+2
    f_enX = enX-2
  else if(actualOrder == 8) then
    f_stX = stX+3
    f_enX = enX-3
  else if(actualOrder == 10) then
    f_stX = stX+4
    f_enX = enX-4
  endif
  
 end subroutine calcOneSidedFilterIndices
 
  subroutine calcFilterRHSCoeff(actualOrder,indice,dir,actualCoeff)
  integer,intent(IN) :: actualOrder,indice,dir
  real(mytype),intent(OUT) :: actualCoeff

  if(actualOrder == 2) then
    actualCoeff = f_rhs2(indice,dir)
  else if(actualOrder == 4) then
    actualCoeff = f_rhs4(indice,dir)
  else if(actualOrder == 6) then
    actualCoeff = f_rhs6(indice,dir)
  else if(actualOrder == 8) then
    actualCoeff = f_rhs8(indice,dir)
  else if(actualOrder == 10) then
    actualCoeff = f_rhs10(indice,dir)  
  endif

 end subroutine calcFilterRHSCoeff
  
 subroutine calcActualOrder(actualOrder,i,j,k,sti,eni,stj,enj,stk,enk,dir)
  integer,intent(IN) :: dir,sti,eni,stj,enj,stk,enk,i,j,k 
  integer,intent(OUT) :: actualOrder
  
  if (f_order == 2) then
    if((dir==1).and.(i==sti)) then
      actualOrder = 2
    else if((dir==1).and.(i==eni)) then
      actualOrder = 2
    else if((dir==2).and.(j==stj)) then
      actualOrder = 2
    else if((dir==2).and.(j==enj)) then
      actualOrder = 2  
    else if((dir==3).and.(k==stk)) then
      actualOrder = 2 
    else if((dir==3).and.(k==enk)) then
      actualOrder = 2
    else
      actualOrder=f_order
    endif
  else if(f_order == 4) then
    if((dir==1).and.(i==sti)) then
      actualOrder = 2
    else if((dir==1).and.(i==eni)) then
      actualOrder = 2
    else if((dir==2).and.(j==stj)) then
      actualOrder = 2
    else if((dir==2).and.(j==enj)) then
      actualOrder = 2  
    else if((dir==3).and.(k==stk)) then
      actualOrder = 2 
    else if((dir==3).and.(k==enk)) then
      actualOrder = 2
    else
      actualOrder=f_order
    endif
  else if(f_order == 6) then
    if((dir==1).and.(i==sti)) then
      actualOrder = 2
    else if((dir==1).and.(i==sti+1)) then
      actualOrder = 4
    else if((dir==1).and.(i==eni)) then
      actualOrder = 2
    else if((dir==1).and.(i==eni-1)) then
      actualOrder = 4
    else if((dir==2).and.(j==stj)) then
      actualOrder = 2
    else if((dir==2).and.(j==stj+1)) then
      actualOrder = 4
    else if((dir==2).and.(j==enj)) then
      actualOrder = 2  
    else if((dir==2).and.(j==enj-1)) then
      actualOrder = 4
    else if((dir==3).and.(k==stk)) then
      actualOrder = 2 
    else if((dir==3).and.(k==stk+1)) then
      actualOrder = 4
    else if((dir==3).and.(k==enk)) then
      actualOrder = 2
    else if((dir==3).and.(k==enk-1)) then
      actualOrder = 4
    else
      actualOrder=f_order
    endif
  else if(f_order == 8) then
    if((dir==1).and.(i==sti)) then
      actualOrder = 2
    else if((dir==1).and.(i==sti+1)) then
      actualOrder = 4
    else if((dir==1).and.(i==sti+2)) then
      actualOrder = 6
    else if((dir==1).and.(i==eni)) then
      actualOrder = 2
    else if((dir==1).and.(i==eni-1)) then
      actualOrder = 4
    else if((dir==1).and.(i==eni-2)) then
      actualOrder = 6
    else if((dir==2).and.(j==stj)) then
      actualOrder = 2
    else if((dir==2).and.(j==stj+1)) then
      actualOrder = 4
    else if((dir==2).and.(j==stj+2)) then
      actualOrder = 6
    else if((dir==2).and.(j==enj)) then
      actualOrder = 2  
    else if((dir==2).and.(j==enj-1)) then
      actualOrder = 4
    else if((dir==2).and.(j==enj-2)) then
      actualOrder = 6
    else if((dir==3).and.(k==stk)) then
      actualOrder = 2 
    else if((dir==3).and.(k==stk+1)) then
      actualOrder = 4
    else if((dir==3).and.(k==stk+2)) then
      actualOrder = 6
    else if((dir==3).and.(k==enk)) then
      actualOrder = 2
    else if((dir==3).and.(k==enk-1)) then
      actualOrder = 4
    else if((dir==3).and.(k==enk-2)) then
      actualOrder = 6
    else
      actualOrder=f_order
    endif
  else if(f_order == 10) then
    if((dir==1).and.(i==sti)) then
      actualOrder = 2
    else if((dir==1).and.(i==sti+1)) then
      actualOrder = 4
    else if((dir==1).and.(i==sti+2)) then
      actualOrder = 6
    else if((dir==1).and.(i==sti+3)) then
      actualOrder = 8
    else if((dir==1).and.(i==eni)) then
      actualOrder = 2
    else if((dir==1).and.(i==eni-1)) then
      actualOrder = 4
    else if((dir==1).and.(i==eni-2)) then
      actualOrder = 6
    else if((dir==1).and.(i==eni-3)) then
      actualOrder = 8
    else if((dir==2).and.(j==stj)) then
      actualOrder = 2
    else if((dir==2).and.(j==stj+1)) then
      actualOrder = 4
    else if((dir==2).and.(j==stj+2)) then
      actualOrder = 6
    else if((dir==2).and.(j==stj+3)) then
      actualOrder = 8
    else if((dir==2).and.(j==enj)) then
      actualOrder = 2  
    else if((dir==2).and.(j==enj-1)) then
      actualOrder = 4
    else if((dir==2).and.(j==enj-2)) then
      actualOrder = 6
    else if((dir==2).and.(j==enj-3)) then
      actualOrder = 8
    else if((dir==3).and.(k==stk)) then
      actualOrder = 2 
    else if((dir==3).and.(k==stk+1)) then
      actualOrder = 4
    else if((dir==3).and.(k==stk+2)) then
      actualOrder = 6
    else if((dir==3).and.(k==stk+3)) then
      actualOrder = 8
    else if((dir==3).and.(k==enk)) then
      actualOrder = 2
    else if((dir==3).and.(k==enk-1)) then
      actualOrder = 4
    else if((dir==3).and.(k==enk-2)) then
      actualOrder = 6
    else if((dir==3).and.(k==enk-3)) then
      actualOrder = 8
    else
      actualOrder=f_order
    endif
  endif
    
 end subroutine calcActualOrder
 

 subroutine setFilterIndices()
  !Sets indices for the filters. Used when any of the surfaces have to be
  !filtered as well.
! call skipBoundaryIndices()	
  call calcActualFilterIndices(f_stiX,f_stjX,f_stkX,f_eniX,f_enjX,f_enkX,1)
  call calcActualFilterIndices(f_stiY,f_stjY,f_stkY,f_eniY,f_enjY,f_enkY,2)
  call calcActualFilterIndices(f_stiZ,f_stjZ,f_stkZ,f_eniZ,f_enjZ,f_enkZ,3)
 end subroutine setFilterIndices
 
  subroutine calcActualFilterIndices(f_sti,f_stj,f_stk,f_eni,f_enj,f_enk,dir)
  
  integer,intent(OUT) :: f_sti,f_stj,f_stk,f_eni,f_enj,f_enk
  integer,intent(IN) :: dir
  
  if(dir == 1) then
    f_sti = bstiX
    f_stj = bstjX
    f_stk = bstkX
    f_eni = beniX
    f_enj = benjX
    f_enk = benkX
    
    if(.not. ncly) then
      if((f_topSurfState == 1).and.(xend(2) == ny_global)) f_enj = xsize(2)
      if((f_bottomSurfState == 1).and.(xstart(2) == 1)) f_stj = 1
    endif
    if(.not. nclz) then
      if((f_frontSurfState == 1).and.(xend(3) == nz_global)) f_enk = xsize(3)
      if((f_backSurfState == 1).and.(xstart(3) == 1)) f_stk = 1
    endif
  else if(dir == 2) then
    f_sti = bstiY
    f_stj = bstjY
    f_stk = bstkY
    f_eni = beniY
    f_enj = benjY
    f_enk = benkY
    
    if(.not. nclx) then
      if((f_rightSurfState == 1).and.(yend(1) == nx_global)) f_eni = ysize(1)
      if((f_leftSurfState == 1).and.(ystart(1) == 1)) f_sti = 1
    endif
    if(.not. nclz) then
      if((f_frontSurfState == 1).and.(yend(3) == nz_global)) f_enk = ysize(3)
      if((f_backSurfState == 1).and.(ystart(3) == 1)) f_stk = 1
    endif
    
  else if(dir == 3) then
    f_sti = bstiZ
    f_stj = bstjZ
    f_stk = bstkZ
    f_eni = beniZ
    f_enj = benjZ
    f_enk = benkZ
      
    if(.not. nclx) then
      if((f_rightSurfState == 1).and.(zend(1) == nx_global)) f_eni = zsize(1)
      if((f_leftSurfState == 1).and.(zstart(1) == 1)) f_sti = 1
    endif
    if(.not. ncly) then
      if((f_topSurfState == 1).and.(zend(2) == ny_global)) f_enj = zsize(2)
      if((f_bottomSurfState == 1).and.(zstart(2) == 1)) f_stj = 1
    endif
  endif
  
 end subroutine calcActualFilterIndices
   
subroutine triInvFilX(tx,ux,rx,sx,triC,triA_Bp,tri_Bp,nx,ny,nz,isPeriodic)
  
  integer :: nx,ny,nz,i,j,k 
  real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
  real(mytype), dimension(ny,nz):: sx
  real(mytype), dimension(nx):: triC,triA_Bp,tri_Bp 
  logical,intent(IN) :: isPeriodic
    
  if(isPeriodic) then
    call filterRHS(tx,ux,nx,ny,nz,f_stiX,f_eniX,f_stjX,f_enjX,f_stkX,f_enkX,isPeriodic,1)
    call filterPerturbRHS(rx,nx,ny,nz,f_stiX,f_eniX,f_stjX,f_enjX,f_stkX,f_enkX,1)
    do k=f_stkX,f_enkX
      do j=f_stjX,f_enjX 
	
	  do i=2, nx
	    tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*triA_Bp(i) 
	    rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*triA_Bp(i) 
	  enddo
	  tx(nx,j,k)=tx(nx,j,k)*tri_Bp(nx) 
	  rx(nx,j,k)=rx(nx,j,k)*tri_Bp(nx) 
	  do i=nx-1,1,-1
	    tx(i,j,k)=(tx(i,j,k)-triC(i)*tx(i+1,j,k))*tri_Bp(i) 
	    rx(i,j,k)=(rx(i,j,k)-triC(i)*rx(i+1,j,k))*tri_Bp(i) 
	  enddo
	  sx(j,k)=(tx(1,j,k)-f_alphaX*tx(nx,j,k))&
	      /(1.+rx(1,j,k)-f_alphaX*rx(nx,j,k)) 
	  do i=1,nx 
	    tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
	  enddo
	  
      enddo
    enddo
! do k=1,nz 
!     do j=1,ny 
! 	tx(1,j,k)=f_rhs(1)*ux(1,j,k) + 0.5*f_rhs(2)*(ux(2,j,k)+ux(nx,j,k)) &
! 		  + 0.5*f_rhs(3)*(ux(3,j,k)+ux(nx-1,j,k))
! 	rx(1,j,k)=-1. 
! 	tx(2,j,k)=f_rhs(1)*ux(2,j,k) + 0.5*f_rhs(2)*(ux(3,j,k)+ux(1,j,k)) &
! 		  + 0.5*f_rhs(3)*(ux(4,j,k)+ux(nx,j,k))
! 	rx(2,j,k)=0. 
! 	do i=3,nx-2
! 	  tx(i,j,k)=f_rhs(1)*ux(i,j,k) + 0.5*f_rhs(2)*(ux(i+1,j,k)+ux(i-1,j,k)) &
! 		  + 0.5*f_rhs(3)*(ux(i+2,j,k)+ux(i-2,j,k))
! 	  rx(i,j,k)=0. 
! 	enddo
! 	tx(nx-1,j,k)=f_rhs(1)*ux(nx-1,j,k) + 0.5*f_rhs(2)*(ux(nx,j,k)+ux(nx-2,j,k)) &
! 		  + 0.5*f_rhs(3)*(ux(1,j,k)+ux(nx-3,j,k))
! 	rx(nx-1,j,k)=0. 
! 	tx(nx,j,k)=f_rhs(1)*ux(nx,j,k) + 0.5*f_rhs(2)*(ux(1,j,k)+ux(nx-1,j,k)) &
! 		  + 0.5*f_rhs(3)*(ux(2,j,k)+ux(nx-2,j,k)) 
! 	rx(nx,j,k)=f_alpha           
! 	do i=2, nx
! 	  tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*triA_Bp(i) 
! 	  rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*triA_Bp(i) 
! 	enddo
! 	tx(nx,j,k)=tx(nx,j,k)*tri_Bp(nx) 
! 	rx(nx,j,k)=rx(nx,j,k)*tri_Bp(nx) 
! 	do i=nx-1,1,-1
! 	  tx(i,j,k)=(tx(i,j,k)-triC(i)*tx(i+1,j,k))*tri_Bp(i) 
! 	  rx(i,j,k)=(rx(i,j,k)-triC(i)*rx(i+1,j,k))*tri_Bp(i) 
! 	enddo
! 	sx(j,k)=(tx(1,j,k)-f_alpha*tx(nx,j,k))&
! 	    /(1.+rx(1,j,k)-f_alpha*rx(nx,j,k)) 
! 	do i=1,nx 
! 	  tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
! 	enddo
!     enddo
!     enddo
   else
    call filterRHS(tx,ux,nx,ny,nz,f_stiX,f_eniX,f_stjX,f_enjX,f_stkX,f_enkX,isPeriodic,1)
    do k=f_stkX,f_enkX
      do j=f_stjX,f_enjX
	do i=3,nx-1 
	  tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*triA_Bp(i) 
	enddo
	tx(nx-1,j,k)=tx(nx-1,j,k)*tri_Bp(nx-1) 
	do i=nx-2,2,-1
	  tx(i,j,k)=(tx(i,j,k)-triC(i)*tx(i+1,j,k))*tri_Bp(i) 
	enddo
      enddo 
    enddo 

    
!     tx = ux
!     do k=f_stk,f_enk
!       do j=f_stj,f_enj
! 	  tx(1,j,k)=ux(1,j,k) 
! ! 	  tx(2,j,k)= f_rhs_red1(1)*ux(2,j,k) + 0.5*f_rhs_red1(2)*(ux(1,j,k)+ux(3,j,k)) &
! ! 			-f_alpha*ux(1,j,k)
! 	  tx(2,j,k)= f_one(1)*ux(1,j,k) + f_one(2)*ux(2,j,k) + f_one(3)*ux(3,j,k) + &
! 		     f_one(4)*ux(4,j,k) + f_one(5)*ux(5,j,k) - f_alpha*ux(1,j,k) 
! 	  do i=3,nx-2
! 	    tx(i,j,k)=f_rhs(1)*ux(i,j,k) + 0.5*f_rhs(2)*(ux(i+1,j,k)+ux(i-1,j,k)) &
! 		  + 0.5*f_rhs(3)*(ux(i+2,j,k)+ux(i-2,j,k))
! 	  enddo
! ! 	  tx(nx-1,j,k)= f_rhs_red1(1)*ux(nx-1,j,k) + 0.5*f_rhs_red1(2)*(ux(nx-2,j,k)+ux(nx,j,k)) &
! ! 			-f_alpha*ux(nx,j,k)
! 	  tx(nx-1,j,k)= f_one(1)*ux(nx,j,k) + f_one(2)*ux(nx-1,j,k) + f_one(3)*ux(nx-2,j,k) + &
! 		     f_one(4)*ux(nx-3,j,k) + f_one(5)*ux(nx-4,j,k) - f_alpha*ux(nx,j,k)
! 	  tx(nx,j,k)=ux(nx,j,k)
! 
! 	  do i=3,nx-1 
! 	    tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*triA_Bp(i) 
! 	  enddo
! 	  tx(nx-1,j,k)=tx(nx-1,j,k)*tri_Bp(nx-1) 
! 	  do i=nx-2,2,-1
! 	    tx(i,j,k)=(tx(i,j,k)-triC(i)*tx(i+1,j,k))*tri_Bp(i) 
! 	  enddo
!       enddo 
!     enddo 
   
   endif

  end subroutine triInvFilX
  
 subroutine triInvFilY(ty,uy,ry,sy,triC,triA_Bp,tri_Bp,nx,ny,nz,isPeriodic) 
  
  integer :: nx,ny,nz,i,j,k
  real(mytype), dimension(nx,ny,nz) :: ty,uy 
  real(mytype), dimension(nx,ny,nz) :: ry
  real(mytype), dimension(nx,nz)  :: sy
  real(mytype), dimension(ny) :: triC,triA_Bp,tri_Bp 
  logical,intent(IN) :: isPeriodic
    
  if(isPeriodic) then
    call filterRHS(ty,uy,nx,ny,nz,f_stiY,f_eniY,f_stjY,f_enjY,f_stkY,f_enkY,isPeriodic,2)
    call filterPerturbRHS(ry,nx,ny,nz,f_stiY,f_eniY,f_stjY,f_enjY,f_stkY,f_enkY,2)
   do k=f_stkY,f_enkY
    do j=2,ny  
      do i=f_stiY,f_eniY 
	  ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*triA_Bp(j) 
	  ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*triA_Bp(j) 
      enddo
    enddo 
   enddo
   do k=f_stkY,f_enkY
    do i=f_stiY,f_eniY
	ty(i,ny,k)=ty(i,ny,k)*tri_Bp(ny) 
	ry(i,ny,k)=ry(i,ny,k)*tri_Bp(ny) 
    enddo
   enddo 
   do k=f_stkY,f_enkY
    do j=ny-1,1,-1  
      do i=f_stiY,f_eniY  
	  ty(i,j,k)=(ty(i,j,k)-triC(j)*ty(i,j+1,k))*tri_Bp(j) 
	  ry(i,j,k)=(ry(i,j,k)-triC(j)*ry(i,j+1,k))*tri_Bp(j) 
      enddo
    enddo 
   enddo
   do k=f_stkY,f_enkY
    do i=f_stiY,f_eniY  
	sy(i,k)=(ty(i,1,k)-f_alphaY*ty(i,ny,k))&
	    /(1.+ry(i,1,k)-f_alphaY*ry(i,ny,k)) 
    enddo
   enddo 
   do k=f_stkY,f_enkY
    do j=1,ny  
      do i=f_stiY,f_eniY  
	  ty(i,j,k)=ty(i,j,k)-sy(i,k)*ry(i,j,k) 
      enddo
    enddo
   enddo
 
  else
!     tmp = beniY
!     if(tmp == ysize(1)-1) tmp=tmp+1
    call filterRHS(ty,uy,nx,ny,nz,f_stiY,f_eniY,f_stjY,f_enjY,f_stkY,f_enkY,isPeriodic,2)
    do k=f_stkY,f_enkY
      do j=3,ny-1  
	do i=f_stiY,f_eniY  
	    ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*triA_Bp(j) 
	enddo 
      enddo 
    enddo 
    do k=f_stkY,f_enkY
      do i=f_stiY,f_eniY  
	  ty(i,ny-1,k)=ty(i,ny-1,k)*tri_Bp(ny-1) 
      enddo 
    enddo 
    do k=f_stkY,f_enkY
      do j=ny-2,2,-1  
	do i=f_stiY,f_eniY  
	    ty(i,j,k)=(ty(i,j,k)-triC(j)*ty(i,j+1,k))*tri_Bp(j) 
	enddo 
      enddo 
    enddo
  endif
 end subroutine triInvFilY
  
subroutine triInvFilZ(tz,uz,rz,sz,triC,triA_Bp,tri_Bp,nx,ny,nz,isPeriodic) 

  integer :: nx,ny,nz,i,j,k
  real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
  real(mytype), dimension(nx,ny) :: sz 
  real(mytype), dimension(nz) :: triC,triA_Bp,tri_Bp
  logical,intent(IN) :: isPeriodic
    
  if(isPeriodic) then
    call filterRHS(tz,uz,nx,ny,nz,f_stiZ,f_eniZ,f_stjZ,f_enjZ,f_stkZ,f_enkZ,isPeriodic,3)
    call filterPerturbRHS(rz,nx,ny,nz,f_stiZ,f_eniZ,f_stjZ,f_enjZ,f_stkZ,f_enkZ,3)
    do k=2,nz
      do j=f_stjZ,f_enjZ
	do i=f_stiZ,f_eniZ
	    tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*triA_Bp(k)
	    rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*triA_Bp(k)
	enddo
      enddo
    enddo
    do j=f_stjZ,f_enjZ
      do i=f_stiZ,f_eniZ
	  tz(i,j,nz)=tz(i,j,nz)*tri_Bp(nz)
	  rz(i,j,nz)=rz(i,j,nz)*tri_Bp(nz)
      enddo
    enddo
    do k=nz-1,1,-1
      do j=f_stjZ,f_enjZ
	do i=f_stiZ,f_eniZ
	    tz(i,j,k)=(tz(i,j,k)-triC(k)*tz(i,j,k+1))*tri_Bp(k)
	    rz(i,j,k)=(rz(i,j,k)-triC(k)*rz(i,j,k+1))*tri_Bp(k)
	enddo
      enddo
    enddo
    do j=f_stjZ,f_enjZ
      do i=f_stiZ,f_eniZ
	  sz(i,j)=(tz(i,j,1)-f_alphaZ*tz(i,j,nz))/&
	      (1.+rz(i,j,1)-f_alphaZ*rz(i,j,nz))
      enddo
    enddo
    do k=1,nz
      do j=f_stjZ,f_enjZ
	do i=f_stiZ,f_eniZ
	    tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
	enddo
      enddo
    enddo
    
    else
      call filterRHS(tz,uz,nx,ny,nz,f_stiZ,f_eniZ,f_stjZ,f_enjZ,f_stkZ,f_enkZ,isPeriodic,3)
      do k=3,nz-1
	do j=f_stjZ,f_enjZ
	  do i=f_stiZ,f_eniZ
	    tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*triA_Bp(k)
	  enddo
	enddo
      enddo
      do j=f_stjZ,f_enjZ
	do i=f_stiZ,f_eniZ
	  tz(i,j,nz-1)=tz(i,j,nz-1)*tri_Bp(nz-1)
	enddo
      enddo
      do k=nz-2,2,-1
	do j=f_stjZ,f_enjZ
	  do i=f_stiZ,f_eniZ
	    tz(i,j,k)=(tz(i,j,k)-triC(k)*tz(i,j,k+1))*tri_Bp(k)
	  enddo
	enddo
      enddo
    endif 
  end subroutine triInvFilZ
 
 
 subroutine computeTriFilteredVar(unfilteredVar)
 !This subroutine computes the filtered variable in an MPI communication effective way
 
     USE var, ONLY : ta2
     USE var, ONLY : ta3
     USE var, ONLY : td1,di1 
     USE var, ONLY : tb2,di2
     USE var, ONLY : tb3,di3
     USE var, ONLY : sx,sy,sz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(INOUT) :: unfilteredVar
      
      if((f_stateX==1).and.(f_stateY==1).and.(f_stateZ==1)) then
	
	call triInvFilX(td1,unfilteredVar,di1,sx,f_triCx,f_triA_Bpx,&
	      f_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	
	call transpose_x_to_y(td1,ta2)
	call triInvFilY(tb2,ta2,di2,sy,f_triCy,f_triA_Bpy,&
	      f_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)

	call transpose_y_to_z(tb2,ta3)
	call triInvFilZ(tb3,ta3,di3,sz,f_triCz,f_triA_Bpz,&
	      f_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
			
      else if((f_stateX==1).and.(f_stateY==1)) then
	
	call triInvFilX(td1,unfilteredVar,di1,sx,f_triCx,f_triA_Bpx,&
	      f_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	
	call transpose_x_to_y(td1,ta2)
	call triInvFilY(tb2,ta2,di2,sy,f_triCy,f_triA_Bpy,&
	      f_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)
	call transpose_y_to_x(tb2,unfilteredVar)      
      
      else if((f_stateX==1).and.(f_stateZ==1)) then
	
	call triInvFilX(td1,unfilteredVar,di1,sx,f_triCx,f_triA_Bpx,&
	      f_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	
	call transpose_x_to_y(td1,tb2)
	call transpose_y_to_z(tb2,ta3)
	call triInvFilZ(tb3,ta3,di3,sz,f_triCz,f_triA_Bpz,&
	      f_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
      
      else if((f_stateY==1).and.(f_stateZ==1)) then
	
	call transpose_x_to_y(unfilteredVar,ta2)
	call triInvFilY(tb2,ta2,di2,sy,f_triCy,f_triA_Bpy,&
	      f_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)

	call transpose_y_to_z(tb2,ta3)
	call triInvFilZ(tb3,ta3,di3,sz,f_triCz,f_triA_Bpz,&
	      f_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
      
      else if(f_stateX==1) then
      
      call triInvFilX(td1,unfilteredVar,di1,sx,f_triCx,f_triA_Bpx,&
	      f_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	unfilteredVar = td1
      
      else if(f_stateY==1) then
	
	call transpose_x_to_y(unfilteredVar,ta2)
	call triInvFilY(tb2,ta2,di2,sy,f_triCy,f_triA_Bpy,&
	      f_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)
	call transpose_y_to_x(tb2,unfilteredVar)
	
      else if(f_stateZ==1) then
	
	call transpose_x_to_y(unfilteredVar,tb2)
	call transpose_y_to_z(tb2,ta3)
	call triInvFilZ(tb3,ta3,di3,sz,f_triCz,f_triA_Bpz,&
	      f_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
      
      endif
                  
 end subroutine computeTriFilteredVar
 
  subroutine triFilterSchemes()
    
   integer :: i 
   real(mytype) :: alpha
   
   do i=1,3
    
    if(i==1) alpha = f_alphaX
    if(i==2) alpha = f_alphaY
    if(i==3) alpha = f_alphaZ
    
    !filter coefficients for filters which reduce order near the boundary
    if(f_order .ge. 2) then
      f_rhs2(1,i) = 0.5 + alpha
      f_rhs2(2,i) = 0.5 + alpha
    endif
    
    if(f_order .ge. 4) then
      f_rhs4(1,i) = 5./8. + 3.*alpha/4.
      f_rhs4(2,i) = 0.5 + alpha
      f_rhs4(3,i) = -1./8. + alpha/4.
    endif
    
    if(f_order .ge. 6) then
      f_rhs6(1,i) = 11./16. + 5.*alpha/8.
      f_rhs6(2,i) = 15./32. + 17.*alpha/16.
      f_rhs6(3,i) = -3./16. + 3.*alpha/8.
      f_rhs6(4,i) = 1./32. - alpha/16.
    endif
    
    if(f_order .ge. 8) then
      f_rhs8(1,i) = (93.+70.*alpha)/128.
      f_rhs8(2,i) = (7.+18.*alpha)/16.
      f_rhs8(3,i) = (-7.+14.*alpha)/32.
      f_rhs8(4,i) = 1./16. - alpha/8.
      f_rhs8(5,i) = -1./128. + alpha/64.
    endif
    
    if(f_order .ge. 10) then
      f_rhs10(1,i) = (193.+126.*alpha)/256.
      f_rhs10(2,i) = (105.+302.*alpha)/256.
      f_rhs10(3,i) = 15.*(-1.+2.*alpha)/64.
      f_rhs10(4,i) = 45.*(1.-2.*alpha)/512.
      f_rhs10(5,i) = 5.*(-1.+2.*alpha)/256.
      f_rhs10(6,i) = (1.-2.*alpha)/512.
    endif
    
    !filter coefficients for one-sided formulation
    if(f_oneSideState == 1) then
      
      if(f_order == 10) then
	!point-2 coefficients
	f_rhs_one(1,1,i) = (1.+1022.*alpha)/1024.
	f_rhs_one(2,1,i) = (507.+10.*alpha)/512.
	f_rhs_one(3,1,i) = (45.+934.*alpha)/1024.
	f_rhs_one(4,1,i) = 15.*(-1.+2.*alpha)/128.
	f_rhs_one(5,1,i) = 105.*(1.-2.*alpha)/512.
	f_rhs_one(6,1,i) = 63.*(-1.+2.*alpha)/256.
	f_rhs_one(7,1,i) = 105.*(1.-2.*alpha)/512.
	f_rhs_one(8,1,i) = 15.*(-1.+2.*alpha)/128.
	f_rhs_one(9,1,i) = 45.*(1.-2.*alpha)/1024.
	f_rhs_one(10,1,i) = 5.*(-1.+2.*alpha)/512.
	f_rhs_one(11,1,i) = (1.-2.*alpha)/1024.
	
	!point-3 coefficients
	f_rhs_one(1,2,i) = (-1.+2.*alpha)/1024.
	f_rhs_one(2,2,i) = (5.+502.*alpha)/512.
	f_rhs_one(3,2,i) = (979.+90.*alpha)/1024.
	f_rhs_one(4,2,i) = (15.+98.*alpha)/128.
	f_rhs_one(5,2,i) = 105.*(-1.+2.*alpha)/512.
	f_rhs_one(6,2,i) = 63.*(1.-2.*alpha)/256.
	f_rhs_one(7,2,i) = 105.*(-1.+2.*alpha)/512.
	f_rhs_one(8,2,i) = 15.*(1.-2.*alpha)/128.
	f_rhs_one(9,2,i) = 45.*(-1.+2.*alpha)/1024.
	f_rhs_one(10,2,i) = 5.*(1.-2.*alpha)/512.
	f_rhs_one(11,2,i) = (-1.+2.*alpha)/1024.
	
	!point-4 coefficients
	f_rhs_one(1,3,i) = (1.-2.*alpha)/1024.
	f_rhs_one(2,3,i) = 5.*(-1.+2.*alpha)/512.
	f_rhs_one(3,3,i) = (45.+934.*alpha)/1024.
	f_rhs_one(4,3,i) = (113.+30.*alpha)/128.
	f_rhs_one(5,3,i) = (105.+302.*alpha)/512.
	f_rhs_one(6,3,i) = 63.*(-1.+2.*alpha)/256.
	f_rhs_one(7,3,i) = 105.*(1.-2.*alpha)/512.
	f_rhs_one(8,3,i) = 15.*(-1.+2.*alpha)/128.
	f_rhs_one(9,3,i) = 45.*(1.-2.*alpha)/1024.
	f_rhs_one(10,3,i) = 5.*(-1.+2.*alpha)/512.
	f_rhs_one(11,3,i) = (1.-2.*alpha)/1024.
	
	!point-5 coefficients
	f_rhs_one(1,4,i) = (-1.+2.*alpha)/1024.
	f_rhs_one(2,4,i) = 5.*(1.-2.*alpha)/512.
	f_rhs_one(3,4,i) = 45.*(-1.+2.*alpha)/1024.
	f_rhs_one(4,4,i) = (15.+98.*alpha)/128.
	f_rhs_one(5,4,i) = (407.+210.*alpha)/512.
	f_rhs_one(6,4,i) = (63.+130.*alpha)/256.
	f_rhs_one(7,4,i) = 105.*(-1.+2.*alpha)/512.
	f_rhs_one(8,4,i) = 15.*(1.-2.*alpha)/128.
	f_rhs_one(9,4,i) = 45.*(-1.+2.*alpha)/1024.
	f_rhs_one(10,4,i) = 5.*(1.-2.*alpha)/512.
	f_rhs_one(11,4,i) = (-1.+2.*alpha)/1024.
      else if(f_order == 8) then
	!point-2 coefficients
	f_rhs_one(1,1,i) = 1./256. + (127.*alpha)/128.
	f_rhs_one(2,1,i) = 31./32. + (alpha)/16.
	f_rhs_one(3,1,i) = 7./64. + (25.*alpha)/32.
	f_rhs_one(4,1,i) = -7./32. + (7.*alpha)/16.
	f_rhs_one(5,1,i) = 7.*(5.-10*alpha)/128.
	f_rhs_one(6,1,i) = -7./32. + (7.*alpha)/16.
	f_rhs_one(7,1,i) = 7./64. - (7.*alpha)/32.
	f_rhs_one(8,1,i) = -1./32. + (alpha)/16.
	f_rhs_one(9,1,i) = 1./256. - (alpha)/128.
	f_rhs_one(10,1,i) = 0.
	f_rhs_one(11,1,i) = 0.
	
	!point-3 coefficients
	f_rhs_one(1,2,i) = -1./256. + (alpha)/128.
	f_rhs_one(2,2,i) = 1./32. + (15.*alpha)/16.
	f_rhs_one(3,2,i) = 57./64. + (7.*alpha)/32.
	f_rhs_one(4,2,i) = 7./32. + (9.*alpha)/16.
	f_rhs_one(5,2,i) = 7.*(-5.+10.*alpha)/128.
	f_rhs_one(6,2,i) = 7./32. - (7.*alpha)/16.
	f_rhs_one(7,2,i) = -7./64. + (7.*alpha)/32.
	f_rhs_one(8,2,i) = 1./32. - (alpha)/16.
	f_rhs_one(9,2,i) = -1./256. + (alpha)/128.
	f_rhs_one(10,2,i) = 0.
	f_rhs_one(11,2,i) = 0.
	
	!point-4 coefficients
	f_rhs_one(1,3,i) = 1./256. - (alpha)/128.
	f_rhs_one(2,3,i) = -1./32. + (alpha)/16.
	f_rhs_one(3,3,i) = 7./64. + (25.*alpha)/32.
	f_rhs_one(4,3,i) = 25./32. + (7.*alpha)/16.
	f_rhs_one(5,3,i) = 35./128. + (29.*alpha)/64.
	f_rhs_one(6,3,i) = -7./32. + (7.*alpha)/16.
	f_rhs_one(7,3,i) = 7./64. - (7.*alpha)/32.
	f_rhs_one(8,3,i) = -1./32. + (alpha)/16.
	f_rhs_one(9,3,i) = 1./256. - (alpha)/128.
	f_rhs_one(10,3,i) = 0.
	f_rhs_one(11,3,i) = 0.
      else if(f_order == 6) then
	!point-2 coefficients
	f_rhs_one(1,1,i) = 1./64. + (31.*alpha)/32.
	f_rhs_one(2,1,i) = 29./32. + (3.*alpha)/16.
	f_rhs_one(3,1,i) = 15./64. + (17.*alpha)/32.
	f_rhs_one(4,1,i) = -5./16. + (5.*alpha)/8.
	f_rhs_one(5,1,i) = 15./64. - (15.*alpha)/32.
	f_rhs_one(6,1,i) = -3./32. + (3.*alpha)/16.
	f_rhs_one(7,1,i) = 1./64. - (alpha)/32.
	f_rhs_one(8,1,i) = 0.
	f_rhs_one(9,1,i) = 0.
	f_rhs_one(10,1,i) = 0.
	f_rhs_one(11,1,i) = 0.
	
	!point-3 coefficients
	f_rhs_one(1,2,i) = -1./64. + (alpha)/32.
	f_rhs_one(2,2,i) = 3./32. + (13.*alpha)/16.
	f_rhs_one(3,2,i) = 49./64. + (15.*alpha)/32.
	f_rhs_one(4,2,i) = 5./16. + (3.*alpha)/8.
	f_rhs_one(5,2,i) = -15./64. + (15.*alpha)/32.
	f_rhs_one(6,2,i) = 3./32. - (3.*alpha)/16.
	f_rhs_one(7,2,i) = -1./64. + (alpha)/32.
	f_rhs_one(8,2,i) = 0.
	f_rhs_one(9,2,i) = 0.
	f_rhs_one(10,2,i) = 0.
	f_rhs_one(11,2,i) = 0.
      else if(f_order == 4) then
	!point-2 coefficients
	f_rhs_one(1,1,i) = 1./16. + 7.*alpha/8.
	f_rhs_one(2,1,i) = 3./4. + alpha/2.
	f_rhs_one(3,1,i) = 3./8. + alpha/4.
	f_rhs_one(4,1,i) =-1./4. + alpha/2.
	f_rhs_one(5,1,i) = 1./16. - alpha/8.
	f_rhs_one(6,1,i) = 0. 
	f_rhs_one(7,1,i) = 0.
	f_rhs_one(8,1,i) = 0.
	f_rhs_one(9,1,i) = 0.
	f_rhs_one(10,1,i) = 0.
	f_rhs_one(11,1,i) = 0.
      endif

    endif
  enddo  
    
    if(nclx) then
      f_triAx = f_alphaX
      f_triBx = 1.
      f_triCx = f_alphaX
      
      f_triBx(1) = 2.*f_triBx(1)
      f_triBx(nx) = f_triBx(nx) + f_alphaX*f_alphaX
      
      f_triAx(1) = 0.
      f_triCx(nx) = 0.
    else
      f_triAx = f_alphaX
      f_triBx = 1.
      f_triCx = f_alphaX
      
      f_triAx(1) = 0.
      f_triAx(2) = 0.
      f_triCx(nx) = 0.
      f_triCx(nx-1) = 0.
    endif
    
    if(ncly) then
      f_triAy = f_alphaY
      f_triBy = 1.
      f_triCy = f_alphaY
      
      f_triBy(1) = 2.*f_triBy(1)
      f_triBy(ny) = f_triBy(ny) + f_alphaY*f_alphaY
      
      f_triAy(1) = 0.
      f_triCy(ny) = 0.
    else
      f_triAy = f_alphaY
      f_triBy = 1.
      f_triCy = f_alphaY
      
      f_triAy(1) = 0.
      f_triAy(2) = 0.
      f_triCy(ny) = 0.
      f_triCy(ny-1) = 0.
    endif
    
    if(nclz) then
      f_triAz = f_alphaZ
      f_triBz = 1.
      f_triCz = f_alphaZ
      
      f_triBz(1) = 2.*f_triBz(1)
      f_triBz(nz) = f_triBz(nz) + f_alphaZ*f_alphaZ
      
      f_triAz(1) = 0.
      f_triCz(nz) = 0.
    else
      f_triAz = f_alphaZ
      f_triBz = 1.
      f_triCz = f_alphaZ
      
      f_triAz(1) = 0.
      f_triAz(2) = 0.
      f_triCz(nz) = 0.
      f_triCz(nz-1) = 0.
    endif
    
    call filterSetUp(f_triAx,f_triBx,f_triCx,f_triA_Bpx,f_tri_Bpx,nx,nclx)
    call filterSetUp(f_triAy,f_triBy,f_triCy,f_triA_Bpy,f_tri_Bpy,ny,ncly)
    call filterSetUp(f_triAz,f_triBz,f_triCz,f_triA_Bpz,f_tri_Bpz,nz,nclz)
    
  end subroutine triFilterSchemes

subroutine filterSetUp(triA,triB,triC,triA_Bp,tri_Bp,n,isPeriodic)
    use decomp_2d, only : mytype
    implicit none
    integer :: i,n
    logical,intent(IN) :: isPeriodic
    real(mytype), dimension(n) :: triA,triB,triC,triA_Bp,tri_Bp
    
    if(isPeriodic)then
      tri_Bp = triB
      do i=2,n
	triA_Bp(i)=triA(i)/tri_Bp(i-1)
	tri_Bp(i)=tri_Bp(i)-triC(i-1)*triA_Bp(i)
      enddo
      tri_Bp = 1.0_mytype/tri_Bp
    else
      tri_Bp = triB
      do i=3,n-1
	triA_Bp(i)=triA(i)/tri_Bp(i-1)
	tri_Bp(i)=tri_Bp(i)-triC(i-1)*triA_Bp(i)
      enddo
      tri_Bp = 1.0_mytype/tri_Bp
    endif
  end subroutine filterSetUp
  
  subroutine skipBoundaryIndices()
 !This code is used in case you need to ensure that any loop indice
 !does not touch the boundaries(for non periodic domains) 
 !This function returns the actual indices that must be used in the 
 !loops which run over the entire domain, except the boundaries.
 !This is used for example, in filters and in flux calculation.
 !This gives values only for x-pencils
    
    !X-Pencils
    bstiX = 2
    beniX = xsize(1)-1
    bstjX = 1
    benjX = xsize(2)
    bstkX = 1
    benkX = xsize(3)
    
    if(nclx) then
    bstiX = 1
    beniX = xsize(1)
    endif
    
    if(xstart(2)==1) then
      if(ncly) then
      bstjX=1
      else
      bstjX=2
      endif
    endif
    
    if(xend(2)==ny_global) then
      if(ncly) then
      benjX=xsize(2)
      else
      benjX=xsize(2)-1
      endif
    endif
    
    if(xstart(3)==1) then
      if(nclz) then
      bstkX=1
      else
      bstkX=2
      endif
    endif
    
    if(xend(3)==nz_global) then
      if(nclz) then
      benkX=xsize(3)
      else
      benkX=xsize(3)-1
      endif
    endif
    
    
    !Y-Pencils
    bstiY = 1
    beniY = ysize(1)
    bstjY = 2
    benjY = ysize(2)-1
    bstkY = 1
    benkY = ysize(3)
    
    if(ncly) then
    bstjY = 1
    benjY = ysize(2)
    endif
    
    if(ystart(1)==1) then
      if(nclx) then
      bstiY=1
      else
      bstiY=2
      endif
    endif
    
    if(yend(1)==nx_global) then
      if(nclx) then
      beniY=ysize(1)
      else
      beniY=ysize(1)-1
      endif
    endif
    
    if(ystart(3)==1) then
      if(nclz) then
      bstkY=1
      else
      bstkY=2
      endif
    endif
    
    if(yend(3)==nz_global) then
      if(nclz) then
      benkY=ysize(3)
      else
      benkY=ysize(3)-1
      endif
    endif
    
    !Z-Pencils
    bstiZ = 1
    beniZ = zsize(1)
    bstjZ = 1
    benjZ = zsize(2)
    bstkZ = 2
    benkZ = zsize(3)-1
    
    if(nclz) then
    bstkZ = 1
    benkZ = zsize(3)
    endif
    
    if(zstart(1)==1) then
      if(nclx) then
      bstiZ=1
      else
      bstiZ=2
      endif
    endif
    
    if(zend(1)==nx_global) then
      if(nclx) then
      beniZ=zsize(1)
      else
      beniZ=zsize(1)-1
      endif
    endif
    
    if(zstart(2)==1) then
      if(ncly) then
      bstjZ=1
      else
      bstjZ=2
      endif
    endif
    
    if(zend(2)==ny_global) then
      if(ncly) then
      benjZ=zsize(2)
      else
      benjZ=zsize(2)-1
      endif
    endif
   
  end subroutine skipBoundaryIndices
  
    subroutine findCorrectedPeriodicIndices(i,nx)
  !This function is used to set the indices for periodic cases 
  !where the calling indice is greater than the domain size
    integer,intent(IN) :: nx
    integer,intent(INOUT) :: i
    
    if(i.gt.nx) i=i-nx
    if(i.lt.1) i=nx+i
  
 end subroutine findCorrectedPeriodicIndices
  
end module trifilter
