module buffer

USE variables, only: nx,ny,nz
USE decomp_2d
USE param

implicit none

! Rohan's Filter Parameters used buffer
real(mytype),save :: b_alphaX,b_alphaY,b_alphaZ
integer,save :: b_oneSideState,b_stateX,b_stateY,b_stateZ
integer,save :: b_order
integer,save :: b_leftSurfState,b_rightSurfState,b_topSurfState, &
                     b_bottomSurfState,b_frontSurfState,b_backSurfState

integer,save :: bstiX,beniX,bstjX,benjX,bstkX,benkX
integer,save :: bstiY,beniY,bstjY,benjY,bstkY,benkY
integer,save :: bstiZ,beniZ,bstjZ,benjZ,bstkZ,benkZ

integer,save :: b_stiX,b_stjX,b_stkX,b_eniX,b_enjX,b_enkX
integer,save :: b_stiY,b_stjY,b_stkY,b_eniY,b_enjY,b_enkY
integer,save :: b_stiZ,b_stjZ,b_stkZ,b_eniZ,b_enjZ,b_enkZ

! Rohan's filter used for buffer
real(mytype),save,allocatable,dimension(:,:) :: b_rhs2,b_rhs4,b_rhs6,b_rhs8,b_rhs10 
real(mytype),allocatable,save, dimension(:) :: b_triAx,b_triAy,b_triAz
real(mytype),allocatable,save, dimension(:) :: b_triBx,b_triBy,b_triBz
real(mytype),allocatable,save, dimension(:) :: b_triCx,b_triCy,b_triCz
real(mytype),allocatable,save, dimension(:) :: b_triA_Bpx,b_tri_Bpx
real(mytype),allocatable,save, dimension(:) :: b_triA_Bpy,b_tri_Bpy
real(mytype),allocatable,save, dimension(:) :: b_triA_Bpz,b_tri_Bpz
!Stores one sided filter coefficients
real(mytype),save,allocatable,dimension(:,:,:) :: b_rhs_one !second indice indicates the point number away from the


contains

!#########################
 subroutine setupBuffer_p
  
  call allocBuffVars_p
  call bufftriFilterSchemes
  call buffskipBoundaryIndices
  call buffsetFilterIndices

 end subroutine setupBuffer_p
!#########################

 
!#########################
 subroutine allocBuffVars_p


    if(b_order == 10) then
      allocate(b_rhs2(2,3),b_rhs4(3,3),b_rhs6(4,3),b_rhs8(5,3),b_rhs10(6,3))
    else if(b_order == 8) then
      allocate(b_rhs2(2,3),b_rhs4(3,3),b_rhs6(4,3),b_rhs8(5,3))
    else if(b_order == 6) then
      allocate(b_rhs2(2,3),b_rhs4(3,3),b_rhs6(4,3))
    else if(b_order == 4) then
      allocate(b_rhs2(2,3),b_rhs4(3,3))
    else if(b_order == 2) then
      allocate(b_rhs2(2,3))
    endif

    if(b_oneSideState == 1) then
      if(b_order == 10) then
        allocate(b_rhs_one(11,4,3))
      else if(b_order == 8) then
        allocate(b_rhs_one(11,3,3))
      else if(b_order == 6) then
        allocate(b_rhs_one(11,2,3))
      else if(b_order == 4) then
        allocate(b_rhs_one(11,1,3))
      endif
    endif

    allocate(b_triAx(nx),b_triBx(nx),b_triCx(nx))
    allocate(b_triAy(ny),b_triBy(ny),b_triCy(ny))
    allocate(b_triAz(nz),b_triBz(nz),b_triCz(nz))
    allocate(b_triA_Bpx(nx),b_tri_Bpx(nx))
    allocate(b_triA_Bpy(ny),b_tri_Bpy(ny))
    allocate(b_triA_Bpz(nz),b_tri_Bpz(nz))

!#########################
  end subroutine allocBuffVars_p
!#########################
                                                  


subroutine bufffilterRHS(rhs,fX,nx,ny,nz,sti,eni,stj,enj,stk,enk,isPeriodic,dir)
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
	  do l=0,b_order/2
	    
	    call buffcalcFilterRHSCoeff(b_order,l+1,dir,actualCoeff)
	    
	    if(dir == 1) then!X-direction
	      ipl = i+l
	      iml = i-l
	      call bufffindCorrectedPeriodicIndices(ipl,nx) 
	      call bufffindCorrectedPeriodicIndices(iml,nx)
	      rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(ipl,j,k)+fX(iml,j,k))
	    else if(dir == 2) then!Y-direction
	      ipl = j+l
	      iml = j-l
	      call bufffindCorrectedPeriodicIndices(ipl,ny) 
	      call bufffindCorrectedPeriodicIndices(iml,ny)
	      rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,ipl,k)+fX(i,iml,k))
	    else if(dir == 3) then!Z-direction
	      ipl = k+l
	      iml = k-l
	      call bufffindCorrectedPeriodicIndices(ipl,nz) 
	      call bufffindCorrectedPeriodicIndices(iml,nz)
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
	  if(b_oneSideState == 0) then !order reduces near the boundary
	    
	    call buffcalcActualOrder(actualOrder,i,j,k,sti,eni,stj,enj,stk,enk,dir)
	    
	    do l=0,actualOrder/2
		
		call buffcalcFilterRHSCoeff(actualOrder,l+1,dir,actualCoeff)
		
		if(dir == 1) then!X-direction
		  ipl = i+l
		  iml = i-l
		  rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(ipl,j,k)+fX(iml,j,k))
		  if((i==sti).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaX*fX(i-1,j,k)
		  if((i==eni).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaX*fX(i+1,j,k)
		else if(dir == 2) then!Y-direction
		  ipl = j+l
		  iml = j-l
		  rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,ipl,k)+fX(i,iml,k))
		  if((j==stj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaY*fX(i,j-1,k)
		  if((j==enj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaY*fX(i,j+1,k)
		else if(dir == 3) then!Z-direction
		  ipl = k+l
		  iml = k-l
		  rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,j,ipl)+fX(i,j,iml))
		  if((k==stk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaZ*fX(i,j,k-1)
		  if((k==enk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaZ*fX(i,j,k+1)
		endif
	      enddo  
	    else !One sided schemes used near the boundary

	      do l=0,b_order/2
	
		call buffcalcFilterRHSCoeff(b_order,l+1,dir,actualCoeff)
		
		if(dir == 1) then!X-direction
		  ipl = i+l
		  iml = i-l
		  
		  if(l==0) call buffcalcOneSidedFilterIndices(sti,eni,b_order,fOneStX,fOneEnX)
		  
		  if((i.lt.fOneStX).and.(l==0)) then
		    rhs(i,j,k) = b_rhs_one(1,i-1,dir)*fX(1,j,k) + b_rhs_one(2,i-1,dir)*fX(2,j,k) + b_rhs_one(3,i-1,dir)*fX(3,j,k) + &
				 b_rhs_one(4,i-1,dir)*fX(4,j,k) + b_rhs_one(5,i-1,dir)*fX(5,j,k) + b_rhs_one(6,i-1,dir)*fX(6,j,k) + &
				 b_rhs_one(7,i-1,dir)*fX(7,j,k) + b_rhs_one(8,i-1,dir)*fX(8,j,k) + b_rhs_one(9,i-1,dir)*fX(9,j,k) + &
				 b_rhs_one(10,i-1,dir)*fX(10,j,k) + b_rhs_one(11,i-1,dir)*fX(11,j,k)
		    if(i==sti) rhs(i,j,k) = rhs(i,j,k) - b_alphaX*fX(i-1,j,k) 
		  else if((i.gt.fOneEnX).and.(l==0)) then
		    rhs(i,j,k) = b_rhs_one(1,eni-i+1,dir)*fX(eni+1,j,k) + b_rhs_one(2,eni-i+1,dir)*fX(eni,j,k) + &
				 b_rhs_one(3,eni-i+1,dir)*fX(eni-1,j,k) + b_rhs_one(4,eni-i+1,dir)*fX(eni-2,j,k) + &
				 b_rhs_one(5,eni-i+1,dir)*fX(eni-3,j,k) + b_rhs_one(6,eni-i+1,dir)*fX(eni-4,j,k) + &
				 b_rhs_one(7,eni-i+1,dir)*fX(eni-5,j,k) + b_rhs_one(8,eni-i+1,dir)*fX(eni-6,j,k) + &
				 b_rhs_one(9,eni-i+1,dir)*fX(eni-7,j,k) + b_rhs_one(10,eni-i+1,dir)*fX(eni-8,j,k) + &
				 b_rhs_one(11,eni-i+1,dir)*fX(eni-9,j,k)
		    if(i==eni) rhs(i,j,k) = rhs(i,j,k) - b_alphaX*fX(i+1,j,k)
		  else if((i.ge.fOneStX).and.(i.le.fOneEnX)) then
		    rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(ipl,j,k)+fX(iml,j,k))
		    if((i==eni).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaX*fX(i+1,j,k)!needed for b_order == 2
		    if((i==sti).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaX*fX(i-1,j,k)!needed for b_order == 2
		  endif
		else if(dir == 2) then!Y-direction
		  ipl = j+l
		  iml = j-l
		  
		  if(l==0) call buffcalcOneSidedFilterIndices(stj,enj,b_order,fOneStX,fOneEnX)
		  
		  if((j.lt.fOneStX).and.(l==0)) then
		    rhs(i,j,k) = b_rhs_one(1,j-1,dir)*fX(i,1,k) + b_rhs_one(2,j-1,dir)*fX(i,2,k) + b_rhs_one(3,j-1,dir)*fX(i,3,k) + &
				 b_rhs_one(4,j-1,dir)*fX(i,4,k) + b_rhs_one(5,j-1,dir)*fX(i,5,k) + b_rhs_one(6,j-1,dir)*fX(i,6,k) + &
				 b_rhs_one(7,j-1,dir)*fX(i,7,k) + b_rhs_one(8,j-1,dir)*fX(i,8,k) + b_rhs_one(9,j-1,dir)*fX(i,9,k) + &
				 b_rhs_one(10,j-1,dir)*fX(i,10,k) + b_rhs_one(11,j-1,dir)*fX(i,11,k)
		  if(j==stj) rhs(i,j,k) = rhs(i,j,k) - b_alphaY*fX(i,j-1,k) 
		  else if((j.gt.fOneEnX).and.(l==0)) then
		    rhs(i,j,k) = b_rhs_one(1,enj-j+1,dir)*fX(i,enj+1,k) + b_rhs_one(2,enj-j+1,dir)*fX(i,enj,k) + &
				 b_rhs_one(3,enj-j+1,dir)*fX(i,enj-1,k) + b_rhs_one(4,enj-j+1,dir)*fX(i,enj-2,k) + &
				 b_rhs_one(5,enj-j+1,dir)*fX(i,enj-3,k) + b_rhs_one(6,enj-j+1,dir)*fX(i,enj-4,k) + &
				 b_rhs_one(7,enj-j+1,dir)*fX(i,enj-5,k) + b_rhs_one(8,enj-j+1,dir)*fX(i,enj-6,k) + &
				 b_rhs_one(9,enj-j+1,dir)*fX(i,enj-7,k) + b_rhs_one(10,enj-j+1,dir)*fX(i,enj-8,k) + &
				 b_rhs_one(11,enj-j+1,dir)*fX(i,enj-9,k)
		    if(j==enj) rhs(i,j,k) = rhs(i,j,k) - b_alphaY*fX(i,j+1,k)
		  else if((j.ge.fOneStX).and.(j.le.fOneEnX)) then
		    rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,ipl,k)+fX(i,iml,k))
		    if((j==enj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaY*fX(i,j+1,k)!needed for b_order == 2
		    if((j==stj).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaY*fX(i,j-1,k)!needed for b_order == 2 
		  endif
		else if(dir == 3) then!Z-direction
		  ipl = k+l
		  iml = k-l
		  
		  if(l==0) call buffcalcOneSidedFilterIndices(stk,enk,b_order,fOneStX,fOneEnX)
		  
		  if((k.lt.fOneStX).and.(l==0)) then
		    rhs(i,j,k) = b_rhs_one(1,k-1,dir)*fX(i,j,1) + b_rhs_one(2,k-1,dir)*fX(i,j,2) + b_rhs_one(3,k-1,dir)*fX(i,j,3) + &
				 b_rhs_one(4,k-1,dir)*fX(i,j,4) + b_rhs_one(5,k-1,dir)*fX(i,j,5) + b_rhs_one(6,k-1,dir)*fX(i,j,6) + &
				 b_rhs_one(7,k-1,dir)*fX(i,j,7) + b_rhs_one(8,k-1,dir)*fX(i,j,8) + b_rhs_one(9,k-1,dir)*fX(i,j,9) + &
				 b_rhs_one(10,k-1,dir)*fX(i,j,10) + b_rhs_one(11,k-1,dir)*fX(i,j,11)
		    if(k==stk) rhs(i,j,k) = rhs(i,j,k) - b_alphaZ*fX(i,j,k-1) 
		  else if((k.gt.fOneEnX).and.(l==0)) then
		    rhs(i,j,k) = b_rhs_one(1,enk-k+1,dir)*fX(i,j,enk+1) + b_rhs_one(2,enk-k+1,dir)*fX(i,j,enk) + &
				 b_rhs_one(3,enk-k+1,dir)*fX(i,j,enk-1) + b_rhs_one(4,enk-k+1,dir)*fX(i,j,enk-2) + &
				 b_rhs_one(5,enk-k+1,dir)*fX(i,j,enk-3) + b_rhs_one(6,enk-k+1,dir)*fX(i,j,enk-4) + &
				 b_rhs_one(7,enk-k+1,dir)*fX(i,j,enk-5) + b_rhs_one(8,enk-k+1,dir)*fX(i,j,enk-6) + &
				 b_rhs_one(9,enk-k+1,dir)*fX(i,j,enk-7) + b_rhs_one(10,enk-k+1,dir)*fX(i,j,enk-8) + &
				 b_rhs_one(11,enk-k+1,dir)*fX(i,j,enk-9)
		    if(k==enk) rhs(i,j,k) = rhs(i,j,k) - b_alphaZ*fX(i,j,k+1)
		  else if((k.ge.fOneStX).and.(k.le.fOneEnX)) then
		    rhs(i,j,k) = rhs(i,j,k) + 0.5*actualCoeff*(fX(i,j,ipl)+fX(i,j,iml))
		    if((k==enk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaZ*fX(i,j,k+1)!needed for b_order == 2
		    if((k==stk).and.(l==0)) rhs(i,j,k) = rhs(i,j,k) - b_alphaZ*fX(i,j,k-1)!needed for b_order == 2 
		  endif
		endif
	      
	      enddo
	    
	    endif
	  
	enddo
      enddo
    enddo
  endif
  
 end subroutine bufffilterRHS
 
 subroutine bufffilterPerturbRHS(pert,nx,ny,nz,sti,eni,stj,enj,stk,enk,dir)
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
	pert(nx,j,k) = b_alphaX
      enddo
    enddo  
  else if(dir == 2) then
    do k=stk,enk 
      do i=sti,eni
	pert(i,1,k) = -1.
	do j=2,ny-1
	  pert(i,j,k) = 0.
	enddo
	pert(i,ny,k) = b_alphaY
      enddo
    enddo
  else if(dir == 3) then
    do j=stj,enj 
      do i=sti,eni
	pert(i,j,1) = -1.
	do k=2,nz-1
	  pert(i,j,k) = 0.
	enddo
	pert(i,j,nz) = b_alphaZ
      enddo
    enddo
  endif
 end subroutine bufffilterPerturbRHS
 
 subroutine buffcalcOneSidedFilterIndices(stX,enX,actualOrder,b_stX,b_enX)
  !stX and enX are boundary corrected parameters i.e. bstiX,beniX,bstjY,benjY etc
  integer,intent(IN) :: stX,enX,actualOrder
  integer,intent(OUT) :: b_stX,b_enX
  
  if(actualOrder == 2) then
    b_stX = stX
    b_enX = enX
  else if(actualOrder == 4) then
    b_stX = stX+1
    b_enX = enX-1
  else if(actualOrder == 6) then
    b_stX = stX+2
    b_enX = enX-2
  else if(actualOrder == 8) then
    b_stX = stX+3
    b_enX = enX-3
  else if(actualOrder == 10) then
    b_stX = stX+4
    b_enX = enX-4
  endif
  
 end subroutine buffcalcOneSidedFilterIndices
 
  subroutine buffcalcFilterRHSCoeff(actualOrder,indice,dir,actualCoeff)
  integer,intent(IN) :: actualOrder,indice,dir
  real(mytype),intent(OUT) :: actualCoeff

  if(actualOrder == 2) then
    actualCoeff = b_rhs2(indice,dir)
  else if(actualOrder == 4) then
    actualCoeff = b_rhs4(indice,dir)
  else if(actualOrder == 6) then
    actualCoeff = b_rhs6(indice,dir)
  else if(actualOrder == 8) then
    actualCoeff = b_rhs8(indice,dir)
  else if(actualOrder == 10) then
    actualCoeff = b_rhs10(indice,dir)  
  endif

 end subroutine buffcalcFilterRHSCoeff
  
 subroutine buffcalcActualOrder(actualOrder,i,j,k,sti,eni,stj,enj,stk,enk,dir)
  integer,intent(IN) :: dir,sti,eni,stj,enj,stk,enk,i,j,k 
  integer,intent(OUT) :: actualOrder
  
  if (b_order == 2) then
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
      actualOrder=b_order
    endif
  else if(b_order == 4) then
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
      actualOrder=b_order
    endif
  else if(b_order == 6) then
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
      actualOrder=b_order
    endif
  else if(b_order == 8) then
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
      actualOrder=b_order
    endif
  else if(b_order == 10) then
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
      actualOrder=b_order
    endif
  endif
    
 end subroutine buffcalcActualOrder
 

 subroutine buffsetFilterIndices()
  !Sets indices for the filters. Used when any of the surfaces have to be
  !filtered as well.
! call buffskipBoundaryIndices()	
  call buffcalcActualFilterIndices(b_stiX,b_stjX,b_stkX,b_eniX,b_enjX,b_enkX,1)
  call buffcalcActualFilterIndices(b_stiY,b_stjY,b_stkY,b_eniY,b_enjY,b_enkY,2)
  call buffcalcActualFilterIndices(b_stiZ,b_stjZ,b_stkZ,b_eniZ,b_enjZ,b_enkZ,3)
 end subroutine buffsetFilterIndices
 
  subroutine buffcalcActualFilterIndices(b_sti,b_stj,b_stk,b_eni,b_enj,b_enk,dir)
  
  integer,intent(OUT) :: b_sti,b_stj,b_stk,b_eni,b_enj,b_enk
  integer,intent(IN) :: dir
  
  if(dir == 1) then
    b_sti = bstiX
    b_stj = bstjX
    b_stk = bstkX
    b_eni = beniX
    b_enj = benjX
    b_enk = benkX
    
    if(.not. ncly) then
      if((b_topSurfState == 1).and.(xend(2) == ny_global)) b_enj = xsize(2)
      if((b_bottomSurfState == 1).and.(xstart(2) == 1)) b_stj = 1
    endif
    if(.not. nclz) then
      if((b_frontSurfState == 1).and.(xend(3) == nz_global)) b_enk = xsize(3)
      if((b_backSurfState == 1).and.(xstart(3) == 1)) b_stk = 1
    endif
  else if(dir == 2) then
    b_sti = bstiY
    b_stj = bstjY
    b_stk = bstkY
    b_eni = beniY
    b_enj = benjY
    b_enk = benkY
    
    if(.not. nclx) then
      if((b_rightSurfState == 1).and.(yend(1) == nx_global)) b_eni = ysize(1)
      if((b_leftSurfState == 1).and.(ystart(1) == 1)) b_sti = 1
    endif
    if(.not. nclz) then
      if((b_frontSurfState == 1).and.(yend(3) == nz_global)) b_enk = ysize(3)
      if((b_backSurfState == 1).and.(ystart(3) == 1)) b_stk = 1
    endif
    
  else if(dir == 3) then
    b_sti = bstiZ
    b_stj = bstjZ
    b_stk = bstkZ
    b_eni = beniZ
    b_enj = benjZ
    b_enk = benkZ
      
    if(.not. nclx) then
      if((b_rightSurfState == 1).and.(zend(1) == nx_global)) b_eni = zsize(1)
      if((b_leftSurfState == 1).and.(zstart(1) == 1)) b_sti = 1
    endif
    if(.not. ncly) then
      if((b_topSurfState == 1).and.(zend(2) == ny_global)) b_enj = zsize(2)
      if((b_bottomSurfState == 1).and.(zstart(2) == 1)) b_stj = 1
    endif
  endif
  
 end subroutine buffcalcActualFilterIndices
   
subroutine bufftriInvFilX(tx,ux,rx,sx,triC,triA_Bp,tri_Bp,nx,ny,nz,isPeriodic)
  
  integer :: nx,ny,nz,i,j,k 
  real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
  real(mytype), dimension(ny,nz):: sx
  real(mytype), dimension(nx):: triC,triA_Bp,tri_Bp 
  logical,intent(IN) :: isPeriodic
    
  if(isPeriodic) then
    call bufffilterRHS(tx,ux,nx,ny,nz,b_stiX,b_eniX,b_stjX,b_enjX,b_stkX,b_enkX,isPeriodic,1)
    call bufffilterPerturbRHS(rx,nx,ny,nz,b_stiX,b_eniX,b_stjX,b_enjX,b_stkX,b_enkX,1)
    do k=b_stkX,b_enkX
      do j=b_stjX,b_enjX 
	
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
	  sx(j,k)=(tx(1,j,k)-b_alphaX*tx(nx,j,k))&
	      /(1.+rx(1,j,k)-b_alphaX*rx(nx,j,k)) 
	  do i=1,nx 
	    tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
	  enddo
	  
      enddo
    enddo
! do k=1,nz 
!     do j=1,ny 
! 	tx(1,j,k)=b_rhs(1)*ux(1,j,k) + 0.5*b_rhs(2)*(ux(2,j,k)+ux(nx,j,k)) &
! 		  + 0.5*b_rhs(3)*(ux(3,j,k)+ux(nx-1,j,k))
! 	rx(1,j,k)=-1. 
! 	tx(2,j,k)=b_rhs(1)*ux(2,j,k) + 0.5*b_rhs(2)*(ux(3,j,k)+ux(1,j,k)) &
! 		  + 0.5*b_rhs(3)*(ux(4,j,k)+ux(nx,j,k))
! 	rx(2,j,k)=0. 
! 	do i=3,nx-2
! 	  tx(i,j,k)=b_rhs(1)*ux(i,j,k) + 0.5*b_rhs(2)*(ux(i+1,j,k)+ux(i-1,j,k)) &
! 		  + 0.5*b_rhs(3)*(ux(i+2,j,k)+ux(i-2,j,k))
! 	  rx(i,j,k)=0. 
! 	enddo
! 	tx(nx-1,j,k)=b_rhs(1)*ux(nx-1,j,k) + 0.5*b_rhs(2)*(ux(nx,j,k)+ux(nx-2,j,k)) &
! 		  + 0.5*b_rhs(3)*(ux(1,j,k)+ux(nx-3,j,k))
! 	rx(nx-1,j,k)=0. 
! 	tx(nx,j,k)=b_rhs(1)*ux(nx,j,k) + 0.5*b_rhs(2)*(ux(1,j,k)+ux(nx-1,j,k)) &
! 		  + 0.5*b_rhs(3)*(ux(2,j,k)+ux(nx-2,j,k)) 
! 	rx(nx,j,k)=b_alpha           
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
! 	sx(j,k)=(tx(1,j,k)-b_alpha*tx(nx,j,k))&
! 	    /(1.+rx(1,j,k)-b_alpha*rx(nx,j,k)) 
! 	do i=1,nx 
! 	  tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
! 	enddo
!     enddo
!     enddo
   else
    call bufffilterRHS(tx,ux,nx,ny,nz,b_stiX,b_eniX,b_stjX,b_enjX,b_stkX,b_enkX,isPeriodic,1)
    do k=b_stkX,b_enkX
      do j=b_stjX,b_enjX
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
!     do k=b_stk,b_enk
!       do j=b_stj,b_enj
! 	  tx(1,j,k)=ux(1,j,k) 
! ! 	  tx(2,j,k)= b_rhs_red1(1)*ux(2,j,k) + 0.5*b_rhs_red1(2)*(ux(1,j,k)+ux(3,j,k)) &
! ! 			-b_alpha*ux(1,j,k)
! 	  tx(2,j,k)= b_one(1)*ux(1,j,k) + b_one(2)*ux(2,j,k) + b_one(3)*ux(3,j,k) + &
! 		     b_one(4)*ux(4,j,k) + b_one(5)*ux(5,j,k) - b_alpha*ux(1,j,k) 
! 	  do i=3,nx-2
! 	    tx(i,j,k)=b_rhs(1)*ux(i,j,k) + 0.5*b_rhs(2)*(ux(i+1,j,k)+ux(i-1,j,k)) &
! 		  + 0.5*b_rhs(3)*(ux(i+2,j,k)+ux(i-2,j,k))
! 	  enddo
! ! 	  tx(nx-1,j,k)= b_rhs_red1(1)*ux(nx-1,j,k) + 0.5*b_rhs_red1(2)*(ux(nx-2,j,k)+ux(nx,j,k)) &
! ! 			-b_alpha*ux(nx,j,k)
! 	  tx(nx-1,j,k)= b_one(1)*ux(nx,j,k) + b_one(2)*ux(nx-1,j,k) + b_one(3)*ux(nx-2,j,k) + &
! 		     b_one(4)*ux(nx-3,j,k) + b_one(5)*ux(nx-4,j,k) - b_alpha*ux(nx,j,k)
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

  end subroutine bufftriInvFilX
  
 subroutine bufftriInvFilY(ty,uy,ry,sy,triC,triA_Bp,tri_Bp,nx,ny,nz,isPeriodic) 
  
  integer :: nx,ny,nz,i,j,k
  real(mytype), dimension(nx,ny,nz) :: ty,uy 
  real(mytype), dimension(nx,ny,nz) :: ry
  real(mytype), dimension(nx,nz)  :: sy
  real(mytype), dimension(ny) :: triC,triA_Bp,tri_Bp 
  logical,intent(IN) :: isPeriodic
    
  if(isPeriodic) then
    call bufffilterRHS(ty,uy,nx,ny,nz,b_stiY,b_eniY,b_stjY,b_enjY,b_stkY,b_enkY,isPeriodic,2)
    call bufffilterPerturbRHS(ry,nx,ny,nz,b_stiY,b_eniY,b_stjY,b_enjY,b_stkY,b_enkY,2)
   do k=b_stkY,b_enkY
    do j=2,ny  
      do i=b_stiY,b_eniY 
	  ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*triA_Bp(j) 
	  ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*triA_Bp(j) 
      enddo
    enddo 
   enddo
   do k=b_stkY,b_enkY
    do i=b_stiY,b_eniY
	ty(i,ny,k)=ty(i,ny,k)*tri_Bp(ny) 
	ry(i,ny,k)=ry(i,ny,k)*tri_Bp(ny) 
    enddo
   enddo 
   do k=b_stkY,b_enkY
    do j=ny-1,1,-1  
      do i=b_stiY,b_eniY  
	  ty(i,j,k)=(ty(i,j,k)-triC(j)*ty(i,j+1,k))*tri_Bp(j) 
	  ry(i,j,k)=(ry(i,j,k)-triC(j)*ry(i,j+1,k))*tri_Bp(j) 
      enddo
    enddo 
   enddo
   do k=b_stkY,b_enkY
    do i=b_stiY,b_eniY  
	sy(i,k)=(ty(i,1,k)-b_alphaY*ty(i,ny,k))&
	    /(1.+ry(i,1,k)-b_alphaY*ry(i,ny,k)) 
    enddo
   enddo 
   do k=b_stkY,b_enkY
    do j=1,ny  
      do i=b_stiY,b_eniY  
	  ty(i,j,k)=ty(i,j,k)-sy(i,k)*ry(i,j,k) 
      enddo
    enddo
   enddo
 
  else
!     tmp = beniY
!     if(tmp == ysize(1)-1) tmp=tmp+1
    call bufffilterRHS(ty,uy,nx,ny,nz,b_stiY,b_eniY,b_stjY,b_enjY,b_stkY,b_enkY,isPeriodic,2)
    do k=b_stkY,b_enkY
      do j=3,ny-1  
	do i=b_stiY,b_eniY  
	    ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*triA_Bp(j) 
	enddo 
      enddo 
    enddo 
    do k=b_stkY,b_enkY
      do i=b_stiY,b_eniY  
	  ty(i,ny-1,k)=ty(i,ny-1,k)*tri_Bp(ny-1) 
      enddo 
    enddo 
    do k=b_stkY,b_enkY
      do j=ny-2,2,-1  
	do i=b_stiY,b_eniY  
	    ty(i,j,k)=(ty(i,j,k)-triC(j)*ty(i,j+1,k))*tri_Bp(j) 
	enddo 
      enddo 
    enddo
  endif
 end subroutine bufftriInvFilY
  
subroutine bufftriInvFilZ(tz,uz,rz,sz,triC,triA_Bp,tri_Bp,nx,ny,nz,isPeriodic) 

  integer :: nx,ny,nz,i,j,k
  real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
  real(mytype), dimension(nx,ny) :: sz 
  real(mytype), dimension(nz) :: triC,triA_Bp,tri_Bp
  logical,intent(IN) :: isPeriodic
    
  if(isPeriodic) then
    call bufffilterRHS(tz,uz,nx,ny,nz,b_stiZ,b_eniZ,b_stjZ,b_enjZ,b_stkZ,b_enkZ,isPeriodic,3)
    call bufffilterPerturbRHS(rz,nx,ny,nz,b_stiZ,b_eniZ,b_stjZ,b_enjZ,b_stkZ,b_enkZ,3)
    do k=2,nz
      do j=b_stjZ,b_enjZ
	do i=b_stiZ,b_eniZ
	    tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*triA_Bp(k)
	    rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*triA_Bp(k)
	enddo
      enddo
    enddo
    do j=b_stjZ,b_enjZ
      do i=b_stiZ,b_eniZ
	  tz(i,j,nz)=tz(i,j,nz)*tri_Bp(nz)
	  rz(i,j,nz)=rz(i,j,nz)*tri_Bp(nz)
      enddo
    enddo
    do k=nz-1,1,-1
      do j=b_stjZ,b_enjZ
	do i=b_stiZ,b_eniZ
	    tz(i,j,k)=(tz(i,j,k)-triC(k)*tz(i,j,k+1))*tri_Bp(k)
	    rz(i,j,k)=(rz(i,j,k)-triC(k)*rz(i,j,k+1))*tri_Bp(k)
	enddo
      enddo
    enddo
    do j=b_stjZ,b_enjZ
      do i=b_stiZ,b_eniZ
	  sz(i,j)=(tz(i,j,1)-b_alphaZ*tz(i,j,nz))/&
	      (1.+rz(i,j,1)-b_alphaZ*rz(i,j,nz))
      enddo
    enddo
    do k=1,nz
      do j=b_stjZ,b_enjZ
	do i=b_stiZ,b_eniZ
	    tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
	enddo
      enddo
    enddo
    
    else
      call bufffilterRHS(tz,uz,nx,ny,nz,b_stiZ,b_eniZ,b_stjZ,b_enjZ,b_stkZ,b_enkZ,isPeriodic,3)
      do k=3,nz-1
	do j=b_stjZ,b_enjZ
	  do i=b_stiZ,b_eniZ
	    tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*triA_Bp(k)
	  enddo
	enddo
      enddo
      do j=b_stjZ,b_enjZ
	do i=b_stiZ,b_eniZ
	  tz(i,j,nz-1)=tz(i,j,nz-1)*tri_Bp(nz-1)
	enddo
      enddo
      do k=nz-2,2,-1
	do j=b_stjZ,b_enjZ
	  do i=b_stiZ,b_eniZ
	    tz(i,j,k)=(tz(i,j,k)-triC(k)*tz(i,j,k+1))*tri_Bp(k)
	  enddo
	enddo
      enddo
    endif 
  end subroutine bufftriInvFilZ
 
 
 subroutine buffcomputeTriFilteredVar(unfilteredVar)
 !This subroutine buffcomputes the filtered variable in an MPI communication effective way
 
     USE var, ONLY : ta2
     USE var, ONLY : ta3
     USE var, ONLY : td1,di1 
     USE var, ONLY : tb2,di2
     USE var, ONLY : tb3,di3
     USE var, ONLY : sx,sy,sz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(INOUT) :: unfilteredVar
      
      if((b_stateX==1).and.(b_stateY==1).and.(b_stateZ==1)) then
	
	call bufftriInvFilX(td1,unfilteredVar,di1,sx,b_triCx,b_triA_Bpx,&
	      b_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	
	call transpose_x_to_y(td1,ta2)
	call bufftriInvFilY(tb2,ta2,di2,sy,b_triCy,b_triA_Bpy,&
	      b_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)

	call transpose_y_to_z(tb2,ta3)
	call bufftriInvFilZ(tb3,ta3,di3,sz,b_triCz,b_triA_Bpz,&
	      b_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
			
      else if((b_stateX==1).and.(b_stateY==1)) then
	
	call bufftriInvFilX(td1,unfilteredVar,di1,sx,b_triCx,b_triA_Bpx,&
	      b_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	
	call transpose_x_to_y(td1,ta2)
	call bufftriInvFilY(tb2,ta2,di2,sy,b_triCy,b_triA_Bpy,&
	      b_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)
	call transpose_y_to_x(tb2,unfilteredVar)      
      
      else if((b_stateX==1).and.(b_stateZ==1)) then
	
	call bufftriInvFilX(td1,unfilteredVar,di1,sx,b_triCx,b_triA_Bpx,&
	      b_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	
	call transpose_x_to_y(td1,tb2)
	call transpose_y_to_z(tb2,ta3)
	call bufftriInvFilZ(tb3,ta3,di3,sz,b_triCz,b_triA_Bpz,&
	      b_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
      
      else if((b_stateY==1).and.(b_stateZ==1)) then
	
	call transpose_x_to_y(unfilteredVar,ta2)
	call bufftriInvFilY(tb2,ta2,di2,sy,b_triCy,b_triA_Bpy,&
	      b_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)

	call transpose_y_to_z(tb2,ta3)
	call bufftriInvFilZ(tb3,ta3,di3,sz,b_triCz,b_triA_Bpz,&
	      b_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
      
      else if(b_stateX==1) then
      
      call bufftriInvFilX(td1,unfilteredVar,di1,sx,b_triCx,b_triA_Bpx,&
	      b_tri_Bpx,xsize(1),xsize(2),xsize(3),nclx)
	unfilteredVar = td1
      
      else if(b_stateY==1) then
	
	call transpose_x_to_y(unfilteredVar,ta2)
	call bufftriInvFilY(tb2,ta2,di2,sy,b_triCy,b_triA_Bpy,&
	      b_tri_Bpy,ysize(1),ysize(2),ysize(3),ncly)
	call transpose_y_to_x(tb2,unfilteredVar)
	
      else if(b_stateZ==1) then
	
	call transpose_x_to_y(unfilteredVar,tb2)
	call transpose_y_to_z(tb2,ta3)
	call bufftriInvFilZ(tb3,ta3,di3,sz,b_triCz,b_triA_Bpz,&
	      b_tri_Bpz,zsize(1),zsize(2),zsize(3),nclz)
	call transpose_z_to_y(tb3,tb2)
	call transpose_y_to_x(tb2,unfilteredVar)
      
      endif
                  
 end subroutine buffcomputeTriFilteredVar
 
  subroutine bufftriFilterSchemes()
    
   integer :: i 
   real(mytype) :: alpha
   
   do i=1,3
    
    if(i==1) alpha = b_alphaX
    if(i==2) alpha = b_alphaY
    if(i==3) alpha = b_alphaZ
    
    !filter coefficients for filters which reduce order near the boundary
    if(b_order .ge. 2) then
      b_rhs2(1,i) = 0.5 + alpha
      b_rhs2(2,i) = 0.5 + alpha
    endif
    
    if(b_order .ge. 4) then
      b_rhs4(1,i) = 5./8. + 3.*alpha/4.
      b_rhs4(2,i) = 0.5 + alpha
      b_rhs4(3,i) = -1./8. + alpha/4.
    endif
    
    if(b_order .ge. 6) then
      b_rhs6(1,i) = 11./16. + 5.*alpha/8.
      b_rhs6(2,i) = 15./32. + 17.*alpha/16.
      b_rhs6(3,i) = -3./16. + 3.*alpha/8.
      b_rhs6(4,i) = 1./32. - alpha/16.
    endif
    
    if(b_order .ge. 8) then
      b_rhs8(1,i) = (93.+70.*alpha)/128.
      b_rhs8(2,i) = (7.+18.*alpha)/16.
      b_rhs8(3,i) = (-7.+14.*alpha)/32.
      b_rhs8(4,i) = 1./16. - alpha/8.
      b_rhs8(5,i) = -1./128. + alpha/64.
    endif
    
    if(b_order .ge. 10) then
      b_rhs10(1,i) = (193.+126.*alpha)/256.
      b_rhs10(2,i) = (105.+302.*alpha)/256.
      b_rhs10(3,i) = 15.*(-1.+2.*alpha)/64.
      b_rhs10(4,i) = 45.*(1.-2.*alpha)/512.
      b_rhs10(5,i) = 5.*(-1.+2.*alpha)/256.
      b_rhs10(6,i) = (1.-2.*alpha)/512.
    endif
    
    !filter coefficients for one-sided formulation
    if(b_oneSideState == 1) then
      
      if(b_order == 10) then
	!point-2 coefficients
	b_rhs_one(1,1,i) = (1.+1022.*alpha)/1024.
	b_rhs_one(2,1,i) = (507.+10.*alpha)/512.
	b_rhs_one(3,1,i) = (45.+934.*alpha)/1024.
	b_rhs_one(4,1,i) = 15.*(-1.+2.*alpha)/128.
	b_rhs_one(5,1,i) = 105.*(1.-2.*alpha)/512.
	b_rhs_one(6,1,i) = 63.*(-1.+2.*alpha)/256.
	b_rhs_one(7,1,i) = 105.*(1.-2.*alpha)/512.
	b_rhs_one(8,1,i) = 15.*(-1.+2.*alpha)/128.
	b_rhs_one(9,1,i) = 45.*(1.-2.*alpha)/1024.
	b_rhs_one(10,1,i) = 5.*(-1.+2.*alpha)/512.
	b_rhs_one(11,1,i) = (1.-2.*alpha)/1024.
	
	!point-3 coefficients
	b_rhs_one(1,2,i) = (-1.+2.*alpha)/1024.
	b_rhs_one(2,2,i) = (5.+502.*alpha)/512.
	b_rhs_one(3,2,i) = (979.+90.*alpha)/1024.
	b_rhs_one(4,2,i) = (15.+98.*alpha)/128.
	b_rhs_one(5,2,i) = 105.*(-1.+2.*alpha)/512.
	b_rhs_one(6,2,i) = 63.*(1.-2.*alpha)/256.
	b_rhs_one(7,2,i) = 105.*(-1.+2.*alpha)/512.
	b_rhs_one(8,2,i) = 15.*(1.-2.*alpha)/128.
	b_rhs_one(9,2,i) = 45.*(-1.+2.*alpha)/1024.
	b_rhs_one(10,2,i) = 5.*(1.-2.*alpha)/512.
	b_rhs_one(11,2,i) = (-1.+2.*alpha)/1024.
	
	!point-4 coefficients
	b_rhs_one(1,3,i) = (1.-2.*alpha)/1024.
	b_rhs_one(2,3,i) = 5.*(-1.+2.*alpha)/512.
	b_rhs_one(3,3,i) = (45.+934.*alpha)/1024.
	b_rhs_one(4,3,i) = (113.+30.*alpha)/128.
	b_rhs_one(5,3,i) = (105.+302.*alpha)/512.
	b_rhs_one(6,3,i) = 63.*(-1.+2.*alpha)/256.
	b_rhs_one(7,3,i) = 105.*(1.-2.*alpha)/512.
	b_rhs_one(8,3,i) = 15.*(-1.+2.*alpha)/128.
	b_rhs_one(9,3,i) = 45.*(1.-2.*alpha)/1024.
	b_rhs_one(10,3,i) = 5.*(-1.+2.*alpha)/512.
	b_rhs_one(11,3,i) = (1.-2.*alpha)/1024.
	
	!point-5 coefficients
	b_rhs_one(1,4,i) = (-1.+2.*alpha)/1024.
	b_rhs_one(2,4,i) = 5.*(1.-2.*alpha)/512.
	b_rhs_one(3,4,i) = 45.*(-1.+2.*alpha)/1024.
	b_rhs_one(4,4,i) = (15.+98.*alpha)/128.
	b_rhs_one(5,4,i) = (407.+210.*alpha)/512.
	b_rhs_one(6,4,i) = (63.+130.*alpha)/256.
	b_rhs_one(7,4,i) = 105.*(-1.+2.*alpha)/512.
	b_rhs_one(8,4,i) = 15.*(1.-2.*alpha)/128.
	b_rhs_one(9,4,i) = 45.*(-1.+2.*alpha)/1024.
	b_rhs_one(10,4,i) = 5.*(1.-2.*alpha)/512.
	b_rhs_one(11,4,i) = (-1.+2.*alpha)/1024.
      else if(b_order == 8) then
	!point-2 coefficients
	b_rhs_one(1,1,i) = 1./256. + (127.*alpha)/128.
	b_rhs_one(2,1,i) = 31./32. + (alpha)/16.
	b_rhs_one(3,1,i) = 7./64. + (25.*alpha)/32.
	b_rhs_one(4,1,i) = -7./32. + (7.*alpha)/16.
	b_rhs_one(5,1,i) = 7.*(5.-10*alpha)/128.
	b_rhs_one(6,1,i) = -7./32. + (7.*alpha)/16.
	b_rhs_one(7,1,i) = 7./64. - (7.*alpha)/32.
	b_rhs_one(8,1,i) = -1./32. + (alpha)/16.
	b_rhs_one(9,1,i) = 1./256. - (alpha)/128.
	b_rhs_one(10,1,i) = 0.
	b_rhs_one(11,1,i) = 0.
	
	!point-3 coefficients
	b_rhs_one(1,2,i) = -1./256. + (alpha)/128.
	b_rhs_one(2,2,i) = 1./32. + (15.*alpha)/16.
	b_rhs_one(3,2,i) = 57./64. + (7.*alpha)/32.
	b_rhs_one(4,2,i) = 7./32. + (9.*alpha)/16.
	b_rhs_one(5,2,i) = 7.*(-5.+10.*alpha)/128.
	b_rhs_one(6,2,i) = 7./32. - (7.*alpha)/16.
	b_rhs_one(7,2,i) = -7./64. + (7.*alpha)/32.
	b_rhs_one(8,2,i) = 1./32. - (alpha)/16.
	b_rhs_one(9,2,i) = -1./256. + (alpha)/128.
	b_rhs_one(10,2,i) = 0.
	b_rhs_one(11,2,i) = 0.
	
	!point-4 coefficients
	b_rhs_one(1,3,i) = 1./256. - (alpha)/128.
	b_rhs_one(2,3,i) = -1./32. + (alpha)/16.
	b_rhs_one(3,3,i) = 7./64. + (25.*alpha)/32.
	b_rhs_one(4,3,i) = 25./32. + (7.*alpha)/16.
	b_rhs_one(5,3,i) = 35./128. + (29.*alpha)/64.
	b_rhs_one(6,3,i) = -7./32. + (7.*alpha)/16.
	b_rhs_one(7,3,i) = 7./64. - (7.*alpha)/32.
	b_rhs_one(8,3,i) = -1./32. + (alpha)/16.
	b_rhs_one(9,3,i) = 1./256. - (alpha)/128.
	b_rhs_one(10,3,i) = 0.
	b_rhs_one(11,3,i) = 0.
      else if(b_order == 6) then
	!point-2 coefficients
	b_rhs_one(1,1,i) = 1./64. + (31.*alpha)/32.
	b_rhs_one(2,1,i) = 29./32. + (3.*alpha)/16.
	b_rhs_one(3,1,i) = 15./64. + (17.*alpha)/32.
	b_rhs_one(4,1,i) = -5./16. + (5.*alpha)/8.
	b_rhs_one(5,1,i) = 15./64. - (15.*alpha)/32.
	b_rhs_one(6,1,i) = -3./32. + (3.*alpha)/16.
	b_rhs_one(7,1,i) = 1./64. - (alpha)/32.
	b_rhs_one(8,1,i) = 0.
	b_rhs_one(9,1,i) = 0.
	b_rhs_one(10,1,i) = 0.
	b_rhs_one(11,1,i) = 0.
	
	!point-3 coefficients
	b_rhs_one(1,2,i) = -1./64. + (alpha)/32.
	b_rhs_one(2,2,i) = 3./32. + (13.*alpha)/16.
	b_rhs_one(3,2,i) = 49./64. + (15.*alpha)/32.
	b_rhs_one(4,2,i) = 5./16. + (3.*alpha)/8.
	b_rhs_one(5,2,i) = -15./64. + (15.*alpha)/32.
	b_rhs_one(6,2,i) = 3./32. - (3.*alpha)/16.
	b_rhs_one(7,2,i) = -1./64. + (alpha)/32.
	b_rhs_one(8,2,i) = 0.
	b_rhs_one(9,2,i) = 0.
	b_rhs_one(10,2,i) = 0.
	b_rhs_one(11,2,i) = 0.
      else if(b_order == 4) then
	!point-2 coefficients
	b_rhs_one(1,1,i) = 1./16. + 7.*alpha/8.
	b_rhs_one(2,1,i) = 3./4. + alpha/2.
	b_rhs_one(3,1,i) = 3./8. + alpha/4.
	b_rhs_one(4,1,i) =-1./4. + alpha/2.
	b_rhs_one(5,1,i) = 1./16. - alpha/8.
	b_rhs_one(6,1,i) = 0. 
	b_rhs_one(7,1,i) = 0.
	b_rhs_one(8,1,i) = 0.
	b_rhs_one(9,1,i) = 0.
	b_rhs_one(10,1,i) = 0.
	b_rhs_one(11,1,i) = 0.
      endif

    endif
  enddo  
    
    if(nclx) then
      b_triAx = b_alphaX
      b_triBx = 1.
      b_triCx = b_alphaX
      
      b_triBx(1) = 2.*b_triBx(1)
      b_triBx(nx) = b_triBx(nx) + b_alphaX*b_alphaX
      
      b_triAx(1) = 0.
      b_triCx(nx) = 0.
    else
      b_triAx = b_alphaX
      b_triBx = 1.
      b_triCx = b_alphaX
      
      b_triAx(1) = 0.
      b_triAx(2) = 0.
      b_triCx(nx) = 0.
      b_triCx(nx-1) = 0.
    endif
    
    if(ncly) then
      b_triAy = b_alphaY
      b_triBy = 1.
      b_triCy = b_alphaY
      
      b_triBy(1) = 2.*b_triBy(1)
      b_triBy(ny) = b_triBy(ny) + b_alphaY*b_alphaY
      
      b_triAy(1) = 0.
      b_triCy(ny) = 0.
    else
      b_triAy = b_alphaY
      b_triBy = 1.
      b_triCy = b_alphaY
      
      b_triAy(1) = 0.
      b_triAy(2) = 0.
      b_triCy(ny) = 0.
      b_triCy(ny-1) = 0.
    endif
    
    if(nclz) then
      b_triAz = b_alphaZ
      b_triBz = 1.
      b_triCz = b_alphaZ
      
      b_triBz(1) = 2.*b_triBz(1)
      b_triBz(nz) = b_triBz(nz) + b_alphaZ*b_alphaZ
      
      b_triAz(1) = 0.
      b_triCz(nz) = 0.
    else
      b_triAz = b_alphaZ
      b_triBz = 1.
      b_triCz = b_alphaZ
      
      b_triAz(1) = 0.
      b_triAz(2) = 0.
      b_triCz(nz) = 0.
      b_triCz(nz-1) = 0.
    endif
    
    call bufffilterSetUp(b_triAx,b_triBx,b_triCx,b_triA_Bpx,b_tri_Bpx,nx,nclx)
    call bufffilterSetUp(b_triAy,b_triBy,b_triCy,b_triA_Bpy,b_tri_Bpy,ny,ncly)
    call bufffilterSetUp(b_triAz,b_triBz,b_triCz,b_triA_Bpz,b_tri_Bpz,nz,nclz)
    
  end subroutine bufftriFilterSchemes

subroutine bufffilterSetUp(triA,triB,triC,triA_Bp,tri_Bp,n,isPeriodic)
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
  end subroutine bufffilterSetUp
  
  subroutine buffskipBoundaryIndices()
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
    
    if (nclx) then
    bstiX = 1
    beniX = xsize(1)
    endif
    
    if(xstart(2)==1) then
      if (ncly) then
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
   
  end subroutine buffskipBoundaryIndices
  
    subroutine bufffindCorrectedPeriodicIndices(i,nx)
  !This function is used to set the indices for periodic cases 
  !where the calling indice is greater than the domain size
    integer,intent(IN) :: nx
    integer,intent(INOUT) :: i
    
    if(i.gt.nx) i=i-nx
    if(i.lt.1) i=nx+i
  
 end subroutine bufffindCorrectedPeriodicIndices
  
end module buffer
