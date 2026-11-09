!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module jet_diverging

  use decomp_2d
  use variables
  use param
  !!$use complex_geometry, only: tol

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: geomcomplex_jet_diverging, init_jet_diverging, boundary_conditions_jet_diverging

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: geomcomplex_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !!      Modified by Tolstoy Debbarma
  !! DESCRIPTION: Generates epsilon matrix for pipe geometry
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine geomcomplex_jet_diverging(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  !
  !********************************************************************

    use decomp_2d,only : mytype
    use MPI
    use param,only : zero,one,two,yly,zlz
    use ibm_param
	use param_diverging
	
    implicit none

    integer                                         :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny)                      :: yp
    real(mytype)                                    :: dx,dz,remp,tol
    !LOCALS
    real(mytype)                                    :: r,ym,zm,yc,zc,xm
    integer                                         :: i,j,k,code,ierror
	
	real(mytype)                                    :: angle_rad, r_inner, r_outer
	
    epsi(:,:,:) = zero
    yc = yly/two
    zc = zlz/two
    tol=1e-15
	
	
	! Convert divergence angle from degrees to radians
	angle_rad = div_angle_degree * (acos(-1.0_mytype)/180.0_mytype)
  
    !safety check
    if (nrank == 0) then
       if (rai.le.0) then
           write(*,*) 'SIMULATION IS STOPPED!'
           write(*,*) 'Please specify a valid value for the pipe inner radius in input.i3d (rai)'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif
       if (rao.le.0) then
           write(*,*) 'SIMULATION IS STOPPED!'
           write(*,*) 'Please specify a valid value for the pipe outer radius in input.i3d (rao)'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif
    endif

    do k=nzi,nzf
        zm=real(k-1,mytype)*dz-zc
        do j=nyi,nyf
            ym=yp(j)-yc
            r=sqrt(ym*ym+zm*zm)
            do i=nxi,nxf
				xm=real(i-1,mytype)*dx
				
				! Calculate local inner and outer radius at position xm
				r_inner = rai + xm * tan(angle_rad)
				r_outer = rao + xm * tan(angle_rad)
			  
                if (r.gt.r_inner .and. r.lt.r_outer .and. xm.lt.len_div) then
					epsi(i,j,k)=remp
                elseif (abs(r-r_inner).lt.tol .and. abs(xm-len_div).lt.tol) then
                    epsi(i,j,k)=remp
                elseif (abs(r-r_outer).lt.tol .and. abs(xm-len_div).lt.tol) then
                    epsi(i,j,k)=remp
                endif
            enddo
        enddo
    enddo
	!
    return

  end subroutine geomcomplex_jet_diverging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: init_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Initial laminar conditions in pipe
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine init_jet_diverging (ux1,uy1,uz1,phi1)
  !
  !********************************************************************
    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
    use dbg_schemes, only: exp_prec
	use ibm_param
	
	use param_diverging
	
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
	
	! re5strt problem
	real(mytype) :: ym,zm,r,um,yc,zc
	
    integer :: k,j,i,ii,is,code
	
	yc = yly / two
    zc = zlz / two

    ux1=zero; uy1=zero; uz1=zero
	
	if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

    endif
	
    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)

       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux1(i,j,k)=init_noise*ux1(i,j,k)
                uy1(i,j,k)=init_noise*uy1(i,j,k)
                uz1(i,j,k)=init_noise*uz1(i,j,k)
             enddo
          enddo
       enddo
	   
	   call ecoule_div
       !modulation of the random noise + inlet velocity profile
       do k=1,xsize(3)
		  zm=dz*real(xstart(3)-1+k-1,mytype)-zlz/2d0
          do j=1,xsize(2)
             if (istret.eq.0) ym=(j+xstart(2)-1-1)*dy-yly/2d0
             if (istret.ne.0) ym=yp(j+xstart(2)-1)-yly/2d0
             r=sqrt(ym*ym+zm*zm)
			 um=exp_prec(-zptwo*r*r)
			 do i=1,xsize(1)
                ux1(i,j,k)=um*ux1(i,j,k)+bxx1(j,k) ! + u1
                uy1(i,j,k)=um*uy1(i,j,k)+bxy1(j,k)
                uz1(i,j,k)=um*uz1(i,j,k)+bxz1(j,k)
             enddo
          enddo
       enddo
    endif

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return

  end subroutine init_jet_diverging
  !********************************************************************
  !
  subroutine boundary_conditions_jet_diverging (ux,uy,uz,phi)
  !
  !********************************************************************

    use param
    use variables
    use decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    !
	
	call inflow (ux,uy,uz,phi)
    call outflow (ux,uy,uz,phi)
	
    return
    !
  end subroutine boundary_conditions_jet_diverging
  
  subroutine ecoule_div
	USE param
	use variables
    USE decomp_2d
    USE ibm_param
	
	use param_diverging
	
	implicit none
	
	integer  :: j,k
	
	! re5strt problem
	real(mytype) :: ym,zm,r,yc,zc
	
	real(mytype) :: ux_p,uy_p,uz_p,uth_p    						!@@6A velocity components
	real(mytype) :: sub_fr_r												! Tolstoy Jet - Maxworthy
	! Make changes in ecoule_td too
	
	if (nrank==0) then 
		write(*,*) 'alp_div = ', alp_div
		write(*,*) 'del_div = ', del_div
		write(*,*) 'del_cf_div = ', del_cf_div
		write(*,*) 'S_div = ', S_div
		write(*,*) 'div_angle_degree = ', div_angle_degree
		write(*,*) 'len_div = ', len_div
	endif
	yc = yly / two
    zc = zlz / two
	
	sub_fr_r=1d0-2d0*del_div
	
	bxz1=zero
	bxy1=zero
	bxz1=zero
	! re5strt problem
	do k=1,xsize(3)
	      zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
		do j=1,xsize(2)
	         ym=real(j+xstart(2)-1-1,mytype)*dy-yc
		     r=sqrt(ym*ym+zm*zm)
		     !****************** Maxworthy *************************
		     ux_p = 1d0 - 0.5d0*(1d0+erf((r-sub_fr_r)/del_div))
		     uth_p = S_div*(r/two)*(1d0 - erf((r-sub_fr_r)/del_div))
		     
			 if (r.gt.rao) then
				ux_p = ux_p+erf((r-rao)/del_cf_div)*alp_div
			 endif
			 
			 if (r .ne. 0) then
			    uz_p =  uth_p*ym/r 		!@@P Prev error correct, right-hand coord now!
			    uy_p = -uth_p*zm/r 		!@@P y-z to r-th (uy=urcos-uthsin,uz=ursin+uthcos)
		     else
			    uz_p = 0d0
			    uy_p = 0d0
		     endif
		     bxx1(j,k) = ux_p
		     bxy1(j,k) = uy_p
		     bxz1(j,k) = uz_p
		enddo
	enddo
	! Make changes in ecoule_td too
	
  end subroutine ecoule_div
  !********************************************************************
  subroutine inflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE ibm_param

    implicit none

    integer  :: j,k,is
	real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
	
    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
	
	if (irestart==1 .and. itime==ifirst) then
		call ecoule_div
	endif
	
	do k=1,xsize(3)
       do j=1,xsize(2)
		  bxx1(j,k) = bxx1(j,k)+bxo(j,k)*inflow_noise
		  bxy1(j,k) = bxy1(j,k)+bxo(j,k)*inflow_noise
		  bxz1(j,k) = bxz1(j,k)+bxo(j,k)*inflow_noise	
       enddo
    enddo
	
    if (iscalar.eq.1) then
       do is=1, numscalar
          do k=1,xsize(3)
             do j=1,xsize(2)
                phi(1,j,k,is)=cp(is)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine inflow
  !********************************************************************
  
  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm_param

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609._mytype
    uxmin=1609._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
       enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1 == one) then
       cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
    !elseif (u1 == one) then
    !   cx=uxmax1*gdt(itr)*udx
    elseif (u1 == two) then
       cx=u2*gdt(itr)*udx    !works better
    else
       cx=(half*(u1+u2))*gdt(itr)*udx
    endif

    do k=1,xsize(3)
       do j=1,xsize(2)
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (iscalar==1) then
       if (u2==zero) then
          cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
       elseif (u2==one) then
          cx=uxmax1*gdt(itr)*udx
       elseif (u2==two) then
          cx=u2*gdt(itr)*udx    !works better
       else
          stop
       endif

       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow
  
end module jet_diverging