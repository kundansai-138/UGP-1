!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module jet_PM

  use decomp_2d
  use variables
  use param
  !!$use complex_geometry, only: tol

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_jet_PM, boundary_conditions_jet_PM, visu_jet_PM, visu_jet_PM_init

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: init_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Initial laminar conditions in pipe
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine init_jet_PM (ux1,uy1,uz1,phi1)
  !
  !********************************************************************
    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
    use dbg_schemes, only: exp_prec
	use ibm_param
	
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
	
	! re5strt problem
	real(mytype) :: ym,zm,r,um,yc,zc
	
	real(mytype) :: ux_p,uy_p,uz_p,uth_p    						!@@6A velocity components
	
	! Make changes in ecoule_td too
	
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
		
	   call ecoule_pm

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

  end subroutine init_jet_PM
  !********************************************************************
  !
  subroutine boundary_conditions_jet_PM (ux,uy,uz,phi)
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
  end subroutine boundary_conditions_jet_PM
  
  subroutine ecoule_pm
	USE param
	use variables
    USE decomp_2d
    USE ibm_param
	
	use param_swirling_PM
	
	implicit none
	
	integer  :: j,k
	
	! re5strt problem
	real(mytype) :: ym,zm,rad,yc,zc
	
	real(mytype) :: ux_p,uy_p,uz_p,uth_p    						!@@6A velocity components
	
	if (nrank==0) then 
		write(*,*) 'alp = ', alp_pm
		write(*,*) 'del = ', del_pm
		write(*,*) 'Swirl Number = ', S_pm
	endif
	
	yc = yly / two
    zc = zlz / two
	
	bxz1=zero
	bxy1=zero
	bxz1=zero
	
	! re5strt problem
	do k=1,xsize(3)
	      zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
		do j=1,xsize(2)
	         ym=real(j+xstart(2)-1-1,mytype)*dy-yc
		     rad=sqrt(ym*ym+zm*zm)
		     !****************** Maxworthy *************************
		     ux_p = 1d0 - (alp_pm-1d0)/(2d0*alp_pm)*(1d0+erf((rad-1d0)/del_pm))
		     uth_p = S_pm*rad/2d0*(1d0 - erf((rad-1d0)/del_pm))
			 
			 if (rad .ne. 0) then
			    uz_p =  uth_p*ym/rad 		!@@P Prev error correct, right-hand coord now!
			    uy_p = -uth_p*zm/rad 		!@@P y-z to r-th (uy=urcos-uthsin,uz=ursin+uthcos)
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
	
  end subroutine ecoule_pm
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
		call ecoule_pm
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
  !********************************************************************
  
  subroutine visu_jet_PM_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised
	
	!call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    !call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)
	
    !call decomp_2d_register_variable(io_name, "dux", 1, 0, output2D, mytype)
    !call decomp_2d_register_variable(io_name, "duy", 1, 0, output2D, mytype)
	!call decomp_2d_register_variable(io_name, "dvx", 1, 0, output2D, mytype)
	!call decomp_2d_register_variable(io_name, "dvy", 1, 0, output2D, mytype)
	!call decomp_2d_register_variable(io_name, "dwx", 1, 0, output2D, mytype)
	!call decomp_2d_register_variable(io_name, "dwy", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_jet_PM_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_cyl
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs cylinder-specific visualization
  !!
  !############################################################################
  subroutine visu_jet_PM(ux1, uy1, uz1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
	!du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
	!VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    !call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                  - td1(:,:,:)*tb1(:,:,:) &
                  - tg1(:,:,:)*tc1(:,:,:) &
                  - th1(:,:,:)*tf1(:,:,:)
    !call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush
	
	!call write_field(ta1, ".", "dux", num, flush = .true.)
	!call write_field(td1, ".", "duy", num, flush = .true.)
	!call write_field(tb1, ".", "dvx", num, flush = .true.)
	!call write_field(te1, ".", "dvy", num, flush = .true.)
	!call write_field(tc1, ".", "dwx", num, flush = .true.)
	!call write_field(tf1, ".", "dwy", num, flush = .true.)
	
  end subroutine visu_jet_PM
  
end module jet_PM
