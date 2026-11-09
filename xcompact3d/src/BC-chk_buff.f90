!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module chk_buff

  use decomp_2d
  use variables
  use param
  !!$use complex_geometry, only: tol

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_chk_buff, boundary_conditions_chk_buff

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
  subroutine init_chk_buff (ux1,uy1,uz1,phi1)
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

    real(mytype) :: ym,zm,r,um
    integer :: k,j,i,ii,is,code

    if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

    endif

    ux1=zero; uy1=zero; uz1=zero

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

       !modulation of the random noise + inlet velocity profile
       do k=1,xsize(3)
		  zm=dz*real(xstart(3)-1+k-1,mytype)-zlz/2d0
          do j=1,xsize(2)
             if (istret.eq.0) ym=(j+xstart(2)-1-1)*dy-yly/2d0
             if (istret.ne.0) ym=yp(j+xstart(2)-1)-yly/2d0
             r=sqrt(ym*ym+zm*zm)
			 um=exp_prec(-zptwo*r*r)
			 do i=1,xsize(1)
                ux1(i,j,k)=um*ux1(i,j,k) + u1
                uy1(i,j,k)=um*uy1(i,j,k)
                uz1(i,j,k)=um*uz1(i,j,k)
             enddo
          enddo
       enddo
    endif

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return

  end subroutine init_chk_buff
  !********************************************************************
  !
  subroutine boundary_conditions_chk_buff (ux,uy,uz,phi)
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
	
    return
    !
  end subroutine boundary_conditions_chk_buff
  
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
	
	real(mytype)                    :: r,ym,zm, yc,zc
	
	real(mytype) :: ux_p,uy_p,uz_p,uth_p    						!@@6A velocity components
	real(mytype), parameter :: alp_p=100d0, del_p=0.2d0, S_p=1.5d0		! Tolstoy Jet - Maxworthy
	
	yc = yly / two
    zc = zlz / two
	
    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
    
	do k=1,xsize(3)
       do j=1,xsize(2)
          bxx1(j,k)=u1+bxo(j,k)*inflow_noise
          bxy1(j,k)=zero+byo(j,k)*inflow_noise
          bxz1(j,k)=zero+bzo(j,k)*inflow_noise
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


end module chk_buff
