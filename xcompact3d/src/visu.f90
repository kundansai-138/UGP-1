!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module visu
  
  implicit none

  ! True to activate the XDMF output
  logical, save :: use_xdmf = .true.
  ! True to use the new enumeration
  logical, save :: filenamedigits = .false.
  ! output2D is defined in the input.i3d file
  !        0 for 3D output (default)
  !        1 for 2D output with X average
  !        2 for 2D output with Y average
  !        3 for 2D output with Z average
  integer, save :: output2D
  integer :: ioxdmf
  character(len=9) :: ifilenameformat = '(I3.3)'
  real, save :: tstart, tend

  character(len=*), parameter :: io_name = "solution-io"

  private
  public :: output2D, visu_init, visu_ready, visu_finalise, write_snapshot, end_snapshot, &
       write_field, io_name, gen_filename

contains

  !
  ! Initialize the visu module
  !
  subroutine visu_init()

    use MPI
    use param, only : ilmn, iscalar, ilast, ifirst, ioutput, istret
    use variables, only : numscalar, prec, nvisu
    use param, only : dx, dy, dz
    use decomp_2d, only : nrank, mytype, xszV, yszV, zszV, xsize, ysize, zsize
    use decomp_2d_io, only : decomp_2d_init_io, decomp_2d_open_io, decomp_2d_append_mode
    use decomp_2d_io, only : decomp_2d_register_variable

    
    implicit none

    ! Local variables
    integer :: noutput, nsnapout
    real(mytype) :: memout

    integer :: is

    ! HDD usage of visu module
    if (nrank==0) then
      noutput = (ilast - ifirst +1)/ioutput

      nsnapout = 4
      if (ilmn)         nsnapout = nsnapout + 1
      if (iscalar.ne.0) nsnapout = nsnapout + numscalar

      memout = prec * nsnapout * noutput
      if (output2D.eq.0) then
        memout = memout * xszV(1) * yszV(2) * zszV(3)
      else if (output2D.eq.1) then
        memout = memout *           yszV(2) * zszV(3)
      else if (output2D.eq.2) then
        memout = memout * xszV(1)           * zszV(3)
      else if (output2D.eq.3) then
        memout = memout * xszV(1) * yszV(2)
      endif
      write(*,*)'==========================================================='
      write(*,*)'Visu module requires ',real(memout*1e-9,4),'GB'
      write(*,*)'==========================================================='
    end if

    ! Safety check
    if (output2D < 0 .or. output2D > 3 &
        .or. (output2d == 2.and.istret /= 0)) then
      if (nrank.eq.0) write(*,*) "Visu module: incorrect value for output2D."
      call MPI_ABORT(MPI_COMM_WORLD, 0, noutput)
      stop
    endif

    call decomp_2d_init_io(io_name)

    !! Register variables
    call decomp_2d_register_variable(io_name, "ux", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "uy", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "uz", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "pp", 1, 0, output2D, mytype)
    if (ilmn) then
       call decomp_2d_register_variable(io_name, "rho", 1, 0, output2D, mytype)
    endif
    if (iscalar.ne.0) then
       do is = 1, numscalar
          call decomp_2d_register_variable(io_name, "phi"//char(48+is), 1, 0, output2D, mytype)
       enddo
    endif
    
  end subroutine visu_init

  !
  ! Indicate visu ready for IO. This is only required for ADIOS2 backend, using MPIIO this
  ! subroutine doesn't do anything.
  ! XXX: Call after all visu initialisation (main visu + case visu)
  !
  subroutine visu_ready ()

    use decomp_2d_io, only : decomp_2d_open_io, decomp_2d_append_mode, decomp_2d_write_mode, gen_iodir_name

    use param, only : irestart
    
    implicit none

    integer :: mode
    
#ifdef ADIOS2
    logical, save :: outloc_init
    logical :: dir_exists
    
    mode = decomp_2d_write_mode

    ! XXX: A fix was applied to ADIOS2.7.1 series to prevent corrupting files when appended from Fortran
    if (.not.outloc_init) then
       if (irestart == 1) then
          !! Restarting - is the output already available to write to?
          inquire(file=gen_iodir_name("data", io_name), exist=dir_exists)
          if (dir_exists) then
             outloc_init = .true.
          end if
       end if
       
       if (.not.outloc_init) then !! Yes, yes, but check the restart check above.
          mode = decomp_2d_write_mode
       else
          mode = decomp_2d_append_mode
       end if
       outloc_init = .true.
    else
       mode = decomp_2d_append_mode
    end if

    call decomp_2d_open_io(io_name, "data", mode)
#endif
    
  end subroutine visu_ready
  
  !
  ! Finalise the visu module. When using the ADIOS2 backend this closes the IO which is held open
  ! for the duration of the simulation, otherwise does nothing.
  ! 
  subroutine visu_finalise()

    use decomp_2d_io, only : decomp_2d_close_io
    
    implicit none

#ifdef ADIOS2
    call decomp_2d_close_io(io_name, "data")
#endif
    
  end subroutine visu_finalise

  !
  ! Write a snapshot
  !
  subroutine write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime, num)
	
	use mpi
	use decomp_2d, only : transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_start_io

    use param, only : nrhotime, ilmn, iscalar, ioutput, irestart
	use param_p, only : isubdom, subdom_together, nx_first, nx_last, ny_first, ny_last, nz_first, nz_last
	use param_p, only : iplanez0, iplaney0

    use variables, only : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : numscalar

    use var, only : pp1, ta1, di1, nxmsize
    use var, only : pp2, ppi2, dip2, ph2, nymsize
    use var, only : ppi3, dip3, ph3, nzmsize
    use var, only : npress

    use tools, only : rescale_pressure

    implicit none
	
	integer :: ierr_mpi
    !! inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize,npress), intent(in) :: pp3
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    integer, intent(in) :: itime
    integer, intent(out) :: num

    ! Local variables
    integer :: is
    integer :: ierr
    character(len=30) :: scname
    integer :: mode
    logical, save :: outloc_init = .false.
    logical :: dir_exists

    ! Update log file
    if (nrank.eq.0) then
      call cpu_time(tstart)
      print *,'Writing snapshots =>',itime/ioutput
    end if

#ifdef ADIOS2
    call decomp_2d_start_io(io_name, "data")
#endif
    
    ! Snapshot number
#ifndef ADIOS2
    if (filenamedigits) then
       ! New enumeration system, it works integrated with xcompact3d_toolbox
       num = itime
    else
       ! Classic enumeration system
       num = itime / ioutput
    endif
#else
    ! ADIOS2 is zero-indexed
    num = itime/ioutput - 1
#endif
    
    ! Write XDMF header
    if (use_xdmf) call write_xdmf_header(".", "snapshot", num)

    ! Write velocity
    call write_field(ux1, ".", "ux", num)
	call write_field(uy1, ".", "uy", num)  ! comment nsm
	call write_field(uz1, ".", "uz", num)  ! comment nsm
	
	if (isubdom==1) then
		call write_field_subdom(ux1, "data_subdom", "ux_subdom", num)
		call write_field_subdom(uy1, "data_subdom", "uy_subdom", num)
		call write_field_subdom(uz1, "data_subdom", "uz_subdom", num)
    endif
	
	if (iplanez0==1) then
		call write_field_subdom(ux1, "iplanez0", "ux_z0", num)
		call write_field_subdom(uy1, "iplanez0", "uy_z0", num)
		call write_field_subdom(uz1, "iplanez0", "uz_z0", num)
    endif
	
	if (iplaney0==1) then
		call write_field_subdom(ux1, "iplaney0", "ux_y0", num)
		call write_field_subdom(uy1, "iplaney0", "uy_y0", num)
		call write_field_subdom(uz1, "iplaney0", "uz_y0", num)
    endif
	
    ! Interpolate pressure
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
    (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)


    ! Rescale pressure
    call rescale_pressure(ta1)

    ! Write pressure
    call write_field(ta1, ".", "pp", num, .true., flush=.true.)  ! comment nsm
	if (isubdom==1) then
		call write_field_subdom(ta1, "data_subdom", "pp_subdom", num, .true., flush=.true.)
    endif

	if (iplanez0==1) then
		call write_field_subdom(ta1, "iplanez0", "pp_z0", num, .true., flush=.true.)
    endif

	if (iplaney0==1) then
		call write_field_subdom(ta1, "iplaney0", "pp_y0", num, .true., flush=.true.)
    endif
	
    ! LMN - write density
    if (ilmn) call write_field(rho1(:,:,:,1), ".", "rho", num)
	
    ! Write scalars
    if (iscalar.ne.0) then
      do is = 1, numscalar
        write(scname,"('phi',I2.2)") is
        call write_field(phi1(:,:,:,is), ".", trim(scname), num, .true.)
      enddo
    endif
	
	! Ensure all ranks have finished writing subdomain files
	call MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
	! Only rank 0 merges and cleans up
	if (nrank == 0 .and. isubdom==1 .and. subdom_together==1) then
		print *, 'DEBUG: calling merge_and_cleanup_subdom_files for num =', num
		call merge_and_cleanup_subdom_files("data_subdom", num, nx_last-nx_first+1, ny_last-ny_first+1, nz_last-nz_first+1)
		print *, 'DEBUG: merge_and_cleanup_subdom_files returned'
	end if

  end subroutine write_snapshot

  subroutine end_snapshot(itime, num)

    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_end_io
    use param, only : istret, xlx, yly, zlz, last_snapshot_only
    use param_p, only : iSFD
    use variables, only : nx, ny, nz, beta
    use var, only : dt,t

    implicit none
	
	integer :: ierr_mpi
	logical :: ex_new
	character(len=:), allocatable :: fnew, fold
	character(len=:), allocatable :: fnew2, fold2
	character(len=:), allocatable :: fnew3, fold3
	character(len=:), allocatable :: fnew4, fold4
	character(len=:), allocatable :: fnew5, fold5
	character(len=:), allocatable :: fnew6, fold6
	character(len=:), allocatable :: fnew7, fold7
	character(len=:), allocatable :: fnew8, fold8
	character(len=:), allocatable :: fnew9, fold9
	character(len=:), allocatable :: fnew10, fold10
	character(len=:), allocatable :: fnew11, fold11
	character(len=:), allocatable :: fnew12, fold12
	character(len=:), allocatable :: fnew13, fold13
	
    integer, intent(in) :: itime
    integer, intent(in) :: num

    character(len=:), allocatable :: filename
    
    character(len=32) :: fmt2, fmt3, fmt4
    integer :: is
    integer :: ierr
    
    ! Write XDMF footer
    if (use_xdmf) call write_xdmf_footer()

    ! Add metadata if XDMF is not used
    if (.not. use_xdmf) then
      if (nrank==0) then

        write(fmt2,'("(A,I16)")')
        write(fmt3,'("(A,F16.4)")')
        write(fmt4,'("(A,F16.12)")')

        filename = "./data/snap"//int_to_str(num)//".ini"
        open(newunit=is,file=filename,action='write',status='replace')
        write(is,'(A)')'[domain]'
        write(is,fmt2) 'nx=      ',nx
        write(is,fmt2) 'ny=      ',ny
        write(is,fmt2) 'nz=      ',nz
        write(is,fmt2) 'istret=  ',istret
        write(is,fmt4) 'beta=    ',beta
        write(is,fmt3) 'Lx=      ',xlx
        write(is,fmt3) 'Ly=      ',yly
        write(is,fmt3) 'Lz=      ',zlz
        write(is,'(A)')'[time]'
        write(is,fmt2) 'itime=   ',itime
        write(is,fmt3) 'dt=      ',dt
        write(is,fmt3) 't =      ',t
        close(is)

      endif
    endif

#ifdef ADIOS2
    call decomp_2d_end_io(io_name, "data")
#endif
    
    ! Update log file
    if (nrank.eq.0) then
      call cpu_time(tend)
      write(*,'(" Time for writing snapshots (s): ",F12.8)') tend-tstart
    endif

	if (last_snapshot_only==1) then
		if (nrank == 0 .and. num > 0) then
			! If your files are in ./data (common with this module), use:
			fnew = "./data/" // gen_filename(".", "ux", num  , "bin")
			fold = "./data/" // gen_filename(".", "ux", num-1, "bin")
		
			fnew2 = "./data/" // gen_filename(".", "uy", num  , "bin")
			fold2 = "./data/" // gen_filename(".", "uy", num-1, "bin")
		
			fnew3 = "./data/" // gen_filename(".", "uz", num  , "bin")
			fold3 = "./data/" // gen_filename(".", "uz", num-1, "bin")
		
			fnew4 = "./data/" // gen_filename(".", "pp", num  , "bin")
			fold4 = "./data/" // gen_filename(".", "pp", num-1, "bin")
		
			fnew5 = "./data/" // gen_filename(".", "vort", num  , "bin")
			fold5 = "./data/" // gen_filename(".", "vort", num-1, "bin")
		
			fnew6 = "./data/" // gen_filename(".", "critq", num  , "bin")
			fold6 = "./data/" // gen_filename(".", "critq", num-1, "bin")
		
			fnew7 = "./data/" // gen_filename(".", "snapshot", num  , "xdmf")
			fold7 = "./data/" // gen_filename(".", "snapshot", num-1, "xdmf")
			
			if (iSFD==1) then
				fnew8 = "./data/" // gen_filename(".", "dux", num  , "bin")
				fold8 = "./data/" // gen_filename(".", "dux", num-1, "bin")

				fnew9 = "./data/" // gen_filename(".", "duy", num  , "bin")
				fold9 = "./data/" // gen_filename(".", "duy", num-1, "bin")
		
				fnew10 = "./data/" // gen_filename(".", "dvx", num  , "bin")
				fold10 = "./data/" // gen_filename(".", "dvx", num-1, "bin")
		
				fnew11 = "./data/" // gen_filename(".", "dvy", num  , "bin")
				fold11 = "./data/" // gen_filename(".", "dvy", num-1, "bin")
		
				fnew12 = "./data/" // gen_filename(".", "dwx", num  , "bin")
				fold12 = "./data/" // gen_filename(".", "dwx", num-1, "bin")
		
				fnew13 = "./data/" // gen_filename(".", "dwy", num  , "bin")
				fold13 = "./data/" // gen_filename(".", "dwy", num-1, "bin")
			end if
		
			inquire(file=trim(fnew), exist=ex_new)
			if (ex_new) then
				call execute_command_line('rm -f "'//trim(fold)//'"')
				call execute_command_line('rm -f "'//trim(fold2)//'"')
				call execute_command_line('rm -f "'//trim(fold3)//'"')
				call execute_command_line('rm -f "'//trim(fold4)//'"')
				call execute_command_line('rm -f "'//trim(fold5)//'"')
				call execute_command_line('rm -f "'//trim(fold6)//'"')
				call execute_command_line('rm -f "'//trim(fold7)//'"')
				if (iSFD==1) then
					call execute_command_line('rm -f "'//trim(fold8)//'"')
					call execute_command_line('rm -f "'//trim(fold9)//'"')
					call execute_command_line('rm -f "'//trim(fold10)//'"')
					call execute_command_line('rm -f "'//trim(fold11)//'"')
					call execute_command_line('rm -f "'//trim(fold12)//'"')
					call execute_command_line('rm -f "'//trim(fold13)//'"')
				end if
			end if
		end if
	end if

  end subroutine end_snapshot

  !
  ! Write the header of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_header(pathname, filename, num)

    use variables, only : nvisu, yp
    use param, only : dx,dy,dz,istret, xlx, yly, zlz
    use decomp_2d, only : mytype, nrank, xszV, yszV, zszV, ystV

    implicit none

    ! Arguments
    character(len=*), intent(in) :: pathname, filename
    integer, intent(in) :: num

    ! Local variables
    integer :: i,k
    real(mytype) :: xp(xszV(1)), zp(zszV(3))

    character(len=:), allocatable :: fmt
    
    if (nrank.eq.0) then
      OPEN(newunit=ioxdmf,file="./data/"//gen_snapshotname(pathname, filename, num, "xdmf"))

      write(ioxdmf,'(A22)')'<?xml version="1.0" ?>'
      write(ioxdmf,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(ioxdmf,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(ioxdmf,*)'<Domain>'
      call write_xdmf_topo()
      if (istret.ne.0) then
        do i=1,xszV(1)
          xp(i) = real(i-1,mytype)*dx*nvisu
        enddo
		do i=1,yszV(1)
          yp(i) = real(i-1,mytype)*dy*nvisu 			!- yly/2d0
        enddo
        do k=1,zszV(3)
          zp(k) = real(k-1,mytype)*dz*nvisu 			!- zlz/2d0
        enddo
        write(ioxdmf,*)'    <Geometry name="geo" Type="VXVYVZ">'
        if (output2D.ne.1) then
          write(ioxdmf,*)'        <DataItem Dimensions="',xszV(1),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xp(:)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        if (output2D.ne.2) then
          write(ioxdmf,*)'        <DataItem Dimensions="',yszV(2),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',yp(ystV(1)::nvisu)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',yp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        if (output2D.ne.3) then
          write(ioxdmf,*)'        <DataItem Dimensions="',zszV(3),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',zp(:)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',zp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      else
        write(ioxdmf,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
        write(ioxdmf,*)'        <!-- Origin -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        write(ioxdmf,*)'        ',-zlz/2d0," ",-yly/2d0," ",0.0						!	0.0 0.0 0.0'
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'        <!-- DxDyDz -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        if (mytype == kind(0.0d0)) then
           fmt = "(A, E24.17, A, E24.17, A, E24.17)"
        else
           fmt = "(A, E16.9, A, E16.9, A, E16.9)"
        end if
        if (output2D.eq.0) then
          write(ioxdmf,fmt)'        ',nvisu*dz," ",nvisu*dy," ",nvisu*dx
        else if (output2D.eq.1) then
          write(ioxdmf,fmt)'        ',dz," ",dy," ",1.
        else if (output2D.eq.2) then
          write(ioxdmf,fmt)'        ',dz," ",1.," ",dx
        else if (output2D.eq.3) then
          write(ioxdmf,fmt)'        ',1.," ",dy," ",dx
        endif
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      endif
      write(ioxdmf,'(A, I0, A)')'    <Grid Name="', num, '" GridType="Uniform">'
      write(ioxdmf,*)'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(ioxdmf,*)'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    endif
  end subroutine write_xdmf_header

  subroutine write_xdmf_topo()

    use decomp_2d, only : xszV, yszV, zszV
    use param, only : istret
    
    implicit none

    character(len=:), allocatable :: topo_type
    character(len=:), allocatable :: fmt
    
    if (istret /= 0) then
       topo_type = "3DRectMesh"
    else
       topo_type = "3DCoRectMesh"
    end if
    
    write(ioxdmf,'(A)')'    <Topology name="topo" TopologyType="'//topo_type//'"'

    fmt = "(A, I0, A, I0, A, I0, A)"
    if (output2D.eq.0) then
       write(ioxdmf,fmt)'        Dimensions="',zszV(3)," ",yszV(2)," ",xszV(1),'">'
    else if (output2D.eq.1) then
       write(ioxdmf,fmt)'        Dimensions="',zszV(3)," ",yszV(2)," ",1,'">'
    else if (output2D.eq.2) then
       write(ioxdmf,fmt)'        Dimensions="',zszV(3)," ",1," ",xszV(1),'">'
    else if (output2D.eq.3) then
       write(ioxdmf,fmt)'        Dimensions="',1," ",yszV(2)," ",xszV(1),'">'
    endif
    write(ioxdmf,'(A)')'    </Topology>'
  end subroutine write_xdmf_topo
  
  !
  ! Write the footer of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_footer()

    use decomp_2d, only : nrank

    implicit none

    if (nrank.eq.0) then
      write(ioxdmf,'(/)')
      write(ioxdmf,*)'    </Grid>'
      write(ioxdmf,*)'</Domain>'
      write(ioxdmf,'(A7)')'</Xdmf>'
      close(ioxdmf)
    endif

  end subroutine write_xdmf_footer

  !
  ! Write the given field for visualization
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !

  subroutine write_field(f1, pathname, filename, num, skip_ibm, flush)

    use mpi
    
    use var, only : ep1
    use var, only : zero, one
    use var, only : uvisu
    use param, only : iibm
	use param_p, only : isubdom
    use decomp_2d, only : mytype, xsize, xszV, yszV, zszV
    use decomp_2d, only : nrank, fine_to_coarseV
    use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_write_plane

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: f1
    character(len=*), intent(in) :: pathname, filename
    integer, intent(in) :: num
    logical, optional, intent(in) :: skip_ibm, flush

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: local_array
    
    logical :: mpiio, force_flush
    
    integer :: ierr

    integer :: precision
    character(len=:), allocatable :: fmt
    
#ifndef ADIOS2
    mpiio = .true.
#else
    mpiio = .false.
#endif

    if (present(flush)) then
       force_flush = flush
    else
       force_flush = .false.
    end if
 
    if (use_xdmf) then
       if (nrank.eq.0) then
          write(ioxdmf,*)'        <Attribute Name="'//filename//'" Center="Node">'
#ifndef ADIOS2
          write(ioxdmf,*)'           <DataItem Format="Binary"'
#else
          write(ioxdmf,*)'           <DataItem Format="HDF"'
#endif

#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
          if (output2D.eq.0) then
             precision = 4
          else
             precision = 8
          endif
#else
          precision = 8
#endif
#else
          precision = 4
#endif
          write(ioxdmf,"(A,I0,A)")'            DataType="Float" Precision="', precision, '" Endian="little" Seek="0"'

          fmt = "(A, I0, A, I0, A, I0, A)"
          if (output2D.eq.0) then
             write(ioxdmf,fmt)'            Dimensions="',zszV(3)," ",yszV(2)," ",xszV(1),'">'
          else if (output2D.eq.1) then
             write(ioxdmf,fmt)'            Dimensions="',zszV(3)," ",yszV(2)," ",1,'">'
          else if (output2D.eq.2) then
             write(ioxdmf,fmt)'            Dimensions="',zszV(3)," ",1," ",xszV(1),'">'
          else if (output2D.eq.3) then
             write(ioxdmf,fmt)'            Dimensions="',1," ",yszV(2)," ",xszV(1),'">'
          endif

          write(ioxdmf,*)'              '//gen_h5path(gen_filename(pathname, filename, num, 'bin'), num)

          write(ioxdmf,*)'           </DataItem>'
          write(ioxdmf,*)'        </Attribute>'
       endif
    endif

    if ((iibm == 2) .and. .not.present(skip_ibm)) then
       local_array(:,:,:) = (one - ep1(:,:,:)) * f1(:,:,:)
    else
       local_array(:,:,:) = f1(:,:,:)
    endif

    if (output2D.eq.0) then
       if (mpiio .or. (iibm == 2) .or. force_flush) then
          !! XXX: This (re)uses a temporary array for data - need to force synchronous writes.
          uvisu = zero
          
          call fine_to_coarseV(1,local_array,uvisu)
          call decomp_2d_write_one(1,uvisu,"data",gen_filename(pathname, filename, num, 'bin'),2,io_name,&
               opt_deferred_writes=.false.)
       else
          call decomp_2d_write_one(1,f1,"data",gen_filename(pathname, filename, num, 'bin'),0,io_name,opt_deferred_writes=.false.)
       end if
    else
       call decomp_2d_write_plane(1,local_array,output2D,-1,"data",gen_filename(pathname, filename, num, 'bin'),io_name)
    endif

  end subroutine write_field

  subroutine write_field_subdom(f1, pathname, filename, num, skip_ibm, flush)

    use mpi
    
    use var, only : ep1
    use var, only : zero, one
    use var, only : uvisu
    use param, only : iibm
	use param_p, only : isubdom, nx_first, nx_last, ny_first, ny_last, nz_first, nz_last
	
	use param_p, only : iplanez0, nx_first_z0, nx_last_z0, ny_first_z0, ny_last_z0, nz_first_z0, nz_last_z0
	use param_p, only : iplaney0, nx_first_y0, nx_last_y0, ny_first_y0, ny_last_y0, nz_first_y0, nz_last_y0
	
    use decomp_2d, only : mytype, xsize, xszV, yszV, zszV
    use decomp_2d, only : nrank, fine_to_coarseV
    use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_write_plane, decomp_2d_write_subdomain

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: f1
    character(len=*), intent(in) :: pathname, filename
    integer, intent(in) :: num
    logical, optional, intent(in) :: skip_ibm, flush

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: local_array
    
    logical :: mpiio, force_flush
    
    integer :: ierr

    integer :: precision
    character(len=:), allocatable :: fmt
    
#ifndef ADIOS2
    mpiio = .true.
#else
    mpiio = .false.
#endif

    if (present(flush)) then
       force_flush = flush
    else
       force_flush = .false.
    end if

    if ((iibm == 2) .and. .not.present(skip_ibm)) then
       local_array(:,:,:) = (one - ep1(:,:,:)) * f1(:,:,:)
    else
       local_array(:,:,:) = f1(:,:,:)
    endif

   if (mpiio .or. (iibm == 2) .or. (isubdom == 1) .or. force_flush) then
      !! XXX: This (re)uses a temporary array for data - need to force synchronous writes.
      uvisu = zero
          
      call fine_to_coarseV(1,local_array,uvisu)
      !call decomp_2d_write_one(1,uvisu,"data_subdom",gen_filename(pathname, filename, num, 'bin'),2,io_name,&
      !     opt_deferred_writes=.false.)
	  call decomp_2d_write_subdomain(1,uvisu,nx_first,nx_last,ny_first,ny_last,nz_first,nz_last,gen_filename(pathname, filename, num, 'bin'))

   end if
   
   if (mpiio .or. (iibm == 2) .or. (iplanez0 == 1)  .or. (iplaney0 == 1) .or. force_flush) then
      !! XXX: This (re)uses a temporary array for data - need to force synchronous writes.
      uvisu = zero
          
      call fine_to_coarseV(1,local_array,uvisu)
      !call decomp_2d_write_one(1,uvisu,"data_subdom",gen_filename(pathname, filename, num, 'bin'),2,io_name,&
      !     opt_deferred_writes=.false.)
	  if (iplanez0 == 1) then
		call decomp_2d_write_subdomain(1,uvisu,nx_first_z0,nx_last_z0,ny_first_z0,ny_last_z0,nz_first_z0,nz_last_z0,gen_filename(pathname, filename, num, 'bin'))
      endif
	  if (iplaney0 == 1) then
		call decomp_2d_write_subdomain(1,uvisu,nx_first_y0,nx_last_y0,ny_first_y0,ny_last_y0,nz_first_y0,nz_last_y0,gen_filename(pathname, filename, num, 'bin'))
      endif
   end if

  end subroutine write_field_subdom


subroutine merge_and_cleanup_subdom_files(pathname, num, nx_td, ny_td, nz_td)
  use mpi
  use var, only : mytype

  implicit none

  character(len=*), intent(in) :: pathname
  integer, intent(in) :: num, nx_td, ny_td, nz_td

  character(len=:), allocatable :: fname_ux, fname_uy, fname_uz, fname_pp, fname_out
  character(len=512) :: tmp
  integer :: unit_in, unit_out
  integer :: iostat_val
  logical :: exist_flag

  integer :: N
  integer(kind=8) :: fb
  real(mytype), allocatable :: buffer(:)

  ! Total number of elements per field
  N = nx_td * ny_td * nz_td

  ! Byte size of one field
  fb = int(N,8) * storage_size(1.0_mytype)/8_8

  ! Allocate buffer for one field
  allocate(buffer(N))

  ! Build filenames (adjust pattern if your gen_filename differs)
  write(tmp,'(A,"/ux_subdom-",I0,".bin")') trim(pathname), num
  fname_ux = trim(tmp)

  write(tmp,'(A,"/uy_subdom-",I0,".bin")') trim(pathname), num
  fname_uy = trim(tmp)

  write(tmp,'(A,"/uz_subdom-",I0,".bin")') trim(pathname), num
  fname_uz = trim(tmp)

  write(tmp,'(A,"/pp_subdom-",I0,".bin")') trim(pathname), num
  fname_pp = trim(tmp)

  write(tmp,'(A,"/",I0,".bin")') trim(pathname), num
  fname_out = trim(tmp)

  ! Check that all four files exist
  inquire(file=fname_ux, exist=exist_flag)
  if (.not. exist_flag) then
     print *, 'merge_and_cleanup_subdom_files: missing ', fname_ux
     deallocate(buffer)
     return
  end if
  inquire(file=fname_uy, exist=exist_flag)
  if (.not. exist_flag) then
     print *, 'merge_and_cleanup_subdom_files: missing ', fname_uy
     deallocate(buffer)
     return
  end if
  inquire(file=fname_uz, exist=exist_flag)
  if (.not. exist_flag) then
     print *, 'merge_and_cleanup_subdom_files: missing ', fname_uz
     deallocate(buffer)
     return
  end if
  inquire(file=fname_pp, exist=exist_flag)
  if (.not. exist_flag) then
     print *, 'merge_and_cleanup_subdom_files: missing ', fname_pp
     deallocate(buffer)
     return
  end if

  ! Open output file
  open(newunit=unit_out, file=fname_out, form='unformatted', &
       access='stream', status='replace', action='write', iostat=iostat_val)
  if (iostat_val /= 0) then
     print *, 'merge_and_cleanup_subdom_files: cannot create ', fname_out
     deallocate(buffer)
     return
  end if

  ! --- ux ---
  open(newunit=unit_in, file=fname_ux, form='unformatted', &
       access='stream', status='old', action='read', iostat=iostat_val)
  if (iostat_val /= 0) stop 'cannot open ux file'
  read(unit_in) buffer
  close(unit_in)
  write(unit_out, pos=1_8) buffer

  ! --- uy ---
  open(newunit=unit_in, file=fname_uy, form='unformatted', &
       access='stream', status='old', action='read', iostat=iostat_val)
  if (iostat_val /= 0) stop 'cannot open uy file'
  read(unit_in) buffer
  close(unit_in)
  write(unit_out, pos=1_8 + fb) buffer

  ! --- uz ---
  open(newunit=unit_in, file=fname_uz, form='unformatted', &
       access='stream', status='old', action='read', iostat=iostat_val)
  if (iostat_val /= 0) stop 'cannot open uz file'
  read(unit_in) buffer
  close(unit_in)
  write(unit_out, pos=1_8 + 2*fb) buffer

  ! --- pp ---
  open(newunit=unit_in, file=fname_pp, form='unformatted', &
       access='stream', status='old', action='read', iostat=iostat_val)
  if (iostat_val /= 0) stop 'cannot open pp file'
  read(unit_in) buffer
  close(unit_in)
  write(unit_out, pos=1_8 + 3*fb) buffer

  close(unit_out)

  deallocate(buffer)

  ! Delete the four intermediate files
  call system('rm -f ' // fname_ux)
  call system('rm -f ' // fname_uy)
  call system('rm -f ' // fname_uz)
  call system('rm -f ' // fname_pp)

end subroutine merge_and_cleanup_subdom_files




  function gen_snapshotname(pathname, varname, num, ext)
    character(len=*), intent(in) :: pathname, varname, ext
    integer, intent(in) :: num
    character(len=:), allocatable :: gen_snapshotname
#ifndef ADIOS2
    gen_snapshotname = gen_filename(pathname, varname, num, ext)
#else
    gen_snapshotname = varname//'-'//int_to_str(num)//'.'//ext
#endif
  end function gen_snapshotname
  
  function gen_filename(pathname, varname, num, ext)

    character(len=*), intent(in) :: pathname, varname, ext
    integer, intent(in) :: num
#ifndef ADIOS2
    character(len=:), allocatable :: gen_filename
    gen_filename = pathname//'/'//varname//'-'//int_to_str(num)//'.'//ext
#else
    character(len=len(varname)) :: gen_filename
    write(gen_filename, "(A)") varname
#endif
    
  end function gen_filename

  function gen_h5path(filename, num)

    character(len=*), intent(in) :: filename
    integer, intent(in) :: num
#ifndef ADIOS2
    character(len=*), parameter :: path_to_h5file = "./"
    character(len=(len(path_to_h5file) + len(filename))) :: gen_h5path
    write(gen_h5path, "(A)") path_to_h5file//filename
#else
    character(len=*), parameter :: path_to_h5file = "../data.hdf5:/Step"
    character(len=:), allocatable :: gen_h5path
    gen_h5path = path_to_h5file//int_to_str(num)//"/"//filename
#endif
    
  end function gen_h5path

  ! Function converting an integer to a string.
  function int_to_str(i)

    integer, intent(in) :: i ! Integer input.

    ! String return value. The string must be long enough to contain all the characters required to
    ! represent the integer, i.e 1 + log_10(i). To protect against calling with integer 0 the value
    ! passed to log_10 must be >= 1.
    character(len=:), allocatable :: int_to_str
    character(len=32) :: buf

    write(buf, "(I0)") i
    int_to_str = trim(buf)
  end function int_to_str
  
end module visu
