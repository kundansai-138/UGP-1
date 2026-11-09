! =====================================================================
! spod-adios2-nodelocal.F90
!
! Streaming SPOD reading Xcompact3d data through ADIOS2 (SST), with the
! node-local input layout: each SPOD rank reads only data written by
! Xcompact3d ranks on its own node. The SPOD maths is the same as in
! spod-adios2.txt.
!
! One file, three program units, in this order:
!   module node_local_map       which boxes each SPOD rank reads
!   module streamingSPODModule  SPOD (PETSc/SLEPc) + ADIOS2 input
!   program streamingSPOD
!
! Build: exactly like your working ADIOS2 version, with this file as the
! source, e.g. in the Makefile:   SRC = spod-adios2-nodelocal.o
! Needs PETSc/SLEPc (complex scalars), FFTW with MPI, ADIOS2 (Fortran).
!
! Before running, edit Subroutine UserInput:
!   x3dProw, x3dPcol : p_row, p_col of the Xcompact3d run. Set them
!                      explicitly in its input.i3d (not 0) and copy them
!                      here. Left at 0, the run stops at start-up.
!   gridNx/Ny/Nz, streamName, adiosConfigFile, numberOfSnapShot, nDFT,
!   k, percentageOfOverall, weightMatrixFilePath, nFreqToSave
! Xcompact3d must write full 3-D fields (output2D = 0, nvisu = 1).
! Launch both codes on the same nodes, ranks filled node by node, e.g.
!   mpiexec -f hosts -ppn 32 -n 128 ./xcompact3d &
!   mpiexec -f hosts -ppn 16 -n 64  ./streamingSPOD
! =====================================================================

! ======================================================================
! node_local_map
!
! Decides which boxes of the Xcompact3d (X3D) output each SPOD rank reads,
! so that a SPOD rank only reads data written by X3D ranks on its own node.
!
! Plain Fortran: no MPI, PETSc or ADIOS2, so it can be tested on its own
! (see test_node_local_map.f90). Each routine mirrors the function with the
! same name in tools/node_local_model.py.
!
! Background (2decomp-fft v2.0.4, used by Xcompact3d):
!   * X3D rank b sits at (coord1, coord2) = (b / p_col, mod(b, p_col)) of a
!     p_row x p_col grid (MPI_Cart_create, no reordering).
!   * Pencil 1 (X, what Xcompact3d visu writes): x is whole, y is split over
!     p_row by coord1, z is split over p_col by coord2.
!   * Each rank writes its pencil as one block of the 3-D global array,
!     start = xstart-1, count = xsize.
! ======================================================================
module node_local_map
  implicit none
  private

  integer, parameter, public :: i8 = 8     ! same kind as ADIOS2's integer(kind=8) start/count

  ! One box to read: 0-based global start / count in Fortran order (x,y,z)
  ! and the X3D rank whose pencil contains it.
  type, public :: box_t
    integer(i8) :: start(3) = 0
    integer(i8) :: count(3) = 0
    integer     :: writer   = -1
  end type box_t

  public :: distribute_2decomp, writer_box, writer_node, split_range, &
            line_range_to_boxes, node_local_boxes, box_points

contains

  ! 1-D split of n points over p ranks, as 2decomp does it.
  ! extraFirst = .false.: the n mod p extra points go to the LAST ranks
  ! (always the case in v2.0.4, and the default in newer versions).
  subroutine distribute_2decomp(n, p, extraFirst, st, sz)
    integer(i8), intent(in)  :: n
    integer,     intent(in)  :: p
    logical,     intent(in)  :: extraFirst
    integer(i8), intent(out) :: st(0:p-1), sz(0:p-1)
    integer :: i, r

    sz = n / p
    r  = int(mod(n, int(p, i8)))
    if (extraFirst) then
      do i = 0, r - 1
        sz(i) = sz(i) + 1
      end do
    else
      do i = p - r, p - 1
        sz(i) = sz(i) + 1
      end do
    end if
    st(0) = 0
    do i = 1, p - 1
      st(i) = st(i-1) + sz(i-1)
    end do
  end subroutine distribute_2decomp


  ! Box written by X3D rank b (0-based start, count).
  subroutine writer_box(b, nx, ny, nz, prow, pcol, pencil, extraFirst, start, count)
    integer,     intent(in)  :: b, prow, pcol, pencil
    integer(i8), intent(in)  :: nx, ny, nz
    logical,     intent(in)  :: extraFirst
    integer(i8), intent(out) :: start(3), count(3)
    integer(i8) :: st1(0:prow-1), sz1(0:prow-1), st2(0:pcol-1), sz2(0:pcol-1)
    integer :: c1, c2

    c1 = b / pcol
    c2 = mod(b, pcol)
    select case (pencil)
    case (1)          ! x whole, y over p_row, z over p_col
      call distribute_2decomp(ny, prow, extraFirst, st1, sz1)
      call distribute_2decomp(nz, pcol, extraFirst, st2, sz2)
      start = [0_i8,    st1(c1), st2(c2)]
      count = [nx,      sz1(c1), sz2(c2)]
    case (2)          ! x over p_row, y whole, z over p_col
      call distribute_2decomp(nx, prow, extraFirst, st1, sz1)
      call distribute_2decomp(nz, pcol, extraFirst, st2, sz2)
      start = [st1(c1), 0_i8,    st2(c2)]
      count = [sz1(c1), ny,      sz2(c2)]
    case (3)          ! x over p_row, y over p_col, z whole
      call distribute_2decomp(nx, prow, extraFirst, st1, sz1)
      call distribute_2decomp(ny, pcol, extraFirst, st2, sz2)
      start = [st1(c1), st2(c2), 0_i8]
      count = [sz1(c1), sz2(c2), nz]
    case default
      start = -1
      count = 0
    end select
  end subroutine writer_box


  ! Node index of X3D rank b.
  ! perNode > 0 : ranks 0..perNode-1 on node 0, the next perNode on node 1, ...
  ! perNode = 0 : sequential (block) placement of nWriters ranks on nNodes
  !               nodes, the first mod(nWriters,nNodes) nodes get one more
  !               (what Slurm/mpirun do by default).
  ! A result >= nNodes is wrapped, so every writer is read by some SPOD node.
  integer function writer_node(b, nWriters, nNodes, perNode)
    integer, intent(in) :: b, nWriters, nNodes, perNode
    integer :: base, extra

    if (perNode > 0) then
      writer_node = mod(b / perNode, nNodes)
      return
    end if
    base  = nWriters / nNodes
    extra = mod(nWriters, nNodes)
    if (b < extra * (base + 1)) then
      writer_node = b / (base + 1)
    else
      writer_node = extra + (b - extra * (base + 1)) / base
    end if
  end function writer_node


  ! Balanced contiguous split of total items over nParts: part gets [lo, hi).
  ! Same rule as PETSC_DECIDE.
  subroutine split_range(total, nParts, part, lo, hi)
    integer(i8), intent(in)  :: total
    integer,     intent(in)  :: nParts, part
    integer(i8), intent(out) :: lo, hi
    integer(i8) :: base, extra

    base  = total / nParts
    extra = mod(total, int(nParts, i8))
    lo = part * base + min(int(part, i8), extra)
    hi = lo + base
    if (part < extra) hi = hi + 1
  end subroutine split_range


  ! Lines [a, e) of a writer box -> at most 3 boxes, appended to boxes(nb+1:).
  ! A line is one run along x (wcount(1) points). Line l of the box sits at
  ! y = mod(l, wcount(2)), z = l / wcount(2). Boxes come out in line order:
  ! partial first plane, whole planes, partial last plane.
  subroutine line_range_to_boxes(wstart, wcount, a, e, writer, boxes, nb)
    integer(i8), intent(in)    :: wstart(3), wcount(3), a, e
    integer,     intent(in)    :: writer
    type(box_t), intent(inout) :: boxes(:)
    integer,     intent(inout) :: nb
    integer(i8) :: c2, za, ya, zb, yb
    logical :: haveLast

    c2 = wcount(2)
    za = a / c2
    ya = mod(a, c2)
    zb = (e - 1) / c2
    yb = mod(e - 1, c2)

    if (za == zb) then
      call add(ya, za, yb - ya + 1, 1_i8)
    else
      if (ya > 0) then                       ! partial first plane
        call add(ya, za, c2 - ya, 1_i8)
        za = za + 1
      end if
      haveLast = (yb < c2 - 1)
      if (haveLast) zb = zb - 1
      if (za <= zb) call add(0_i8, za, c2, zb - za + 1)   ! whole planes
      if (haveLast) call add(0_i8, zb + 1, yb + 1, 1_i8)  ! partial last plane
    end if

  contains

    subroutine add(y0, z0, cy, cz)
      integer(i8), intent(in) :: y0, z0, cy, cz
      nb = nb + 1
      boxes(nb)%start  = wstart + [0_i8, y0, z0]
      boxes(nb)%count  = [wcount(1), cy, cz]
      boxes(nb)%writer = writer
    end subroutine add

  end subroutine line_range_to_boxes


  ! Boxes read by the SPOD rank with index nodeRank (0..nodeSize-1) on node
  ! myNode: all lines written by the X3D ranks of that node, concatenated in
  ! X3D rank order and split evenly over the nodeSize SPOD ranks of the node.
  ! Every returned box lies inside one X3D pencil (boxes(i)%writer).
  subroutine node_local_boxes(myNode, nodeRank, nodeSize, nNodes, nx, ny, nz, &
                              prow, pcol, pencil, extraFirst, perNode, boxes, nb)
    integer,     intent(in)  :: myNode, nodeRank, nodeSize, nNodes
    integer(i8), intent(in)  :: nx, ny, nz
    integer,     intent(in)  :: prow, pcol, pencil, perNode
    logical,     intent(in)  :: extraFirst
    type(box_t), allocatable, intent(out) :: boxes(:)
    integer,     intent(out) :: nb
    integer :: b, nW, nMine
    integer(i8) :: total, lo, hi, off, nl, a, e, ws(3), wc(3)

    nW = prow * pcol

    ! pass 1: how many lines do the writers on my node hold?
    total = 0
    nMine = 0
    do b = 0, nW - 1
      if (writer_node(b, nW, nNodes, perNode) /= myNode) cycle
      call writer_box(b, nx, ny, nz, prow, pcol, pencil, extraFirst, ws, wc)
      total = total + wc(2) * wc(3)
      nMine = nMine + 1
    end do

    allocate(boxes(max(1, 3 * nMine)))
    nb = 0
    if (nMine == 0) return
    call split_range(total, nodeSize, nodeRank, lo, hi)

    ! pass 2: cut my share [lo, hi) out of the concatenated line list
    off = 0
    do b = 0, nW - 1
      if (writer_node(b, nW, nNodes, perNode) /= myNode) cycle
      call writer_box(b, nx, ny, nz, prow, pcol, pencil, extraFirst, ws, wc)
      nl = wc(2) * wc(3)
      a  = max(lo, off)
      e  = min(hi, off + nl)
      if (a < e) call line_range_to_boxes(ws, wc, a - off, e - off, b, boxes, nb)
      off = off + nl
    end do
  end subroutine node_local_boxes


  integer(i8) function box_points(bx)
    type(box_t), intent(in) :: bx
    box_points = bx%count(1) * bx%count(2) * bx%count(3)
  end function box_points

end module node_local_map


Module streamingSPODModule
#include <slepc/finclude/slepcsvd.h>

                              use, intrinsic :: iso_c_binding
                              use slepc
                              use slepcsvd
                              use adios2
                              use node_local_map, only: box_t, node_local_boxes, writer_node, box_points

                              implicit none

                            ! ================================================================
                            ! ADIOS2 SST INPUT, NODE-LOCAL LAYOUT
                            !
                            ! Snapshot vector q = [ux ; uy ; uz ; pp]. In the NATURAL order each
                            ! block is one Xcompact3d field flattened in Fortran order
                            ! p = i + nx*(j + ny*kz) (i fastest, 0-based).
                            !
                            ! The SPOD maths only does row-by-row updates plus sums over all rows
                            ! (U^H x, norms), so it gives the same result for ANY ordering of the
                            ! rows, as long as data and weights use the same ordering.
                            ! This version therefore stores the rows in a NODE-LOCAL order:
                            !   * each SPOD rank reads only boxes that lie inside the pencil of one
                            !     Xcompact3d rank running on the SAME node (no inter-node traffic
                            !     for snapshot data),
                            !   * its PETSc rows are exactly those points (local size nLocalRows,
                            !     no longer PETSC_DECIDE),
                            !   * weights (natural order on disk) and saved modes (natural order on
                            !     disk) are moved between the two orders with one VecScatter.
                            ! ================================================================
                            integer, parameter :: nStreamVars = 4
                            character(len=8)   :: streamVarNames(nStreamVars)
                            character(len=500) :: adiosConfigFile   ! XML with IO 'spod-input'
                            character(len=500) :: streamName        ! SST stream (contact file <name>.sst)
                            integer(kind=8)    :: gridNx, gridNy, gridNz
                            real               :: beginStepTimeout  ! seconds, -1 = wait forever

                            ! Xcompact3d decomposition of the run that writes the stream
                            integer :: x3dProw, x3dPcol             ! p_row, p_col from input.i3d (explicit, not 0)
                            integer :: x3dPencil                    ! pencil of the written fields (visu: 1 = X)
                            logical :: x3dExtraPointsFirst          ! 2decomp split rule (.false. for v2.0.4)
                            integer :: x3dRanksPerNode              ! 0 = sequential over the SPOD nodes

                            type(adios2_adios)  :: adios
                            type(adios2_io)     :: ioIn
                            type(adios2_engine) :: engineIn

                            ! Where this rank sits and what it reads
                            integer, parameter :: hostLen = 256
                            character(len=hostLen) :: myHost
                            integer :: myNode, nodeRank, nodeSize, nNodes
                            type(box_t), allocatable :: myBoxes(:)  ! boxes read every step (same for all vars)
                            integer :: nMyBoxes
                            integer(kind=8) :: nMyPoints            ! grid points per variable on this rank
                            PetscInt :: nLocalRows                  ! = nStreamVars * nMyPoints

                            type :: rbuf_t
                              real(kind=8), allocatable :: v(:)
                            end type rbuf_t
                            type(rbuf_t), allocatable :: getBuf(:,:) ! (box, variable) ADIOS2 targets

                            ! Natural order <-> node-local order (weights in, modes out)
                            Vec :: natVec                           ! natural order, PETSC_DECIDE layout
                            VecScatter :: natToLocal                ! forward: natural -> node-local

                            double precision :: tWaitStep = 0.0d0, tRead = 0.0d0, tCompute = 0.0d0


                            ! ================================================================
                            ! USER-DEFINED CONSTANTS
                            ! ================================================================

                            integer :: numberOfDataInSnapShot
                            integer :: numberOfSnapShot
                            integer :: nDFT
                            integer :: k
                            real :: percentageOfOverall
                            logical :: useHamming
                            integer :: nFreqToSave

                            ! User can change these values as needed

                            character(len=500) :: filenameFormat
                            character(len=500) :: weightMatrixFilePath

                            ! User can change this path and pattern as needed


                            ! ================================================================
                            ! VARIABLES SET INTERNALLY
                            ! ================================================================

                            integer :: nFreq
                            integer :: nOvlp                 ! Number of overlapping points
                            integer :: currentSnapshotNumber ! Current snapshot number
                            integer :: fileIndex             ! File index for output
                            character(len=200) :: filename   ! Output filename

                            ! ================================================================
                            ! VARIABLES FOR WINDOW
                            ! ================================================================
                            real :: winWeight
                            PetscScalar, allocatable :: window(:)


                            ! ================================================================
                            !  MPI VARIABLES AND ERROR VARIABLES
                            ! ================================================================
                            ! Declare PETSc and MPI-related variables
                            PetscMPIInt    :: rank, size   ! rank: process rank, size: total number of processes
                            PetscErrorCode :: ierr         ! ierr: error code returned by PETSc/MPI routines

                            ! ================================================================
                            !  VARIABLES USED TO DECIDE WHICH BLOCK TO UPDATE
                            ! ================================================================
                            integer, allocatable :: t_idx(:)
                            integer :: dn

                            ! ================================================================
                            !  VARIABLES USED TO DFT
                            ! ================================================================
                            PetscScalar, allocatable :: fourierF90Array(:,:)

                            ! ================================================================
                            !  VARIABLES USED TO MEAN UPDATE
                            ! ================================================================
                            integer :: t_i

                            ! ================================================================
                            !  VARIABLES USED TO CONSTRUCT K MATRIX
                            ! ================================================================
                            integer :: current_Block

                            ! ================================================================
                            !  PETSC MATRIX
                            ! ================================================================
                            Mat :: X_sum_BlockA
                            Mat :: X_sum_BlockB
                            Mat :: mu
                            Mat :: x_new
                            Mat :: U_hat
                            Mat :: S_hat
                            Mat :: FourierRange ! Only used for range

                            ! ================================================================
                            !  PETSC VECTORS
                            ! ================================================================
                            Vec :: sqrtW


                              contains
                              Subroutine UserInput
                                implicit none
                                ! =======================
                                ! USER-DEFINED VARIABLES
                                ! =======================

                                ! =======================
                                ! ADIOS2 SST INPUT FROM XCOMPACT3D
                                ! =======================
                                ! Global grid of each Xcompact3d field, Fortran order (nx, ny, nz).
                                ! Checked against the shape announced in the stream at every step.
                                gridNx = 301
                                gridNy = 301
                                gridNz = 301

                                ! Variables stacked into one snapshot vector, in this order
                                streamVarNames = [character(len=8) :: 'ux', 'uy', 'uz', 'pp']

                                ! Stream name used by Xcompact3d 'solution-io' (contact file data.sst).
                                ! Relative name => run SPOD in the Xcompact3d run folder,
                                ! or put the absolute path of the run folder in front of 'data'.
                                streamName = 'data'
                                adiosConfigFile = 'spod_sst.xml'
                                beginStepTimeout = -1.0     ! wait forever for the next step

                                ! =======================
                                ! XCOMPACT3D DECOMPOSITION (needed for the node-local layout)
                                ! =======================
                                ! p_row and p_col of the Xcompact3d run. Set them explicitly in its
                                ! input.i3d (not 0 = automatic) and copy the same two numbers here.
                                ! 0 stops the run at start-up. Wrong non-zero values do NOT change the
                                ! results, they only make reads remote.
                                x3dProw = 0
                                x3dPcol = 0
                                ! Pencil of the written fields. Xcompact3d visu writes X-pencils.
                                x3dPencil = 1
                                ! 2decomp gives the extra points of an uneven split to the last ranks
                                ! (v2.0.4, used by Xcompact3d). Only change if your 2decomp differs.
                                x3dExtraPointsFirst = .false.
                                ! 0: Xcompact3d ranks are placed sequentially (block) on the same nodes,
                                !    in the same node order, as the SPOD ranks.
                                ! >0: Xcompact3d ranks 0..n-1 on node 0, n..2n-1 on node 1, ...
                                x3dRanksPerNode = 0

                                ! Number of data points in each snapshot = 4 variables x grid
                                numberOfDataInSnapShot = int(nStreamVars * gridNx * gridNy * gridNz)

                                ! Maximum number of snapshots to consume (loop also stops at end of stream)
                                numberOfSnapShot = 31

                                ! Number of data points for DFT window
                                nDFT = 8

                                ! K modes
                                k = 3

                                ! Overlap percentage between DFT windows (0.0 to 1.0)
                                percentageOfOverall = 0.5
                                ! User can change this value as needed

                                ! Snapshots now come from the ADIOS2 SST stream, not from files
                                filenameFormat = '(no files: ADIOS2 SST input)'

                                ! Weight matrix filepath, NATURAL order (same file as before)
                                ! NOTE: IF PATH IS SET AS EMPTY i,e '' THEN 1 WILL BE SET AS WEIGHT i,e EQUAL WEIGHTS
                                weightMatrixFilePath = ''   ! equal weights for now

                                ! Hamming window or not
                                useHamming = .false.   ! User can set this to .false. if desired

                                ! Number of frequencies to save.
                                ! This will be used to save the nFreqToSave frequencies corresponding to the nFreqToSave largest eigenvalues.
                                nFreqToSave = 5


                              End Subroutine UserInput


                                ! ================================================================
                                ! NODE-LOCAL LAYOUT
                                ! ================================================================
                                ! Find this rank's node, then the boxes it reads (node_local_map).
                                ! Sets myBoxes, nMyBoxes, nMyPoints, nLocalRows.
                                subroutine SetupNodeLocalLayout()
                                  implicit none
                                  integer :: nodeComm, leaderComm, color, nameLen, ib, b, m, r, nW, mpiErr
                                  integer :: wFirst, wLast, wCount, sCount
                                  integer(kind=8) :: myInfo(3), sumRows, rowsMin, rowsMax, nodeRows
                                  integer(kind=8), allocatable :: allInfo(:,:)
                                  character(len=hostLen), allocatable :: allHosts(:)
                                  character(len=MPI_MAX_PROCESSOR_NAME) :: hostName
                                  character(len=300) :: msg

                                  ! ---- input checks ----
                                  if (x3dProw <= 0 .or. x3dPcol <= 0) call stream_abort( &
                                    'set x3dProw, x3dPcol in UserInput to p_row, p_col of the Xcompact3d run')
                                  if (x3dPencil < 1 .or. x3dPencil > 3) call stream_abort('x3dPencil must be 1, 2 or 3')
                                  if (x3dProw > min(gridNx, gridNy) .or. x3dPcol > min(gridNy, gridNz)) &
                                    call stream_abort('p_row/p_col too large for the grid')
                                  if (int(nStreamVars, 8) * gridNx * gridNy * gridNz /= int(numberOfDataInSnapShot, 8)) &
                                    call stream_abort('snapshot too large for a default integer (needs 64-bit indices)')
                                  nW = x3dProw * x3dPcol

                                  ! ---- which node am I on? ----
                                  ! nodeComm: ranks sharing my node, ordered by world rank.
                                  ! Node index: order of the lowest rank of each node (leaders comm).
                                  call MPI_Comm_split_type(PETSC_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, &
                                                           MPI_INFO_NULL, nodeComm, mpiErr)
                                  call MPI_Comm_rank(nodeComm, nodeRank, mpiErr)
                                  call MPI_Comm_size(nodeComm, nodeSize, mpiErr)
                                  color = MPI_UNDEFINED
                                  if (nodeRank == 0) color = 0
                                  call MPI_Comm_split(PETSC_COMM_WORLD, color, rank, leaderComm, mpiErr)
                                  if (nodeRank == 0) then
                                    call MPI_Comm_rank(leaderComm, myNode, mpiErr)
                                    call MPI_Comm_size(leaderComm, nNodes, mpiErr)
                                    call MPI_Comm_free(leaderComm, mpiErr)
                                  end if
                                  call MPI_Bcast(myNode, 1, MPI_INTEGER, 0, nodeComm, mpiErr)
                                  call MPI_Bcast(nNodes, 1, MPI_INTEGER, 0, nodeComm, mpiErr)
                                  call MPI_Comm_free(nodeComm, mpiErr)
                                  hostName = ' '
                                  call MPI_Get_processor_name(hostName, nameLen, mpiErr)
                                  myHost = hostName(1:min(nameLen, hostLen))

                                  ! ---- my boxes: lines of the X3D ranks on my node, split over the
                                  !      SPOD ranks of my node ----
                                  call node_local_boxes(myNode, nodeRank, nodeSize, nNodes, gridNx, gridNy, gridNz, &
                                                        x3dProw, x3dPcol, x3dPencil, x3dExtraPointsFirst, &
                                                        x3dRanksPerNode, myBoxes, nMyBoxes)
                                  nMyPoints = 0
                                  do ib = 1, nMyBoxes
                                    nMyPoints = nMyPoints + box_points(myBoxes(ib))
                                  end do
                                  nLocalRows = int(nStreamVars * nMyPoints, kind=kind(nLocalRows))

                                  ! ---- all rows owned exactly once in total? ----
                                  call MPI_Allreduce(int(nLocalRows, 8), sumRows, 1, MPI_INTEGER8, MPI_SUM, &
                                                     PETSC_COMM_WORLD, mpiErr)
                                  if (sumRows /= int(numberOfDataInSnapShot, 8)) then
                                    write(msg, '(A,I0,A,I0)') 'node-local rows add up to ', sumRows, &
                                      ' instead of ', numberOfDataInSnapShot
                                    call stream_abort(msg)
                                  end if

                                  ! ---- report: one line per node on rank 0 ----
                                  myInfo = [int(myNode, 8), int(nodeRank, 8), int(nLocalRows, 8)]
                                  if (rank == 0) then
                                    allocate(allInfo(3, 0:size-1), allHosts(0:size-1))
                                  else
                                    allocate(allInfo(3, 1), allHosts(1))
                                  end if
                                  call MPI_Gather(myInfo, 3, MPI_INTEGER8, allInfo, 3, MPI_INTEGER8, 0, &
                                                  PETSC_COMM_WORLD, mpiErr)
                                  call MPI_Gather(myHost, hostLen, MPI_CHARACTER, allHosts, hostLen, MPI_CHARACTER, 0, &
                                                  PETSC_COMM_WORLD, mpiErr)
                                  if (rank == 0) then
                                    print *, '-------------------------------------------'
                                    print '(A,I0,A,I0,A,I0,A,I0)', ' Node-local input: ', nNodes, ' nodes, ', size, &
                                      ' SPOD ranks, Xcompact3d ', x3dProw, ' x ', x3dPcol
                                    do m = 0, nNodes - 1
                                      sCount = 0
                                      nodeRows = 0
                                      do r = size - 1, 0, -1
                                        if (allInfo(1, r) == m) then
                                          sCount = sCount + 1
                                          nodeRows = nodeRows + allInfo(3, r)
                                          ib = r                      ! ends on the lowest rank of node m
                                        end if
                                      end do
                                      wFirst = -1
                                      wLast = -1
                                      wCount = 0
                                      do b = 0, nW - 1
                                        if (writer_node(b, nW, nNodes, x3dRanksPerNode) == m) then
                                          if (wFirst < 0) wFirst = b
                                          wLast = b
                                          wCount = wCount + 1
                                        end if
                                      end do
                                      print '(A,I4,A,A,A,I5,A,I6,A,I6,A,I5,A,I12)', '  node', m, '  ', &
                                        trim(allHosts(ib)), '  SPOD ranks', sCount, '  reads X3D ranks', wFirst, &
                                        ' ..', wLast, ' (', wCount, ')  rows', nodeRows
                                    end do
                                    rowsMin = minval(allInfo(3, :))
                                    rowsMax = maxval(allInfo(3, :))
                                    print '(A,I0,A,I0)', '  rows per SPOD rank: min ', rowsMin, '  max ', rowsMax
                                    if (rowsMin == 0) print *, &
                                      '  WARNING: some SPOD ranks own no rows (more SPOD ranks than x-lines)'
                                    print *, '  Check: Xcompact3d ranks listed for a node must run on that host.'
                                    print *, '-------------------------------------------'
                                  end if
                                  deallocate(allInfo, allHosts)
                                end subroutine SetupNodeLocalLayout


                                ! Natural index of every local row, then the scatter between the
                                ! natural-order vector (PETSC_DECIDE) and the node-local layout.
                                ! Also checks that the rows form an exact permutation.
                                subroutine BuildNaturalScatter()
                                  implicit none
                                  Vec :: locVec
                                  IS :: isNat, isLoc
                                  PetscInt, allocatable :: natIdx(:)
                                  PetscInt :: rowStart, rowEnd, one
                                  PetscReal :: nrm1, nrmInf
                                  integer :: v, ib
                                  integer(kind=8) :: npts, r, i, j, kk, g0

                                  ! local rows: variable-major, then box, then z, y, x (same order as Load_Snapshot)
                                  npts = gridNx * gridNy * gridNz
                                  allocate(natIdx(max(int(nLocalRows, 8), 1_8)))
                                  r = 0
                                  do v = 1, nStreamVars
                                    do ib = 1, nMyBoxes
                                      associate (s => myBoxes(ib)%start, c => myBoxes(ib)%count)
                                        do kk = 0, c(3) - 1
                                          do j = 0, c(2) - 1
                                            g0 = (v - 1) * npts + s(1) + gridNx * ((s(2) + j) + gridNy * (s(3) + kk))
                                            do i = 0, c(1) - 1
                                              r = r + 1
                                              natIdx(r) = int(g0 + i, kind=kind(natIdx))
                                            end do
                                          end do
                                        end do
                                      end associate
                                    end do
                                  end do

                                  one = 1
                                  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numberOfDataInSnapShot, natVec, ierr)
                                  call VecCreateMPI(PETSC_COMM_WORLD, nLocalRows, numberOfDataInSnapShot, locVec, ierr)
                                  call VecGetOwnershipRange(locVec, rowStart, rowEnd, ierr)
                                  call ISCreateGeneral(PETSC_COMM_WORLD, nLocalRows, natIdx, PETSC_COPY_VALUES, isNat, ierr)
                                  call ISCreateStride(PETSC_COMM_WORLD, nLocalRows, rowStart, one, isLoc, ierr)
                                  ! natVec(natIdx(i)) -> locVec(rowStart + i - 1)
                                  call VecScatterCreate(natVec, isNat, locVec, isLoc, natToLocal, ierr)
                                  call ISDestroy(isNat, ierr)
                                  call ISDestroy(isLoc, ierr)
                                  deallocate(natIdx)

                                  ! Self-check: send 1 from every local row back to its natural row and
                                  ! add. Every natural row must receive exactly 1.
                                  call VecSet(locVec, (1.0d0, 0.0d0), ierr)
                                  call VecSet(natVec, (0.0d0, 0.0d0), ierr)
                                  call VecScatterBegin(natToLocal, locVec, natVec, ADD_VALUES, SCATTER_REVERSE, ierr)
                                  call VecScatterEnd(natToLocal, locVec, natVec, ADD_VALUES, SCATTER_REVERSE, ierr)
                                  call VecNorm(natVec, NORM_1, nrm1, ierr)
                                  call VecNorm(natVec, NORM_INFINITY, nrmInf, ierr)
                                  if (abs(nrm1 - real(numberOfDataInSnapShot, kind(nrm1))) > 0.5d0 .or. &
                                      abs(nrmInf - 1.0d0) > 1.0d-12) then
                                    call stream_abort('node-local rows are not a permutation of the natural rows')
                                  end if
                                  if (rank == 0) print *, 'Node-local layout: permutation check passed'
                                  call VecDestroy(locVec, ierr)
                                end subroutine BuildNaturalScatter


                                subroutine FreeNodeLocalLayout()
                                  implicit none
                                  call VecScatterDestroy(natToLocal, ierr)
                                  call VecDestroy(natVec, ierr)
                                  if (allocated(myBoxes)) deallocate(myBoxes)
                                end subroutine FreeNodeLocalLayout


                                Subroutine InitializeBlocks
                                  implicit none

                                  ! =========================
                                  ! LOCAL VARIABLES OF InitializeBlocks
                                  ! =========================
                                  integer ::  blk_i, nBlks

                                  ! =========================
                                  ! LOCAL VARIABLES FOR SHAT RANGES AND SETTING IT AS 0
                                  ! =========================
                                  PetscInt :: rstart, rend, cstart, cend, row, col





                                  ! Calculate number of overlapping points
                                  nOvlp = int(percentageOfOverall * nDFT)
                                  nFreq = nDFT / 2 + 1

                                  ! =========================
                                  ! OUTPUT FOR CONFIRMATION
                                  ! =========================
                                  if (rank == 0) then

                                      print *, "numberOfDataInSnapShot = ", numberOfDataInSnapShot
                                      print *, "numberOfSnapShot      = ", numberOfSnapShot
                                      print *, "nDFT                 = ", nDFT
                                      print *, "nFreq                 = ", nFreq
                                      print *, "percentageOfOverall  = ", percentageOfOverall
                                      print *, "nOvlp                = ", nOvlp
                                      print *, "filename             = ", filenameFormat

                                  end if

                                  ! =========================
                                  ! LOGIC TO DECIDE WHEN THE BLOCKS WILL BE UPDATES
                                  ! BECAUSE SNAPSHOT CAN BE IN OVERLAP
                                  ! =========================
                                  current_Block = 1
                                  nBlks = 2
                                  dn = nDFT - nOvlp
                                  allocate(t_idx(nBlks))
                                  t_idx = 0
                                  do blk_i = 1, nBlks
                                    t_idx(blk_i) = t_idx(blk_i) - (blk_i - 1)*dn + 1
                                  end do

                                  ! =========================
                                  ! Every matrix with one row per grid point uses the node-local row
                                  ! count nLocalRows (was PETSC_DECIDE). Columns are unchanged.
                                  ! =========================

                                  ! =========================
                                  ! FIRST BLOCK
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, X_sum_BlockA, ierr)
                                  call MatSetSizes(X_sum_BlockA, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, nFreq, ierr)
                                  call MatSetType(X_sum_BlockA, MATMPIDENSE , ierr)
                                  call MatSetUp(X_sum_BlockA, ierr)

                                  ! =========================
                                  ! SECOND BLOCK
                                  ! =========================
                                  if (nOvlp > 0) then
                                    call MatCreate(PETSC_COMM_WORLD, X_sum_BlockB, ierr)
                                    call MatSetSizes(X_sum_BlockB, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, nFreq, ierr)
                                    call MatSetType(X_sum_BlockB, MATMPIDENSE , ierr)
                                    call MatSetUp(X_sum_BlockB, ierr)
                                  end if

                                  ! =========================
                                  ! USED TO STORE THE SNAPSHOT
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, x_new, ierr)
                                  call MatSetSizes(x_new, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, 1, ierr)
                                  call MatSetType(x_new, MATMPIDENSE, ierr)
                                  call MatSetFromOptions(x_new, ierr)
                                  call MatSetUp(x_new, ierr)

                                  ! =========================
                                  ! USED TO STORE THE MEAN
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, mu, ierr)
                                  call MatSetSizes(mu, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, 1, ierr)
                                  call MatSetType(mu, MATDENSE, ierr)
                                  call MatSetFromOptions(mu, ierr)
                                  call MatSetUp(mu, ierr)
                                  call MatZeroEntries(mu, ierr)
                                  call MatAssemblyBegin(mu, MAT_FINAL_ASSEMBLY, ierr)
                                  call MatAssemblyEnd(mu, MAT_FINAL_ASSEMBLY, ierr)

                                  ! =========================
                                  ! USED TO STORE THE LEFT SINGULAR VECTOR
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, U_hat, ierr)
                                  call MatSetSizes(U_hat, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, k*nFreq, ierr)
                                  call MatSetType(U_hat, MATDENSE, ierr)
                                  call MatSetUp(U_hat, ierr)

                                  ! =========================
                                  ! USED TO STORE THE EIGEN VALUES
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, S_hat, ierr)
                                  call MatSetSizes(S_hat, PETSC_DECIDE, PETSC_DECIDE, k, nFreq, ierr)
                                  call MatSetType(S_hat, MATAIJ, ierr)
                                  call MatSetUp(S_hat, ierr)
                                  call MatGetOwnershipRange(S_hat, rstart, rend, ierr)
                                  call MatGetOwnershipRangeColumn(S_hat, cstart, cend, ierr)
                                  do row = rstart, rend - 1
                                    do col = 0, nFreq - 1
                                        call MatSetValue(S_hat, row, col, (0.0d0,0.0d0), INSERT_VALUES, ierr)
                                    end do
                                  end do
                                  call MatAssemblyBegin(S_hat, MAT_FINAL_ASSEMBLY, ierr)
                                  call MatAssemblyEnd(S_hat, MAT_FINAL_ASSEMBLY, ierr)



                                End Subroutine InitializeBlocks


                                Subroutine ReadWeightMatrix
                                  implicit none

                                  ! =========================
                                  ! LOCAL VARIABLES OF ReadWeightMatrix
                                  ! =========================
                                  integer :: unitno
                                  real(kind=8), allocatable :: w(:)
                                  PetscScalar, allocatable :: sqrtWArrayF90(:)
                                  integer, allocatable :: indices(:)
                                  integer :: i

                                  call VecCreate(PETSC_COMM_WORLD, sqrtW, ierr)
                                  call VecSetSizes(sqrtW, nLocalRows, numberOfDataInSnapShot, ierr)
                                  call VecSetFromOptions(sqrtW, ierr)

                                  if (trim(weightMatrixFilePath) /= '') then
                                    ! The file is in the NATURAL order, as before. Rank 0 reads it into
                                    ! the natural-order vector exactly as before ...
                                    if (rank == 0) then
                                        allocate(indices(numberOfDataInSnapShot))
                                        allocate(w(numberOfDataInSnapShot))
                                        allocate(sqrtWArrayF90(numberOfDataInSnapShot))

                                        ! Open and read the weight matrix file
                                        unitno = 10
                                        open(unit=unitno, file=trim(weightMatrixFilePath), status='old')
                                        read(unitno, *) w(:)
                                        close(unitno)

                                        ! Compute the square root of the weights
                                        sqrtWArrayF90 = sqrt(w)

                                        ! Fill indices array (Fortran is 1-based, PETSc is 0-based)
                                        do i = 1, numberOfDataInSnapShot
                                            indices(i) = i - 1
                                        end do

                                        ! Set values in the natural-order PETSc vector
                                        call VecSetValues(natVec, numberOfDataInSnapShot, indices, sqrtWArrayF90, &
                                                          INSERT_VALUES, ierr)

                                        ! Deallocate arrays
                                        deallocate(w, sqrtWArrayF90, indices)
                                      end if
                                      call VecAssemblyBegin(natVec, ierr)
                                      call VecAssemblyEnd(natVec, ierr)
                                      ! ... then each weight moves to the rank and row holding that grid point
                                      call VecScatterBegin(natToLocal, natVec, sqrtW, INSERT_VALUES, SCATTER_FORWARD, ierr)
                                      call VecScatterEnd(natToLocal, natVec, sqrtW, INSERT_VALUES, SCATTER_FORWARD, ierr)
                                    else
                                        ! If no file is provided, set all entries to 1.0 (complex)
                                        call VecSet(sqrtW, cmplx(1.0, 0.0, kind=c_double), ierr)
                                    end if

                                  call VecAssemblyBegin(sqrtW, ierr)
                                  call VecAssemblyEnd(sqrtW, ierr)

                                End Subroutine ReadWeightMatrix


                                subroutine WindowingFunction()
                                  implicit none

                                  integer :: i
                                  real :: sum_window

                                  ! =========================
                                  ! LOCAL VARIABLES OF WindowingFunction
                                  ! =========================

                                  allocate (window(nDFT))


                                  if (rank == 0) then
                                    do i = 1, nDFT
                                      if (useHamming) then
                                        window(i) = 0.54 - 0.46 * cos(2.0 * 3.141592653589793 * (i-1) / (nDFT-1))
                                      else
                                        window(i) = 1.0
                                      end if

                                    end do
                                    sum_window = sum(window)

                                    winWeight = 1.0 / (sum_window / real(nDFT))

                                    print *, 'winWeight on RANK 0', winWeight

                                    call MPI_Bcast(window , nDFT * 2, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                            print *, 'Error broadcasting matrix size.'
                                            call PetscFinalize(ierr)
                                            call SlepcFinalize(ierr)
                                            stop
                                        end if

                                    call MPI_Bcast(winWeight, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                            print *, 'Error broadcasting matrix size.'
                                            call PetscFinalize(ierr)
                                            call SlepcFinalize(ierr)
                                            stop
                                        end if

                                  else
                                    call MPI_Bcast(window , nDFT * 2, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
                                    if (ierr /= 0) then
                                        print *, 'Error broadcasting matrix size.'
                                        call PetscFinalize(ierr)
                                        call SlepcFinalize(ierr)
                                        stop
                                    end if

                                    call MPI_Bcast(winWeight, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                    if (ierr /= 0) then
                                        print *, 'Error broadcasting matrix size.'
                                        call PetscFinalize(ierr)
                                        call SlepcFinalize(ierr)
                                        stop
                                    end if

                                  end if
                                end subroutine WindowingFunction


                                subroutine PerformDFT()
                                  use, intrinsic :: iso_c_binding

                                  implicit none
                                  include 'fftw3-mpi.f03'

                                  ! FFTW variables
                                  complex(c_double), allocatable :: data_in(:,:)
                                  complex(c_double), allocatable :: data_out(:,:)
                                  type(C_PTR) :: plan_ptr


                                  integer :: local_start, local_end
                                  integer :: i, j

                                  call fftw_mpi_init()

                                  call MatCreate(PETSC_COMM_WORLD, FourierRange, ierr)
                                  call MatSetSizes(FourierRange, PETSC_DECIDE, PETSC_DECIDE, nDFT, 1, ierr)
                                  call MatSetType(FourierRange, MATAIJ, ierr)
                                  call MatSetFromOptions(FourierRange, ierr)
                                  call MatSetUp(FourierRange, ierr)
                                  call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)


                                  allocate(fourierF90Array(local_end - local_start,nFreq))

                                  fourierF90Array = 0.0

                                  allocate(data_in(nDFT, local_end - local_start))
                                  allocate(data_out(nDFT, local_end - local_start))
                                  data_in = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                  do i = 1, nDFT
                                      do j = local_start + 1, local_end
                                        if (i == j) then
                                            data_in(i, j - local_start) = 1.0_C_DOUBLE  ! Assigning diagonal elements
                                        end if
                                      end do
                                  end do
                                  do j = 1, local_end - local_start
                                    plan_ptr = fftw_plan_dft_1d(nDFT, data_in(:, j), data_out(:, j), FFTW_FORWARD, FFTW_ESTIMATE)
                                    call dfftw_execute(plan_ptr)
                                    call fftw_destroy_plan(plan_ptr)
                                  end do

                                  do i = local_start, local_end - 1
                                    do j = 1, nFreq

                                       if (j == 1) then

                                          fourierF90Array(i - local_start + 1,j) = cmplx(real(data_out(j,i - local_start + 1), kind=kind(0.0)), &
                                          aimag(data_out(j,i - local_start + 1)), kind=kind((0.0, 0.0)))

                                       else

                                          fourierF90Array(i - local_start + 1,j) = 2 * cmplx(real(data_out(j,i - local_start + 1), kind=kind(0.0)), &
                                          aimag(data_out(j,i - local_start + 1)), kind=kind((0.0, 0.0)))

                                       end if
                                    end do
                                  end do
                                  ! if (rank == 0) then

                                  !   print *, fourierF90Array(1,:)
                                  !   stop
                                  ! end if
                                  deallocate(data_out,data_in)

                                end subroutine PerformDFT



                                ! ================================================================
                                ! ADIOS2 SST INPUT
                                ! ================================================================
                                subroutine stream_abort(msg)
                                  implicit none
                                  character(len=*), intent(in) :: msg
                                  integer :: e
                                  print *, 'SST INPUT ERROR (rank ', rank, '): ', trim(msg)
                                  call MPI_Abort(PETSC_COMM_WORLD, 1, e)
                                end subroutine stream_abort


                                ! Open the stream and allocate one receive buffer per (box, variable).
                                subroutine OpenStream()
                                  implicit none
                                  integer :: v, ib, aerr
                                  logical :: haveConfig

                                  inquire(file=trim(adiosConfigFile), exist=haveConfig)
                                  if (haveConfig) then
                                    call adios2_init(adios, trim(adiosConfigFile), PETSC_COMM_WORLD, aerr)
                                  else
                                    call adios2_init(adios, PETSC_COMM_WORLD, aerr)
                                  end if
                                  if (aerr /= 0) call stream_abort('adios2_init failed')

                                  call adios2_declare_io(ioIn, adios, 'spod-input', aerr)
                                  if (aerr /= 0) call stream_abort('adios2_declare_io failed')
                                  if (.not. haveConfig) then
                                    ! no XML: hard-coded SST reader
                                    call adios2_set_engine(ioIn, 'SST', aerr)
                                    if (aerr /= 0) call stream_abort('adios2_set_engine failed')
                                  end if

                                  if (rank == 0) then
                                    print *, 'ADIOS2 input: config = ', trim(adiosConfigFile), ' (found: ', haveConfig, ')'
                                    print *, 'ADIOS2 input: opening stream ', trim(streamName), ' engine ', trim(ioIn%engine_type)
                                  end if
                                  call adios2_open(engineIn, ioIn, trim(streamName), adios2_mode_read, PETSC_COMM_WORLD, aerr)
                                  if (aerr /= 0) call stream_abort('adios2_open failed for stream '//trim(streamName))

                                  allocate(getBuf(max(1, nMyBoxes), nStreamVars))
                                  do v = 1, nStreamVars
                                    do ib = 1, nMyBoxes
                                      allocate(getBuf(ib, v)%v(box_points(myBoxes(ib))))
                                    end do
                                  end do
                                end subroutine OpenStream


                                ! Read one step (= one snapshot) into x_new.
                                ! gotStep = .false. at end of stream.
                                Subroutine Load_Snapshot(gotStep)
                                  implicit none
                                  logical, intent(out) :: gotStep

                                  type(adios2_variable) :: var
                                  integer :: stepStatus, aerr, v, ndims, ib
                                  integer(kind=8), allocatable :: shp(:)
                                  integer(kind=8) :: m, np, row
                                  PetscScalar, pointer :: local_x_new(:,:)
                                  double precision :: t0, t1, t2
                                  character(len=200) :: msg

                                  gotStep = .false.
                                  t0 = MPI_Wtime()
                                  call adios2_begin_step(engineIn, adios2_step_mode_read, beginStepTimeout, stepStatus, aerr)
                                  if (aerr /= 0) call stream_abort('adios2_begin_step failed')
                                  if (stepStatus == adios2_step_status_end_of_stream) return
                                  if (stepStatus /= adios2_step_status_ok) call stream_abort('begin_step: not ready / timeout')
                                  t1 = MPI_Wtime()

                                  do v = 1, nStreamVars
                                    call adios2_inquire_variable(var, ioIn, trim(streamVarNames(v)), aerr)
                                    if (aerr /= 0 .or. .not. var%valid) &
                                      call stream_abort('variable '//trim(streamVarNames(v))//' not in this step')
                                    if (var%type /= adios2_type_dp) &
                                      call stream_abort('variable '//trim(streamVarNames(v))//' is not real(8)')
                                    call adios2_variable_shape(shp, ndims, var, aerr)
                                    if (ndims /= 3) call stream_abort('variable '//trim(streamVarNames(v))//' is not 3-D')
                                    if (shp(1) /= gridNx .or. shp(2) /= gridNy .or. shp(3) /= gridNz) then
                                      write(msg, '(A,A,A,3I8,A,3I8)') 'shape of ', trim(streamVarNames(v)), ' is', shp(1:3), &
                                        ', hard-coded', gridNx, gridNy, gridNz
                                      call stream_abort(msg)
                                    end if
                                    ! One deferred get per box. Each box lies inside the pencil of one
                                    ! Xcompact3d rank on this node, so ADIOS2 fetches it from that rank only.
                                    ! The selection is captured at each get, so the variable can be reused.
                                    do ib = 1, nMyBoxes
                                      call adios2_set_selection(var, 3, myBoxes(ib)%start, myBoxes(ib)%count, aerr)
                                      if (aerr /= 0) call stream_abort('adios2_set_selection failed')
                                      call adios2_get(engineIn, var, getBuf(ib, v)%v, adios2_mode_deferred, aerr)
                                      if (aerr /= 0) call stream_abort('adios2_get failed')
                                    end do
                                  end do
                                  call adios2_end_step(engineIn, aerr)     ! deferred gets are completed here
                                  if (aerr /= 0) call stream_abort('adios2_end_step failed')

                                  ! real -> complex copy into this rank's rows of x_new
                                  ! (variable-major, then box, same order as BuildNaturalScatter)
                                  call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                  row = 0
                                  do v = 1, nStreamVars
                                    do ib = 1, nMyBoxes
                                      np = box_points(myBoxes(ib))
                                      do m = 1, np
                                        local_x_new(row + m, 1) = cmplx(getBuf(ib, v)%v(m), 0.0d0, kind=8)
                                      end do
                                      row = row + np
                                    end do
                                  end do
                                  call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                                  call MatAssemblyBegin(x_new, MAT_FINAL_ASSEMBLY, ierr)
                                  call MatAssemblyEnd(x_new, MAT_FINAL_ASSEMBLY, ierr)
                                  t2 = MPI_Wtime()

                                  tWaitStep = tWaitStep + (t1 - t0)
                                  tRead     = tRead     + (t2 - t1)
                                  gotStep = .true.
                                End Subroutine Load_Snapshot


                                subroutine CloseStream()
                                  implicit none
                                  integer :: aerr, v, ib
                                  call adios2_close(engineIn, aerr)
                                  call adios2_finalize(adios, aerr)
                                  if (allocated(getBuf)) then
                                    do v = 1, nStreamVars
                                      do ib = 1, nMyBoxes
                                        if (allocated(getBuf(ib, v)%v)) deallocate(getBuf(ib, v)%v)
                                      end do
                                    end do
                                    deallocate(getBuf)
                                  end if
                                end subroutine CloseStream


                                subroutine handle_Latest_Snapshot()

                                  implicit none

                                  PetscScalar, pointer :: local_x_new(:,:)
                                  PetscScalar, pointer :: local_x_sum(:,:)
                                  PetscScalar, pointer :: local_mu(:,:)

                                  PetscInt :: localRows,localCols
                                  PetscInt :: i

                                  ! To decide which process has current DFT row
                                  logical :: i_own_row
                                  integer :: row_owner, global_owner
                                  integer :: global_row_idx, local_row_idx
                                  PetscInt :: local_start, local_end
                                  integer ::  j, ld
                                  PetscInt :: rstart, rend
                                  PetscScalar, allocatable ::  global_row_data(:)

                                  PetscScalar :: alpha, beta

                                  integer ::  min_t_idx


                                  if (rank == 0) then
                                    t_i = t_i + 1

                                    call MPI_Bcast(t_i, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                           print *, 'Error broadcasting matrix size.'
                                           call PetscFinalize(ierr)
                                           call SlepcFinalize(ierr)
                                           stop
                                        end if
                                  else
                                    call MPI_Bcast(t_i, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                          print *, 'Error broadcasting matrix size.'
                                          call PetscFinalize(ierr)
                                          call SlepcFinalize(ierr)
                                          stop
                                        end if
                                  end if

                                  call MatGetLocalSize(mu, localRows, localCols, ierr)

                                  call MatDenseGetArrayF90(mu, local_mu, ierr)
                                  call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                  alpha = t_i - 1.0d0
                                  beta = 1.0d0 / t_i

                                  do i = 1, localRows
                                    local_mu(i, 1) = beta * (alpha * local_mu(i, 1) + local_x_new(i, 1))
                                  end do

                                  call MatDenseRestoreArrayF90(mu, local_mu, ierr)
                                  call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)




                                  if (t_idx(1) > 0) then

                                    call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)

                                    global_row_idx = t_idx(1) - 1

                                    i_own_row = (global_row_idx >= local_start) .and. (global_row_idx < local_end)

                                    allocate(global_row_data(nFreq))
                                    if (i_own_row) then
                                        row_owner = rank
                                    else
                                        row_owner = -1
                                    end if
                                    call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)

                                    if (i_own_row) then
                                      local_row_idx = global_row_idx - local_start + 1
                                      global_row_data = fourierF90Array(local_row_idx, :)
                                    else
                                        global_row_data = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                    end if

                                    call MPI_Bcast(global_row_data, 2*nFreq, MPI_DOUBLE_PRECISION, global_owner, PETSC_COMM_WORLD, ierr)
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)

                                    call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)



                                    call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                    call MatDenseGetArrayF90(X_sum_BlockA, local_x_sum, ierr)

                                    call MatGetOwnershipRange(X_sum_BlockA, rstart, rend, ierr)
                                    ld = rend - rstart

                                    do j = 1, nFreq
                                      do i = 1, ld
                                        local_x_sum(i, j) = local_x_sum(i, j) + (local_x_new(i, 1) * window(t_idx(1))) * global_row_data(j)
                                      end do
                                    end do

                                    ! if (rank == 0) then

                                    !   print *, 'local_x_sum:', local_x_sum(1, 1), local_x_new(1, 1) , window(t_idx(1)), global_row_data(1)

                                    ! end if

                                    call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                                    call MatDenseRestoreArrayF90(X_sum_BlockA, local_x_sum, ierr)
                                    deallocate(global_row_data)


                                  end if


                                  if (nOvlp > 0) then
                                    if (t_idx(2) > 0) then

                                      call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)

                                      global_row_idx = t_idx(2) - 1
                                      i_own_row = (global_row_idx >= local_start) .and. (global_row_idx < local_end)

                                      allocate(global_row_data(nFreq))
                                      if (i_own_row) then
                                          row_owner = rank
                                      else
                                          row_owner = -1
                                      end if

                                      call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)

                                      if (i_own_row) then
                                          ! Calculate local index of the global row
                                          local_row_idx = global_row_idx - local_start + 1
                                          global_row_data = fourierF90Array(local_row_idx, :)

                                      else
                                          global_row_data = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                      end if

                                      call MPI_Bcast(global_row_data, 2*nFreq, MPI_DOUBLE_PRECISION, global_owner, PETSC_COMM_WORLD, ierr)
                                      call MPI_Barrier(PETSC_COMM_WORLD, ierr)


                                      call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                      call MatDenseGetArrayF90(X_sum_BlockB, local_x_sum, ierr)

                                      call MatGetOwnershipRange(X_sum_BlockB, rstart, rend, ierr)
                                      ld = rend - rstart  ! Number of locally owned rows

                                      do j = 1, nFreq
                                        do i = 1, ld
                                          local_x_sum(i, j) = local_x_sum(i, j) + (local_x_new(i, 1) * window(t_idx(2))) * global_row_data(j)
                                        end do
                                      end do


                                      call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                                      call MatDenseRestoreArrayF90(X_sum_BlockB, local_x_sum, ierr)
                                      deallocate(global_row_data)
                                    end if

                                  end if


                                  if (t_idx(1) == nDFT) then
                                    min_t_idx = minval(t_idx)
                                    t_idx(1) = min_t_idx - dn
                                    call handle_Latest_Block(1,X_sum_BlockA)
                                    call MatZeroEntries(X_sum_BlockA, ierr)
                                    call MatAssemblyBegin(X_sum_BlockA, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(X_sum_BlockA, MAT_FINAL_ASSEMBLY, ierr)
                                  else
                                      t_idx(1) = t_idx(1) + 1
                                  end if

                                  if (nOvlp > 0) then
                                    if (t_idx(2) == nDFT) then
                                      min_t_idx = minval(t_idx)
                                      t_idx(2) = min_t_idx - dn
                                      call handle_Latest_Block(2,X_sum_BlockB)
                                      call MatZeroEntries(X_sum_BlockB, ierr)
                                      call MatAssemblyBegin(X_sum_BlockB, MAT_FINAL_ASSEMBLY, ierr)
                                      call MatAssemblyEnd(X_sum_BlockB, MAT_FINAL_ASSEMBLY, ierr)
                                    else
                                        t_idx(2) = t_idx(2) + 1
                                    end if
                                  end if


                                end subroutine handle_Latest_Snapshot


                                subroutine handle_Latest_Block(blk_j,XHat)

                                  implicit none

                                  ! =========================
                                  ! blk_j and XHat WILL BE PASSED BY handle_Latest_Snapshot
                                  ! =========================
                                  integer :: blk_j
                                  Mat, intent(inout) :: XHat

                                  ! =========================
                                  ! LOCAL VARIABLES FOR MEAN REMOVAL
                                  ! =========================
                                  integer ::   ld
                                  PetscScalar, pointer :: local_b(:,:)
                                  PetscScalar, pointer :: local_xhat(:,:)
                                  PetscScalar, allocatable :: global_row_data(:)
                                  integer :: row_idx, i,j , rstart, rend, local_start, local_end
                                  PetscScalar :: scaling_factor
                                  logical :: i_own_row
                                  integer :: row_owner, global_owner
                                  integer :: global_row_idx, local_row_idx

                                  ! =========================
                                  ! LOCAL VARIABLES FOR FREQUENCY
                                  ! =========================
                                  integer :: fi
                                  PetscInt :: col_start_fi, rows_local_fi
                                  PetscInt ::  local_size
                                  PetscInt :: row, col
                                  PetscScalar :: valueOfCmatrixArray(1)
                                  PetscScalar, pointer :: u_hat_array(:,:)
                                  PetscScalar, pointer :: u_array_copy(:,:)
                                  PetscScalar, pointer :: x_hat_col_array(:)
                                  PetscScalar, pointer :: u_p_array(:)
                                  PetscScalar, pointer :: u_new_array(:)

                                  PetscReal abs_up
                                  PetscReal :: factor
                                  PetscScalar ::  scale_factor

                                  ! =========================
                                  ! MATRIX VARIABLES FOR FREQUENCY
                                  ! =========================
                                  Mat :: U
                                  Mat :: Kmode

                                  ! =========================
                                  ! VECTORS VARIABLES FOR FREQUENCY
                                  ! =========================
                                  Vec :: X_hat_column
                                  Vec :: x
                                  Vec :: u_p
                                  Vec :: u_new
                                  Vec :: Ux
                                  Vec :: temp
                                  Vec :: U_col
                                  Vec :: V_col

                                  ! =========================
                                  ! LOCAL VARIABLES FOR SVD
                                  ! =========================
                                  SVD :: svd
                                  PetscReal :: tolerance
                                  PetscInt :: maxIterations
                                  PetscInt :: nconverged, local_array_i
                                  PetscScalar, pointer :: array(:)
                                  PetscReal :: sigma
                                  Mat :: Umat


                                  ! =========================
                                  ! LOCAL VARIABLES FOR RECONSTRUCTION OF SPOD MODE
                                  ! =========================
                                  PetscInt ::  n, m_local, n_local, vec_local_size
                                  Mat :: Utemp
                                  Mat :: UtempMult
                                  PetscInt :: cstart
                                  PetscScalar, pointer :: u_array(:,:), utemp_array(:,:)
                                  PetscScalar, pointer :: vec_array(:)
                                  PetscScalar, pointer ::  uhat_array(:,:)
                                  integer :: nrowsForUtemp


                                  ! =========================
                                  ! VIEWER
                                  ! =========================
                                  PetscViewer :: viewer


                                  current_Block = current_Block + 1
                                  tolerance = 1.0d-30
                                  maxIterations = 10000

                                  if (rank == 0) then
                                    print *, '-------------------------------------------'
                                    print *, 'Processing Block'
                                    print *, '-------------------------------------------'
                                  end if

                                  call MatDenseGetArrayF90(mu, local_b, ierr)

                                  do row_idx = 1, nDFT
                                    if (rank == 0) then
                                      print *, 'Processing Block Realization', row_idx
                                    end if
                                    call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)

                                    global_row_idx = row_idx - 1

                                    i_own_row = (global_row_idx >= local_start) .and. (global_row_idx < local_end)

                                    allocate(global_row_data(nFreq))
                                    if (i_own_row) then
                                        row_owner = rank
                                    else
                                        row_owner = -1
                                    end if

                                    call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)

                                    if (i_own_row) then
                                        ! Calculate local index of the global row
                                        local_row_idx = global_row_idx - local_start + 1
                                        global_row_data = fourierF90Array(local_row_idx, :)

                                    else
                                        global_row_data = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                    end if

                                    call MPI_Bcast(global_row_data, 2*nFreq, MPI_DOUBLE_PRECISION, global_owner, PETSC_COMM_WORLD, ierr)
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)

                                    call MatDenseGetArrayF90(XHat, local_xhat, ierr)
                                    scaling_factor = -window(row_idx)
                                    call MatGetOwnershipRange(XHat, rstart, rend, ierr)
                                    ld = rend - rstart  ! Number of locally owned rows

                                    do j = 1, nFreq
                                      do i = 1, ld
                                        local_xhat(i, j) = local_xhat(i, j) + ((scaling_factor * local_b(i, 1)) * global_row_data(j))
                                      end do
                                    end do

                                    call MatDenseRestoreArrayF90(XHat, local_xhat, ierr)

                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)

                                    deallocate(global_row_data)

                                  end do

                                  call MatDenseRestoreArrayF90(mu, local_b, ierr)
                                  call MatDenseGetArrayF90(XHat, local_xhat, ierr)
                                  scaling_factor = cmplx(winWeight / real(nDFT), 0.0d0, kind=8)

                                  call MatGetOwnershipRange(XHat, rstart, rend, ierr)
                                  ld = rend - rstart  ! Number of locally owned rows

                                  do j = 1, nFreq
                                    do i = 1, ld
                                      local_xhat(i, j) = local_xhat(i, j) * scaling_factor
                                    end do
                                  end do

                                  call MatDenseRestoreArrayF90(XHat, local_xhat, ierr)
                                  call MatDenseGetArrayF90(XHat, local_xhat, ierr)
                                  if (rank == 0) then
                                    print *, '-------------------------------------------'
                                    print *, 'DONE : Processing Block'
                                    print *, '-------------------------------------------'
                                  end if


                                  ! =========================
                                  ! FREQUENCY VARIABLE INITIALISE
                                  ! (objects with one row per grid point use nLocalRows)
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, U, ierr)
                                  call MatSetSizes(U, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, k, ierr)
                                  call MatSetType(U, MATDENSE, ierr)
                                  call MatSetUp(U, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, X_hat_column, ierr)
                                  call VecSetSizes(X_hat_column, nLocalRows, numberOfDataInSnapShot, ierr)
                                  call VecSetType(X_hat_column, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, x, ierr)
                                  call VecSetSizes(x, nLocalRows, numberOfDataInSnapShot, ierr)
                                  call VecSetType(x, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, u_new, ierr)
                                  call VecSetSizes(u_new, nLocalRows, numberOfDataInSnapShot, ierr)
                                  call VecSetType(u_new, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, u_p, ierr)
                                  call VecSetSizes(u_p, nLocalRows, numberOfDataInSnapShot, ierr)
                                  call VecSetType(u_p, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, Ux, ierr)
                                  call VecSetSizes(Ux, PETSC_DECIDE, k, ierr)
                                  call VecSetType(Ux, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, temp, ierr)
                                  call VecSetSizes(temp, nLocalRows, numberOfDataInSnapShot, ierr)
                                  call VecSetType(temp, VECSTANDARD, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, Kmode, ierr)
                                  call MatSetSizes(Kmode, PETSC_DECIDE, PETSC_DECIDE, k + 1, k + 1, ierr)
                                  call MatSetType(Kmode, MATAIJ, ierr)
                                  call MatSetUp(Kmode, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, Umat, ierr)
                                  call MatSetSizes(Umat, PETSC_DECIDE, PETSC_DECIDE, k + 1, k + 1, ierr)
                                  call MatSetType(Umat, MATDENSE, ierr)
                                  call MatSetUp(Umat, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, Utemp, ierr)
                                  call MatSetSizes(Utemp, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, k+1, ierr)
                                  call MatSetType(Utemp, MATDENSE, ierr)
                                  call MatSetUp(Utemp, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, UtempMult, ierr)
                                  call MatSetSizes(UtempMult, nLocalRows, PETSC_DECIDE, numberOfDataInSnapShot, k+1, ierr)
                                  call MatSetType(UtempMult, MATDENSE, ierr)
                                  call MatSetUp(UtempMult, ierr)

                                  ! =========================
                                  ! FREQUENCY CALCULATIONS
                                  ! =========================

                                  if (rank == 0) then
                                    print *, '-------------------------------------------'
                                    print *, 'Processing Each Frequency'
                                    print *, '-------------------------------------------'
                                  end if
                                  call MatDenseGetArrayF90(XHat, local_xhat, ierr)


                                  do fi = 1, nFreq
                                    if (rank == 0) then
                                      print *, 'Processing Frequency', fi
                                    end if
                                  ! =========================
                                  ! Get U from U_hat and X column from Xhat realization
                                  ! =========================
                                    call MatGetOwnershipRange(U_hat, rstart, rend, ierr)
                                    rows_local_fi = rend - rstart
                                    col_start_fi = (fi-1) * k

                                    call MatDenseGetArrayF90(U_hat, u_hat_array, ierr)
                                    call MatDenseGetArrayF90(U, u_array_copy, ierr)

                                    do j = 1, k
                                      do i = 1, rows_local_fi
                                        u_array_copy(i, j) = u_hat_array(i, col_start_fi + j)
                                      end do
                                    end do

                                    call MatDenseRestoreArrayF90(U_hat, u_hat_array, ierr)
                                    call MatDenseRestoreArrayF90(U, u_array_copy, ierr)
                                    call MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY, ierr)

                                    call VecGetArrayF90(X_hat_column, x_hat_col_array, ierr)
                                    call MatGetOwnershipRange(XHat, rstart, rend, ierr)
                                    ld = rend - rstart

                                    do i = 1, ld
                                      x_hat_col_array(i) = local_xhat(i, fi)
                                    end do

                                    call VecRestoreArrayF90(X_hat_column, x_hat_col_array, ierr)


                                    call VecPointwiseMult(x, X_hat_column, sqrtW, ierr)
                                    call MatMultHermitianTranspose(U, x, Ux, ierr)
                                    call MatMult(U, Ux, temp, ierr)
                                    call VecWAXPY(u_p, (-1.0d0, 0.0d0), temp, x, ierr)
                                    call VecNorm(u_p, NORM_2, abs_up, ierr)

                                    call VecGetArrayF90(u_p, u_p_array, ierr)
                                    call VecGetArrayF90(u_new, u_new_array, ierr)
                                    call VecGetLocalSize(u_new, local_size, ierr)

                                    do i = 1, local_size
                                        u_new_array(i) = u_p_array(i)
                                    end do

                                    call VecRestoreArrayF90(u_p, u_p_array, ierr)
                                    call VecRestoreArrayF90(u_new, u_new_array, ierr)

                                    scale_factor = cmplx(1.0d0 / abs_up, 0.0d0, kind=8)
                                    call VecScale(u_new,scale_factor , ierr)


                                  ! ===================================================
                                  ! CONSTRUCT K MATRIX
                                  ! ===================================================
                                    factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * sqrt(real(current_Block + 1))

                                    call MatGetOwnershipRange(S_hat, rstart, rend, ierr)

                                    do row = rstart, rend - 1
                                        if (row < k) then
                                          do col = 0, k - 1
                                            if (row == col) then
                                            call MatGetValues(S_hat, 1, [row], 1, [fi - 1], valueOfCmatrixArray, ierr)

                                              call MatSetValue(Kmode, row, row, factor * valueOfCmatrixArray(1), INSERT_VALUES, ierr)
                                            else
                                              call MatSetValue(Kmode, row, col, cmplx(0.0d0, 0.0d0, kind=8), INSERT_VALUES, ierr)

                                            end if
                                          end do
                                          call VecGetValues(Ux, 1, [row], valueOfCmatrixArray, ierr)

                                          scale_factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * valueOfCmatrixArray(1)
                                          call MatSetValue(Kmode, row, k,scale_factor , INSERT_VALUES, ierr)
                                        else if (row == k) then
                                          do col = 0, k - 1
                                              call MatSetValue(Kmode, row, col,cmplx(0.0d0, 0.0d0, kind=8), INSERT_VALUES, ierr)
                                          end do
                                          scale_factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * abs_up
                                          call MatSetValue(Kmode, row, k,scale_factor , INSERT_VALUES, ierr)
                                        end if
                                    end do

                                    if ((rend) == k) then
                                      do col = 0, k - 1
                                        call MatSetValue(Kmode, rend, col,cmplx(0.0d0, 0.0d0, kind=8), INSERT_VALUES, ierr)
                                      end do
                                      scale_factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * abs_up
                                      call MatSetValue(Kmode, rend, k,scale_factor , INSERT_VALUES, ierr)

                                    end if
                                    call MatAssemblyBegin(Kmode, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(Kmode, MAT_FINAL_ASSEMBLY, ierr)



                                    ! call PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, viewer, ierr)
                                    ! call MatView(Kmode, viewer, ierr)

                                  ! ===================================================
                                  ! PERFORM SVD USING SLEPC
                                  ! ===================================================
                                    call MatCreateVecs(Kmode, U_col, V_col, ierr)
                                    call VecAssemblyBegin(U_col, ierr)
                                    call VecAssemblyEnd(U_col, ierr)

                                    call SVDCreate(PETSC_COMM_WORLD, svd, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in SVDCreate: ", ierr
                                      stop
                                    end if
                                    call SVDSetTolerances(svd,tolerance,maxIterations,ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in SVDSetTolerances: ", ierr
                                      stop
                                    end if
                                    PetscCallA(SVDSetOperators(svd,Kmode,PETSC_NULL_MAT,ierr))
                                    PetscCallA(SVDSetProblemType(svd,SVD_STANDARD,ierr))
                                    PetscCallA(SVDSetType(svd,SVDLAPACK ,ierr))
                                    call SVDSetFromOptions(svd, ierr)
                                    call SVDSolve(svd, ierr)
                                    call SVDGetConverged(svd, nconverged, ierr)

                                  ! ===================================================
                                  ! EXTRACT EIGEN VALUES AND VECTORS
                                  ! ===================================================

                                    do i = 0, nconverged - 1
                                      call PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, viewer, ierr)
                                      call SVDGetSingularTriplet(svd, i, sigma, U_col, V_col, ierr)


                                  ! ===================================================
                                  ! UPDATE SHAT WITH LATEST EIGEN VALUES
                                  ! ===================================================
                                      call MatGetOwnershipRange(S_hat, rstart, rend, ierr)
                                      if (rank == 0) then
                                        print*, sigma
                                      end if
                                      do row = rstart, rend - 1
                                        if (row == i) then
                                            call MatSetValue(S_hat, row, fi - 1, cmplx(sigma, 0.0d0, kind=8), INSERT_VALUES, ierr)
                                        end if
                                      end do
                                      call MatAssemblyBegin(S_hat, MAT_FINAL_ASSEMBLY, ierr)
                                      call MatAssemblyEnd(S_hat, MAT_FINAL_ASSEMBLY, ierr)

                                  ! ===================================================
                                  ! STORE EIGEN VECTORS IN UMAT
                                  ! ===================================================
                                      call VecGetOwnershipRange(U_col, rstart, rend, ierr)

                                      allocate(array(rend - rstart))
                                      call VecGetArrayF90(U_col, array, ierr)

                                      local_array_i = 1
                                      do row = rstart, rend - 1

                                      call MatSetValue(Umat, row, i, array(local_array_i), INSERT_VALUES, ierr)
                                      local_array_i = local_array_i + 1

                                      end do
                                      call VecRestoreArrayF90(U_col, array, ierr)

                                    end do

                                    call MatAssemblyBegin(Umat, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(Umat, MAT_FINAL_ASSEMBLY, ierr)
                                    call SVDDestroy(svd, ierr)





                                  ! ===================================================
                                  ! RECONSTRUCT SPOD MODES AND ASSIGN TO GLOBAL MATRIX
                                  ! ===================================================
                                    call MatGetOwnershipRange(Utemp, rstart, rend, ierr)
                                    call MatGetLocalSize(U, m_local, n_local, ierr)
                                    nrowsForUtemp = rend - rstart
                                    n = k

                                    call MatDenseGetArrayF90(U, u_array, ierr)
                                    call MatDenseGetArrayF90(Utemp, utemp_array, ierr)

                                    utemp_array(1:m_local, 1:n) = u_array(1:m_local, 1:n)

                                    call VecGetArrayReadF90(u_new, vec_array, ierr)
                                    call VecGetLocalSize(u_new, vec_local_size, ierr)

                                    utemp_array(1:vec_local_size, n+1) = vec_array(1:vec_local_size)

                                    call VecRestoreArrayReadF90(u_new, vec_array, ierr)
                                    call MatDenseRestoreArrayF90(U, u_array, ierr)
                                    call MatDenseRestoreArrayF90(Utemp, utemp_array, ierr)

                                    call MatAssemblyBegin(Utemp, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(Utemp, MAT_FINAL_ASSEMBLY, ierr)

                                    call MatMatMult(Utemp, Umat, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, UtempMult, ierr)

                                    call MatGetOwnershipRange(UtempMult, rstart, rend, ierr)
                                    m_local = rend - rstart
                                    cstart = (fi - 1) * k

                                    call MatDenseGetArrayF90(UtempMult, utemp_array, ierr)
                                    call MatDenseGetArrayF90(U_hat, uhat_array, ierr)

                                    do j = 1, k
                                      uhat_array(1:m_local, cstart+j) = utemp_array(1:m_local, j)
                                    end do

                                    call MatDenseRestoreArrayF90(UtempMult, utemp_array, ierr)
                                    call MatDenseRestoreArrayF90(U_hat, uhat_array, ierr)

                                    call MatAssemblyBegin(U_hat, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(U_hat, MAT_FINAL_ASSEMBLY, ierr)


                                    call VecDestroy(U_col, ierr)
                                    call VecDestroy(V_col, ierr)


                                  end do


                                  call MatDenseRestoreArrayF90(XHat, local_xhat, ierr)
                                  call MPI_Barrier(MPI_COMM_WORLD, ierr)

                                  call MatDestroy(U, ierr)
                                  call MatDestroy(Kmode, ierr)
                                  call MatDestroy(Umat, ierr)
                                  call MatDestroy(Utemp, ierr)
                                  call MatDestroy(UtempMult, ierr)

                                  call VecDestroy(X_hat_column, ierr)
                                  call VecDestroy(x, ierr)
                                  call VecDestroy(u_p, ierr)
                                  call VecDestroy(u_new, ierr)
                                  call VecDestroy(Ux, ierr)
                                  call VecDestroy(temp, ierr)



                                end subroutine handle_Latest_Block

                                subroutine SPODMode()
                                  implicit none

                                  Vec :: invSqrtW

                                  call VecDuplicate(sqrtW, invSqrtW, ierr)  ! Create a new PETSc vector
                                  call VecCopy(sqrtW, invSqrtW, ierr)       ! Copy sqrtW to invSqrtW
                                  call VecReciprocal(invSqrtW, ierr)        ! Compute 1/sqrt(w)

                                  call MatDiagonalScale(U_hat, invSqrtW, PETSC_NULL_VEC, ierr)  ! Element-wise row division
                                  call MatAssemblyBegin(U_hat, MAT_FINAL_ASSEMBLY, ierr)
                                  call MatAssemblyEnd(U_hat, MAT_FINAL_ASSEMBLY, ierr)

                                end subroutine SPODMode

                                subroutine SaveEnergySpectraAndSPODMode()
                                  implicit none

                                end subroutine SaveEnergySpectraAndSPODMode

                                ! Subroutine to find and save the z highest SPOD modes
                                ! Subroutine to find and save the z highest SPOD modes
                                subroutine SaveHighestSPODModes()
                                  implicit none

                                  integer :: z

                                  PetscInt :: i, j, nrows, rank, size
                                  PetscInt :: cstart, cend
                                  PetscInt, allocatable :: indices(:)
                                  PetscScalar, allocatable :: values(:), local_value
                                  Vec :: XSPOD_VEC
                                  PetscViewer :: viewer
                                  character(len=300) :: filename
                                  PetscInt :: temp_idx
                                  PetscInt :: actual_col

                                  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
                                  call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

                                  z = nFreqToSave
                                  ! Get MPI rank and size
                                  call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'EnegySpectra.txt', viewer, ierr)
                                  call MatView(S_hat, viewer, ierr)
                                  call PetscViewerDestroy(viewer, ierr)
                                  if (rank == 0) then
                                    print *, 'Saved Enegy spectra to File EnegySpectra.txt '
                                  end if



                                  call MatGetSize(S_hat, nrows, nFreq, ierr)

                                  allocate(values(nFreq))
                                  allocate(indices(nFreq))
                                  values = 0.0


                                  if (rank == 0) then
                                    call MatGetOwnershipRangeColumn(S_hat, cstart, cend, ierr)
                                    do i = 1, nFreq
                                      indices(i) = i - 1  ! 0-based indices for PETSc
                                    end do

                                    do j = 1, nFreq
                                        call MatGetValues(S_hat, 1, [0], 1, [j-1], values(j), ierr)

                                    end do

                                    do i = 1, nFreq
                                      do j = 1, nFreq - i
                                        if (abs(values(j)) < abs(values(j+1))) then
                                          ! Swap values
                                          local_value = values(j)
                                          values(j) = values(j+1)
                                          values(j+1) = local_value

                                          ! Swap indices
                                          temp_idx = indices(j)
                                          indices(j) = indices(j+1)
                                          indices(j+1) = temp_idx
                                        end if
                                      end do
                                    end do



                                  end if


                                  call MPI_Barrier(PETSC_COMM_WORLD, ierr)


                                  ! Broadcast the sorted indices and values from rank 0 to all processors
                                  call MPI_Bcast(values, nFreq, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
                                  call MPI_Bcast(indices, nFreq, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

                                  ! Now all processors have the same sorted list

                                  ! Create a vector to hold the mode (node-local layout, like U_hat rows)
                                  call MatCreateVecs(U_hat, PETSC_NULL_VEC, XSPOD_VEC, ierr)
                                  if (ierr /= 0) then
                                    print *, "Error in XSPOD_VEC MatCreateVecs : ", ierr
                                    return
                                  end if

                                  ! Save the top z modes
                                  do i = 1, z
                                    ! Get the column vector from U_hat corresponding to the ith highest mode
                                    actual_col = indices(i) * k

                                    call MatGetColumnVector(U_hat, XSPOD_VEC, actual_col, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in MatGetColumnVector: ", ierr
                                      return
                                    end if

                                    ! Put the mode back in the NATURAL order (same file layout as before:
                                    ! [ux;uy;uz;pp], i fastest). One-off cost, only at the end.
                                    call VecScatterBegin(natToLocal, XSPOD_VEC, natVec, INSERT_VALUES, SCATTER_REVERSE, ierr)
                                    call VecScatterEnd(natToLocal, XSPOD_VEC, natVec, INSERT_VALUES, SCATTER_REVERSE, ierr)

                                    ! Create filename with mode index and corresponding sigma value
                                    write(filename, '(A,A,I0,A,E12.5,A)') trim('SPODMODE'), "_freq_", i, "_sigma_", abs(values(i)), ".dat"

                                    ! Save the vector to a binary file
                                    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in PetscViewerBinaryOpen: ", ierr
                                      return
                                    end if

                                    call VecView(natVec, viewer, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in VecView: ", ierr
                                      return
                                    end if

                                    call PetscViewerDestroy(viewer, ierr)

                                    ! Print information about the saved mode (only from rank 0 to avoid duplicate messages)
                                    if (rank == 0) then
                                      print *, 'Eigen values at index ',indices(i) , ' Globally in Uhat at ', actual_col
                                      print '(A,I0,A,E12.5,A,A)', "Saved SPOD mode ", i, " with sigma value ", abs(values(i)), " to file: ", trim(filename)
                                    end if
                                  end do

                                  call VecDestroy(XSPOD_VEC, ierr)
                                  deallocate(values)
                                  deallocate(indices)

                                end subroutine SaveHighestSPODModes

                                subroutine performStreamingSPOD()
                                  implicit none

                                  integer :: i, nDone
                                  logical :: gotStep
                                  double precision :: t0, tLoc(3), tMax(3)

                                  call OpenStream()

                                  nDone = 0
                                  do i = 1, numberOfSnapShot
                                    call Load_Snapshot(gotStep)
                                    if (.not. gotStep) then
                                      if (rank == 0) print *, 'End of stream after ', nDone, ' snapshots'
                                      exit
                                    end if
                                    if (rank == 0) then
                                      print *, '-------------------------------------------'
                                      print *, 'Processing Snapshot', i
                                      print *, '-------------------------------------------'
                                    end if
                                    currentSnapshotNumber = i
                                    t0 = MPI_Wtime()
                                    call handle_Latest_Snapshot
                                    tCompute = tCompute + (MPI_Wtime() - t0)
                                    nDone = i
                                  end do

                                  call CloseStream()

                                  tLoc = [tWaitStep, tRead, tCompute]
                                  call MPI_Reduce(tLoc, tMax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, 0, PETSC_COMM_WORLD, ierr)
                                  if (rank == 0) then
                                    print *, 'Snapshots consumed           : ', nDone
                                    print *, 'Max over ranks, wait in begin_step [s]: ', tMax(1)
                                    print *, 'Max over ranks, get + copy      [s]: ', tMax(2)
                                    print *, 'Max over ranks, SPOD compute    [s]: ', tMax(3)
                                  end if

                                end subroutine performStreamingSPOD

                              End Module streamingSPODModule

                              Program streamingSPOD
#include <slepc/finclude/slepcsvd.h>

                                use slepc
                                use slepcsvd
                                use streamingSPODModule

                                implicit none


                                call SlepcInitialize(ierr)
                                if (ierr /= 0) then
                                  print *, 'Error initializing SLEPc.'
                                  call SlepcFinalize(ierr)
                                  stop
                                end if
                                ! Get the rank of the calling process in PETSC_COMM_WORLD
                                call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

                                ! Get the total number of processes in PETSC_COMM_WORLD
                                call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

                                ! Print the number of processes from rank 0
                                if (rank == 0) then
                                  print *, "Number of MPI processes (size) =", size
                                end if

                                call UserInput
                                call SetupNodeLocalLayout      ! which boxes each rank reads -> nLocalRows
                                call InitializeBlocks          ! PETSc objects with nLocalRows rows per rank
                                call BuildNaturalScatter       ! natural <-> node-local (weights, modes)
                                call ReadWeightMatrix
                                call WindowingFunction
                                call PerformDFT
                                call performStreamingSPOD
                                call MPI_Barrier(PETSC_COMM_WORLD, ierr)
                                call SPODMode
                                call SaveHighestSPODModes
                                call FreeNodeLocalLayout


                                call SlepcFinalize(ierr)

                              End Program streamingSPOD
