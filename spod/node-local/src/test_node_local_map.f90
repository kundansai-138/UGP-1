! Standalone test of node_local_map (no MPI / PETSc / ADIOS2 needed).
!
!   gfortran -O2 -o test_node_local_map node_local_map.F90 test_node_local_map.f90
!   ./test_node_local_map
!
! For many grids, X3D rank grids, node counts and SPOD rank counts it checks:
!   1. every grid point is read by exactly one SPOD rank,
!   2. every box lies inside the pencil of the X3D rank it names,
!   3. that X3D rank is on the same node as the SPOD rank reading it.
! SPOD ranks are placed both sequentially (block) and round-robin (cyclic);
! the node index is computed like MPI_Comm_split_type + leaders comm do it.
program test_node_local_map
  use node_local_map
  implicit none

  integer, parameter :: grids(3, 3) = reshape([10, 9, 7,  17, 13, 11,  8, 8, 8], [3, 3])
  integer, parameter :: rowcol(2, 4) = reshape([1, 1,  2, 3,  3, 2,  4, 4], [2, 4])
  integer, parameter :: nodesRanks(2, 5) = reshape([1, 1,  2, 4,  2, 5,  3, 7,  4, 32], [2, 5])
  integer :: ig, irc, inr, pencil, ixf, iplc, nPass, nFail
  logical :: ok

  nPass = 0
  nFail = 0
  do ig = 1, 3
    do irc = 1, 4
      if (rowcol(1, irc) > min(grids(1, ig), grids(2, ig)) .or. &
          rowcol(2, irc) > min(grids(2, ig), grids(3, ig))) cycle
      do inr = 1, 5
        do pencil = 1, 3
          do ixf = 0, 1
            do iplc = 0, 1
              ok = check(int(grids(1, ig), i8), int(grids(2, ig), i8), int(grids(3, ig), i8), &
                         rowcol(1, irc), rowcol(2, irc), nodesRanks(1, inr), nodesRanks(2, inr), &
                         pencil, ixf == 1, iplc == 1)
              if (ok) then
                nPass = nPass + 1
              else
                nFail = nFail + 1
                print '(A,3I4,A,2I3,A,2I4,A,I2,A,L2,A,I2)', ' FAIL grid', grids(:, ig), '  prow/pcol', &
                  rowcol(:, irc), '  nodes/spod', nodesRanks(:, inr), '  pencil', pencil, &
                  '  extraFirst', ixf == 1, '  cyclic', iplc
              end if
            end do
          end do
        end do
      end do
    end do
  end do
  print '(A,I0,A,I0,A)', ' node_local_map: ', nPass, ' configurations PASS, ', nFail, ' FAIL'
  if (nFail > 0) error stop 1

contains

  logical function check(nx, ny, nz, prow, pcol, nNodes, nSpod, pencil, extraFirst, cyclic)
    integer(i8), intent(in) :: nx, ny, nz
    integer, intent(in) :: prow, pcol, nNodes, nSpod, pencil
    logical, intent(in) :: extraFirst, cyclic
    integer :: phys(0:nSpod-1), nodeIdx(0:nSpod-1), nodeRank(0:nSpod-1), nodeSize(0:nSpod-1)
    integer :: firstRankOnPhys(0:nNodes-1), idxOfPhys(0:nNodes-1)
    integer :: r, nUsed, ib, nb, nW
    integer(i8) :: i, j, kk, ws(3), wc(3)
    integer, allocatable :: hits(:, :, :)
    type(box_t), allocatable :: boxes(:)

    check = .true.
    nW = prow * pcol

    ! physical node of each SPOD rank
    do r = 0, nSpod - 1
      if (cyclic) then
        phys(r) = mod(r, nNodes)
      else
        phys(r) = writer_node(r, nSpod, nNodes, 0)
      end if
    end do

    ! node index = order of the lowest SPOD rank on each node (leaders comm)
    firstRankOnPhys = -1
    do r = 0, nSpod - 1
      if (firstRankOnPhys(phys(r)) < 0) firstRankOnPhys(phys(r)) = r
    end do
    nUsed = 0
    idxOfPhys = -1
    do r = 0, nSpod - 1
      if (firstRankOnPhys(phys(r)) == r) then
        idxOfPhys(phys(r)) = nUsed
        nUsed = nUsed + 1
      end if
    end do
    do r = 0, nSpod - 1
      nodeIdx(r)  = idxOfPhys(phys(r))
      nodeRank(r) = count(phys(0:r-1) == phys(r))
      nodeSize(r) = count(phys == phys(r))
    end do

    allocate(hits(0:nx-1, 0:ny-1, 0:nz-1))
    hits = 0
    do r = 0, nSpod - 1
      call node_local_boxes(nodeIdx(r), nodeRank(r), nodeSize(r), nUsed, nx, ny, nz, &
                            prow, pcol, pencil, extraFirst, 0, boxes, nb)
      do ib = 1, nb
        call writer_box(boxes(ib)%writer, nx, ny, nz, prow, pcol, pencil, extraFirst, ws, wc)
        if (any(boxes(ib)%start < ws) .or. any(boxes(ib)%start + boxes(ib)%count > ws + wc)) then
          check = .false.          ! box not inside its writer's pencil
        end if
        ! the X3D rank's physical node (sequential placement over all nodes)
        ! must be the physical node of this SPOD rank
        if (writer_node(boxes(ib)%writer, nW, nNodes, 0) /= phys(r)) check = .false.
        do kk = boxes(ib)%start(3), boxes(ib)%start(3) + boxes(ib)%count(3) - 1
          do j = boxes(ib)%start(2), boxes(ib)%start(2) + boxes(ib)%count(2) - 1
            do i = boxes(ib)%start(1), boxes(ib)%start(1) + boxes(ib)%count(1) - 1
              hits(i, j, kk) = hits(i, j, kk) + 1
            end do
          end do
        end do
      end do
      deallocate(boxes)
    end do
    if (any(hits /= 1)) check = .false.     ! every point read exactly once
    deallocate(hits)
  end function check

end program test_node_local_map
