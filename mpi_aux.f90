module mpi_aux
  use mpi_f08
  use timer
  use accuracy
  implicit none

  public co_init_comms, co_finalize_comms, co_init_distr, co_distr_data, co_write_matrix_distr, co_read_matrix_distr
  public co_create_type

  public send_or_recv, comm_size, mpi_rank
  public co_startdim, co_enddim

  public blacs_size, blacs_rank, blacs_ctxt
  public nprow,npcol,myprow,mypcol, desca,descb,descc

  interface co_sum
    module procedure :: co_sum_double
  end interface

  !interface co_max
  !  module procedure :: co_max_double
  !end interface

  interface co_gather
    module procedure :: co_gather_double
    module procedure :: co_gatherv_double
  end interface

  integer,dimension(:),allocatable  :: proc_sizes, proc_offsets, send_or_recv
  integer                           :: comm_size, mpi_rank
  integer                           :: co_startdim, co_enddim
  logical                           :: comms_inited = .false., distr_inited=.false.
  type(MPI_Datatype) :: mpitype_column
  type(MPI_Datatype),dimension(:), allocatable :: mpi_blocktype

  !blacs/pblas
  integer :: blacs_size, blacs_rank, blacs_ctxt
  integer :: nprow,npcol,myprow,mypcol
  integer :: desca(9)
  integer :: descb(9)
  integer :: descc(9)
  integer :: descd(9)

contains

  subroutine co_init_pblas()
    call blacs_pinfo(blacs_rank, blacs_size)
    if (blacs_rank .lt. 0) return

    call blacs_get(-1, 0, blacs_ctxt)
    call blacs_gridinit(blacs_ctxt, 'R', blacs_size/8, 8)
    call blacs_gridinfo(blacs_ctxt, nprow, npcol, myprow, mypcol)

    write(*,"('BLACS: [',i2,',',i2'](',i4,i4,i4,i4',)')") mpi_rank,blacs_rank,nprow,npcol,myprow,mypcol
  end subroutine co_init_pblas

  subroutine co_sum_double(x, result_image)
    real*8, intent(inout), dimension(:,:) :: x
    integer, optional :: result_image
    integer :: i
    !integer, save :: result_image_mpi[*]

    if (comm_size.eq.1) return
    call TimerStart('co_sum_double')

    !if (present(result_image)) then

      !if (this_image() .eq. 1) then
      !  call mpi_comm_rank(mpi_comm_world, result_image_mpi)
      !  do i = 2, num_images()
      !   result_image_mpi[i] = result_image_mpi
      !  end do
      !end if
      !sync all

      if (mpi_rank .eq. 0) then
        call mpi_reduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world)
      else
        call mpi_reduce(x, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world)
      endif
    !else
    !  call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, mpi_comm_world)
    !end if
    call TimerStop('co_sum_double')
  end subroutine

  !subroutine co_max_double(x, result_image)
  !  real*8, intent(inout), dimension(:,:) :: x
  !  integer, optional :: result_image
  !  integer :: i
  !  integer, save :: result_image_mpi[*]

  !  call TimerStart('co_max_double')

  !  if (present(result_image)) then

  !    if (this_image() .eq. 1) then
  !      call mpi_comm_rank(mpi_comm_world, result_image_mpi)
  !      do i = 2, num_images()
  !       result_image_mpi[i] = result_image_mpi
  !      end do
  !    end if
  !    sync all

  !    if (this_image() .eq. 1) then
  !      call mpi_reduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_max, result_image_mpi, mpi_comm_world)
  !    else
  !      call mpi_reduce(x, x, size(x), mpi_double_precision, mpi_max, result_image_mpi, mpi_comm_world)
  !    endif
  !  else
  !    call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_max, mpi_comm_world)
  !  end if
  !  call TimerStop('co_max_double')
  !end subroutine

  subroutine co_gather_double(x, static)
    real*8, intent(inout), dimension(:,:) :: x
    logical, intent(in) :: static
    integer :: ierr

    if (.not. comms_inited .or. .not. distr_inited) stop "COMMS NOT INITIALISED"
    if (comm_size .eq. 1) return

    call TimerStart('CO_GATHER_DOUBLE')
    call mpi_gather(x, 0, mpi_double_precision, x, proc_sizes(2), mpi_double_precision, 0, mpi_comm_world)
    if (ierr .gt. 0) stop "co_gather_double"
    call TimerStop('CO_GATHER_DOUBLE')

  end subroutine co_gather_double

  subroutine co_gatherv_double(x)
    real*8, intent(inout), dimension(:,:) :: x
    integer :: ierr

    if (.not. comms_inited .or. .not. distr_inited) stop "COMMS NOT INITIALISED"
    if (comm_size .eq. 1) return

    call TimerStart('CO_GATHERV_DOUBLE')
    if (mpi_rank.eq.0) then
      call mpi_gatherv(x, 0, mpi_double_precision, x, proc_sizes, proc_offsets, mpi_double_precision, 0, mpi_comm_world)
    else
      call mpi_gatherv(x, size(x), mpi_double_precision, x, proc_sizes, proc_offsets, mpi_double_precision, 0, mpi_comm_world)
    endif
    if (ierr .gt. 0) stop "co_gatherv_double"
    call TimerStop('CO_GATHERV_DOUBLE')

  end subroutine co_gatherv_double

  subroutine co_init_comms()
    integer :: ierr

    call mpi_init(ierr)
    if (ierr .gt. 0) stop "MPI_INIT"
    call mpi_comm_size(mpi_comm_world, comm_size, ierr)
    if (ierr .gt. 0) stop "MPI_COMM_SIZE"
    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)
    if (ierr .gt. 0) stop "MPI_COMM_RANK"

    comms_inited = .true.

  end subroutine co_init_comms

  subroutine co_finalize_comms()
    integer :: ierr

    if (.not. comms_inited) stop "CO_FINALIZE_COMMS COMMS NOT INITED"

    call mpi_finalize(ierr)

    if (ierr .gt. 0) stop "MPI_FINALIZE"

  end subroutine co_finalize_comms

  subroutine co_init_distr(dimen, startdim, enddim, blocksize)
    integer,intent(in) :: dimen
    integer,intent(out) :: startdim, enddim, blocksize
    integer,dimension(:),allocatable  :: starts, ends
    integer :: localsize, proc_index, localsize_
    integer :: i, ierr, to_calc

    if (.not. comms_inited) stop "COMMS NOT INITIALISED"
    if (distr_inited) stop "DISTRIBUTION ALREADY INITIALISED"

    proc_index = mpi_rank+1

    allocate(proc_sizes(comm_size),proc_offsets(comm_size),send_or_recv(comm_size),starts(comm_size),ends(comm_size),stat=ierr)
    if (ierr .gt. 0) stop "CO_INIT_DISTR ALLOCATION FAILED"

    if (comm_size .eq. 1) then
      startdim = 1
      enddim = dimen
      blocksize = dimen*dimen
      send_or_recv(1) = 0
    else

      if (mpi_rank .eq. 0) then !root

        localsize = dimen/comm_size
        localsize_ = int(1+real(dimen/comm_size))

        starts(1) = 1
        ends(1) = localsize_
        proc_sizes(1) = localsize_*(comm_size*localsize_)
        proc_offsets(1) = 0

        do i=2,comm_size-1
          starts(i) = (i-1)*localsize_+1
          ends(i) = i*localsize_
          proc_sizes(i) = localsize_ * (comm_size*localsize_)!dimen
          proc_offsets(i) = localsize_*(i-1)*(comm_size*localsize_)!dimen
        end do

        starts(comm_size) = (i-1) * localsize_ + 1
        ends(comm_size) = dimen!comm_size*localsize_!dimen
        proc_sizes(comm_size) = localsize_*comm_size*localsize_!dimen

        proc_offsets(comm_size) = (comm_size-1)*localsize_*(comm_size*localsize_)!dimen
      endif

      call mpi_bcast(starts, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(ends, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_sizes, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_offsets, comm_size, mpi_integer, 0, mpi_comm_world)



      blocksize = proc_sizes(proc_index)
      startdim = starts(proc_index)
      enddim = ends(proc_index)

      co_startdim = startdim
      co_enddim = enddim

      allocate(mpi_blocktype(comm_size))

      do i=1,comm_size
        if (mod(comm_size,2).eq.1) then
          to_calc = comm_size/2+1
        else
          if ((mod(comm_size, 4).eq.0) .and. ((mod(i,2).eq.1 .and. i.le.comm_size/2).or.(mod(i,2).eq.0 .and. i.gt.comm_size/2))) then
            to_calc = comm_size/2+1
          elseif (mod(i,2).eq.1 .and. (mod(comm_size, 4).gt.0)) then
            to_calc = comm_size/2+1
          else
            to_calc = comm_size/2
          endif
        endif


        if (i.eq.proc_index) then
          send_or_recv(i) = 0
        elseif ( ((i.gt.(proc_index - to_calc) .and. i.lt.proc_index)) .or. &
            ((proc_index-to_calc).lt.1 .and. (i-comm_size).gt.(proc_index-to_calc))) then
          send_or_recv(i) = 1 ! send
          call co_create_type_subarray(int(1+real(dimen/comm_size)), blocksize, int(1+real(dimen/comm_size)), i, mpi_blocktype(i))
        else
          send_or_recv(i) = -1 ! recv
        endif
      end do

    endif

    call co_create_type(dimen)
    deallocate(starts,ends)

    !call co_init_pblas()

    distr_inited = .true.
  end subroutine co_init_distr

  subroutine co_distr_data(x, tmp, blocksize, lb, ub)
    use accuracy

    implicit none

    real(rk),dimension(:,lb:),intent(inout) :: x
    real(rk),dimension(:,:,:),intent(inout) :: tmp
    integer,intent(in)                :: blocksize, lb, ub

    integer :: i, icoeff, jcoeff, offset, ierr, k
    type(MPI_Request)  :: reqs(comm_size)

    call TimerStart('MPI_transpose')
    call TimerStart('MPI_transpose_sendrecv')

    do i=1,comm_size
      reqs(i)= MPI_REQUEST_NULL
    end do

    do i=1,comm_size
      if (send_or_recv(i).eq.1) then
        call mpi_isend(x,1,mpi_blocktype(i),i-1,0,mpi_comm_world,reqs(i),ierr)
      elseif (send_or_recv(i).eq.-1) then
        call mpi_irecv(tmp(:,:,i),blocksize*blocksize,mpi_double_precision,i-1,mpi_any_tag,mpi_comm_world,reqs(i),ierr)
      else
        reqs(i) = MPI_REQUEST_NULL
      endif
    enddo

    call mpi_waitall(comm_size,reqs,mpi_statuses_ignore,ierr)
    call TimerStop('MPI_transpose_sendrecv')
    call TimerStart('MPI_transpose_local')

    do i=1,comm_size
      if (send_or_recv(i).eq.-1) then
        offset = (i-1)*blocksize
        !$omp parallel do private(icoeff,jcoeff) shared(i,x,tmp,lb,ub,offset,blocksize) schedule(static)
        do icoeff=lb,ub
          do jcoeff=offset+1,offset+blocksize
            x(jcoeff,icoeff) = tmp(icoeff-lb+1,jcoeff-offset,i)
          enddo
        enddo
        !$omp end parallel do
      endif
    enddo
    call TimerStop('MPI_transpose_local')
    call TimerStop('MPI_transpose')

  end subroutine co_distr_data

  subroutine co_read_matrix_distr(x, longdim, lb, ub, infile)
    use accuracy

    real(rk),dimension(:,lb:),intent(in) :: x
    integer,intent(in)                :: longdim, lb, ub

    type(MPI_File),intent(in) :: infile
    type(MPI_Status) :: writestat
    integer(kind=MPI_Offset_kind) :: offset_start,offset_end
    integer :: readcount, mpi_real_size, ierr

    call MPI_Type_size(mpi_double_precision, mpi_real_size,ierr)

    if (mpi_rank .lt. (comm_size-1)) then
      readcount = int(1+real(longdim/comm_size))
    else
      readcount = longdim-((comm_size-1)*int(1+real(longdim/comm_size)))
    endif

    offset_start = mpi_rank * (longdim * int(1+real(longdim/comm_size),mpi_offset_kind) * mpi_real_size)
    offset_end = longdim
    offset_end = (offset_end * offset_end * mpi_real_size) - offset_start

    call MPI_File_seek(infile, offset_start, MPI_SEEK_CUR)
    call MPI_File_read_all(infile,x,readcount,mpitype_column,writestat,ierr)
    call MPI_File_seek(infile, offset_end, MPI_SEEK_CUR)

  end subroutine co_read_matrix_distr

  subroutine co_write_matrix_distr(x, longdim, lb, ub, outfile)
    use accuracy

    real(rk),dimension(:,lb:),intent(in) :: x
    integer,intent(in)                :: longdim, lb, ub
    type(MPI_File),intent(in) :: outfile
    integer :: ierr, mpi_real_size, writecount, mpi_col_size
    !integer(kind=MPI_Offset_kind) :: mpioffset,mpi_write_offsetkind
    integer(kind=MPI_Offset_kind) :: offset_start, offset_end
    type(MPI_Status) :: writestat

    call mpi_barrier(mpi_comm_world, ierr)

    call TimerStart('MPI_write')

    call MPI_Type_size(mpi_double_precision, mpi_real_size,ierr)
    call MPI_Type_size(mpitype_column, mpi_col_size,ierr)

    if (mpi_rank .lt. (comm_size-1)) then
      writecount = int(1+real(longdim/comm_size))
    else
      writecount = longdim-((comm_size-1)*int(1+real(longdim/comm_size)))
    endif

    offset_start = mpi_rank * int(1+real(longdim/comm_size),mpi_offset_kind) * mpi_col_size
    offset_end = 0

    call MPI_File_seek(outfile, offset_start, MPI_SEEK_END)
    call MPI_File_write_all(outfile,x,writecount,mpitype_column,writestat,ierr)
    call MPI_File_seek(outfile, offset_end, MPI_SEEK_END)
    call TimerStop('MPI_write')

  end subroutine co_write_matrix_distr

  subroutine co_create_type(extent)
    integer, intent(in) :: extent
    integer :: ierr

    call MPI_Type_contiguous(extent, mpi_double_precision, mpitype_column, ierr)
    call MPI_Type_commit(mpitype_column, ierr)

  end subroutine co_create_type

  subroutine co_create_type_subarray(extent, coldim, rowdim, blockid, mpi_newtype)
    integer,intent(in) :: extent, coldim, rowdim, blockid
    type(MPI_Datatype),intent(inout) :: mpi_newtype
    integer,dimension(2) :: array_of_sizes, array_of_subsizes, array_of_starts
    integer :: ierr

    array_of_sizes(1) = comm_size * extent!coldim
    array_of_sizes(2) = extent
    array_of_subsizes(:) = extent
    array_of_starts(1) = (blockid - 1) * extent + 0
    array_of_starts(2) = 0


    call MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, mpi_double_precision, mpi_newtype, ierr)
    call MPI_Type_commit(mpi_newtype, ierr)

  end subroutine co_create_type_subarray

end module
