module coarray_aux
  use mpi_f08
  use timer
  implicit none

  public co_init_comms, co_init_distr

  interface co_sum
    module procedure :: co_sum_double
  end interface

  interface co_max
    module procedure :: co_max_double
  end interface

  interface co_gather
    module procedure :: co_gatherv_double
  end interface

  integer,dimension(:),allocatable  :: proc_sizes, proc_offsets
  integer                           :: comm_size, proc_rank
  logical                           :: comms_inited = .false., distr_inited=.false.

contains

  subroutine co_sum_double(x, result_image)
    real*8, intent(inout), dimension(:,:) :: x
    integer, optional :: result_image
    integer :: i
    integer, save :: result_image_mpi[*]

    call TimerStart('co_sum_double')

    if (present(result_image)) then

      if (this_image() .eq. 1) then
        call mpi_comm_rank(mpi_comm_world, result_image_mpi)
        do i = 2, num_images()
         result_image_mpi[i] = result_image_mpi
        end do
      end if
      sync all

      if (this_image() .eq. 1) then
        call mpi_reduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, result_image_mpi, mpi_comm_world)
      else
        call mpi_reduce(x, x, size(x), mpi_double_precision, mpi_sum, result_image_mpi, mpi_comm_world)
      endif
    else
      call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, mpi_comm_world)
    end if
    call TimerStop('co_sum_double')
  end subroutine

  subroutine co_max_double(x, result_image)
    real*8, intent(inout), dimension(:,:) :: x
    integer, optional :: result_image
    integer :: i
    integer, save :: result_image_mpi[*]

    call TimerStart('co_max_double')

    if (present(result_image)) then

      if (this_image() .eq. 1) then
        call mpi_comm_rank(mpi_comm_world, result_image_mpi)
        do i = 2, num_images()
         result_image_mpi[i] = result_image_mpi
        end do
      end if
      sync all

      if (this_image() .eq. 1) then
        call mpi_reduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_max, result_image_mpi, mpi_comm_world)
      else
        call mpi_reduce(x, x, size(x), mpi_double_precision, mpi_max, result_image_mpi, mpi_comm_world)
      endif
    else
      call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_max, mpi_comm_world)
    end if
    call TimerStop('co_max_double')
  end subroutine

  subroutine co_gatherv_double(x)
    real*8, intent(inout), dimension(:,:) :: x
    integer :: ierr

    if (.not. comms_inited .or. .not. distr_inited) stop "COMMS NOT INITIALISED"
    if (comm_size .eq. 1) return

    call TimerStart('CO_GATHERV_DOUBLE')
    if (proc_rank.eq.0) then
      call mpi_gatherv(x, 0, mpi_double_precision, x, proc_sizes, proc_offsets, mpi_double_precision, 0, mpi_comm_world)
    else
      call mpi_gatherv(x, size(x), mpi_double_precision, x, proc_sizes, proc_offsets, mpi_double_precision, 0, mpi_comm_world)
    endif
    if (ierr) stop "co_gatherv_double"
    call TimerStop('CO_GATHERV_DOUBLE')

  end subroutine co_gatherv_double


  subroutine co_init_comms()
    integer :: ierr

    call mpi_init(ierr)
    if (ierr) stop "MPI_INIT"
    call mpi_comm_size(mpi_comm_world, comm_size, ierr)
    if (ierr) stop "MPI_COMM_SIZE"
    call mpi_comm_rank(mpi_comm_world, proc_rank, ierr)
    if (ierr) stop "MPI_COMM_RANK"

    comms_inited = .true.

  end subroutine co_init_comms

  subroutine co_init_distr(dimen, startdim, enddim, blocksize)
    integer,intent(in) :: dimen
    integer,intent(out) :: startdim, enddim, blocksize
    integer,dimension(:),allocatable  :: starts, ends
    integer :: localsize, proc_index
    integer :: i, ierr

    if (.not. comms_inited) stop "COMMS NOT INITIALISED"
    if (distr_inited) stop "DISTRIBUTION ALREADY INITIALISED"

    proc_index = proc_rank+1

    ! While co-arraying
    if (this_image().ne.proc_index) stop "coarray/mpi mixup"

    allocate(proc_sizes(comm_size),proc_offsets(comm_size),starts(comm_size),ends(comm_size),stat=ierr)
    if (ierr) stop "CO_INIT_DISTR ALLOCATION FAILED"

    if (comm_size .eq. 1) then
      startdim = 1
      enddim = dimen
      blocksize = dimen*dimen
    else

      if (proc_rank .eq. 0) then !root

        localsize = dimen/comm_size

        starts(1) = 1
        ends(1) = localsize
        proc_sizes(1) = localsize*dimen
        proc_offsets(1) = 0

        do i=2,comm_size-1
          starts(i) = (i-1)*localsize+1
          ends(i) = i*localsize
          proc_sizes(i) = localsize * dimen
          proc_offsets(i) = localsize*(i-1)*dimen
        end do

        starts(comm_size) = (i-1) * localsize + 1
        ends(comm_size) = dimen
        proc_sizes(comm_size) = (dimen-localsize*(comm_size-1)) * dimen

        proc_offsets(comm_size) = (comm_size-1)*localsize*dimen
      endif

      call mpi_bcast(starts, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(ends, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_sizes, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_offsets, comm_size, mpi_integer, 0, mpi_comm_world)

      blocksize = proc_sizes(proc_index)
      startdim = starts(proc_index)
      enddim = ends(proc_index)

    endif

    deallocate(starts,ends)

    distr_inited = .true.
  end subroutine co_init_distr


end module
