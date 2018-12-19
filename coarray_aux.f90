module coarray_aux
  use mpi_f08
  use timer
  implicit none

  public co_init_comms, co_finalize_comms, co_init_distr, co_distr_data

  public send_or_recv, comm_size, proc_rank

  !interface co_sum
  !  module procedure :: co_sum_double
  !end interface

  !interface co_max
  !  module procedure :: co_max_double
  !end interface

  interface co_gather
    module procedure :: co_gather_double
    module procedure :: co_gatherv_double
  end interface

  integer,dimension(:),allocatable  :: proc_sizes, proc_offsets, send_or_recv
  integer                           :: comm_size, proc_rank
  logical                           :: comms_inited = .false., distr_inited=.false.

contains

  !subroutine co_sum_double(x, result_image)
  !  real*8, intent(inout), dimension(:,:) :: x
  !  integer, optional :: result_image
  !  integer :: i
  !  integer, save :: result_image_mpi[*]

  !  call TimerStart('co_sum_double')

  !  if (present(result_image)) then

  !    if (this_image() .eq. 1) then
  !      call mpi_comm_rank(mpi_comm_world, result_image_mpi)
  !      do i = 2, num_images()
  !       result_image_mpi[i] = result_image_mpi
  !      end do
  !    end if
  !    sync all

  !    if (this_image() .eq. 1) then
  !      call mpi_reduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, result_image_mpi, mpi_comm_world)
  !    else
  !      call mpi_reduce(x, x, size(x), mpi_double_precision, mpi_sum, result_image_mpi, mpi_comm_world)
  !    endif
  !  else
  !    call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, mpi_comm_world)
  !  end if
  !  call TimerStop('co_sum_double')
  !end subroutine

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
    if (ierr) stop "co_gather_double"
    call TimerStop('CO_GATHER_DOUBLE')

  end subroutine co_gather_double

  subroutine co_gatherv_double(x)
    real*8, intent(inout), dimension(:,:) :: x
    integer :: ierr

    if (.not. comms_inited .or. .not. distr_inited) stop "COMMS NOT INITIALISED"
    if (comm_size .eq. 1) return

    !write(*,*) proc_rank, proc_sizes

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

  subroutine co_finalize_comms()
    integer :: ierr

    if (.not. comms_inited) stop "CO_FINALIZE_COMMS COMMS NOT INITED"

    call mpi_finalize(ierr)

    if (ierr) stop "MPI_FINALIZE"
  end subroutine co_finalize_comms

  subroutine co_init_distr(dimen, startdim, enddim, blocksize)
    integer,intent(in) :: dimen
    integer,intent(out) :: startdim, enddim, blocksize
    integer,dimension(:),allocatable  :: starts, ends
    integer :: localsize, proc_index, localsize_
    integer :: i, ierr, to_calc

    if (.not. comms_inited) stop "COMMS NOT INITIALISED"
    if (distr_inited) stop "DISTRIBUTION ALREADY INITIALISED"

    proc_index = proc_rank+1

    ! While co-arraying
    !if (this_image().ne.proc_index) stop "coarray/mpi mixup"

    allocate(proc_sizes(comm_size),proc_offsets(comm_size),send_or_recv(comm_size),starts(comm_size),ends(comm_size),stat=ierr)
    if (ierr) stop "CO_INIT_DISTR ALLOCATION FAILED"

    if (comm_size .eq. 1) then
      startdim = 1
      enddim = dimen
      blocksize = dimen*dimen
      send_or_recv(1) = 0
    else

      if (proc_rank .eq. 0) then !root

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
        !proc_sizes(comm_size) = (dimen-localsize*(comm_size-1)) * dimen
        proc_sizes(comm_size) = localsize_*comm_size*localsize_!dimen

        proc_offsets(comm_size) = (comm_size-1)*localsize_*(comm_size*localsize_)!dimen
      endif

      call mpi_bcast(starts, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(ends, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_sizes, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_offsets, comm_size, mpi_integer, 0, mpi_comm_world)


      !if (proc_rank.eq.0) write(*,*) "PROC_SIZES:", proc_sizes

      blocksize = proc_sizes(proc_index)
      startdim = starts(proc_index)
      enddim = ends(proc_index)


      !if (mod(comm_size,2)) then
      !  do i=1,comm_size
      !    if (i.eq.proc_index) then
      !      send_or_recv(i) = 0
      !    !elseif ( (i .ge. (proc_index - to_calc) .and. i .lt. proc_index) ) .or. &
      !    elseif ( ((i.gt.(proc_index - to_calc) .and. i.lt.proc_index)) .or. &
      !        ((proc_index-to_calc).lt.1 .and. (i-comm_size).gt.(proc_index-to_calc))) then
      !      send_or_recv(i) = 1 ! send
      !    else
      !      send_or_recv(i) = -1 ! recv
      !    endif
      !  end do
      !else
      !  do i=1,comm_size
      !    if (mod(i,2)) then
      !      to_calc = comm_size/2+1
      !    else
      !      to_calc = comm_size/2
      !    endif

      !    if (i.eq.proc_index) then
      !      send_or_recv(i) = 0
      !    !elseif ( (i .ge. (proc_index - to_calc) .and. i .lt. proc_index) ) .or. &
      !    elseif ( ((i.gt.(proc_index - to_calc) .and. i.lt.proc_index)) .or. &
      !        ((proc_index-to_calc).lt.1 .and. (i-comm_size).gt.(proc_index-to_calc))) then
      !      send_or_recv(i) = 1 ! send
      !    else
      !      send_or_recv(i) = -1 ! recv
      !    endif
      !  end do
      !endif

      do i=1,comm_size
        if (mod(comm_size,2)) then
          to_calc = comm_size/2+1
        else
          if ((mod(comm_size, 4).eq.0) .and. ((mod(i,2) .and. i.le.comm_size/2).or.(.not.mod(i,2) .and. i.gt.comm_size/2))) then
            to_calc = comm_size/2+1
          elseif (mod(i,2) .and. (mod(comm_size, 4).gt.0)) then
            to_calc = comm_size/2+1
          else
            to_calc = comm_size/2
          endif
        endif

        if (proc_rank .eq. 0) write(*,*) "TOCALC:", i, to_calc, comm_size, mod(comm_size,4)

        if (i.eq.proc_index) then
          send_or_recv(i) = 0
        !!!!!!!!elseif ( (i .ge. (proc_index - to_calc) .and. i .lt. proc_index) ) .or. &
        elseif ( ((i.gt.(proc_index - to_calc) .and. i.lt.proc_index)) .or. &
            ((proc_index-to_calc).lt.1 .and. (i-comm_size).gt.(proc_index-to_calc))) then
          send_or_recv(i) = 1 ! send
        !elseif (((i.gt.proc_index .and. i.lt.(proc_index+to_calc))) .or. (proc_index+to_calc.gt.comm_size .and. i.lt.mod(proc_index+to_calc,comm_size))) then
          send_or_recv(i) = 1 ! send
        else
          send_or_recv(i) = -1 ! recv
        endif
      end do

      !write(*,*) "SENDRECV:", proc_index, send_or_recv
          

    endif

    deallocate(starts,ends)

    distr_inited = .true.
  end subroutine co_init_distr

  subroutine co_distr_data(x, tmp, blocksize, lb, ub)
    use accuracy


    real(rk),dimension(:,lb:),intent(inout) :: x
    real(rk),dimension(:,:,:),intent(inout) :: tmp
    integer,intent(in)                :: blocksize, lb, ub

    integer :: i, icoeff, jcoeff, offset, ierr
    type(MPI_Request)  :: reqs(comm_size)


    !!!write(*,*) "DISTR1", proc_rank, send_or_recv
    !!!write(*,*) "DISTR2", proc_rank, blocksize, lb, ub
    !!!write(*,*) "DISTR3", proc_rank, shape(x), shape(tmp)

    do i=1,comm_size
      if (send_or_recv(i).eq.1) then
        call mpi_isend(x(((i-1)*blocksize)+1:i*blocksize,:),blocksize*blocksize,mpi_double_precision,i-1,0,mpi_comm_world,reqs(i),ierr)
      elseif (send_or_recv(i).eq.-1) then
        call mpi_irecv(tmp(:,:,i),blocksize*blocksize,mpi_double_precision,i-1,mpi_any_tag,mpi_comm_world,reqs(i),ierr)
      else
        reqs(i) = MPI_REQUEST_NULL
      endif
    enddo

    call mpi_waitall(comm_size,reqs,mpi_statuses_ignore,ierr)

    do i=1,comm_size
      if (send_or_recv(i).eq.-1) then
        offset = (i-1)*blocksize
        !$omp parallel do private(icoeff,jcoeff) shared(i,x) schedule(static)
        do icoeff=lb,ub
          do jcoeff=offset+1,offset+blocksize
            x(jcoeff,icoeff) = tmp(icoeff-lb+1,jcoeff-offset,i)
          enddo
        enddo
        !$omp end parallel do
      endif
    enddo

  end subroutine co_distr_data

end module
