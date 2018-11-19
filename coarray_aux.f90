module coarray_aux
  use mpi_f08
  use timer
  implicit none

  interface co_sum
    module procedure :: co_sum_double
  end interface

  interface co_max
    module procedure :: co_max_double
  end interface

contains

  subroutine co_sum_double(x, result_image)
    real*8, intent(inout), dimension(:,:) :: x
    integer, optional :: result_image
    integer :: i
    integer, save :: result_image_mpi[*]

    call TimerStart('co_sum_double')
    !write(*,*) "CO_SUM_DOUBLE ENTERED @", this_image()

    if (present(result_image)) then

      !write(*,*) "CO_SUM_DOUBLE/RESULT_IMAGE PRESENT"

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
      !write(*,*) "CO_SUM_DOUBLE/RESULT_IMAGE <NOT> PRESENT"
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
    !write(*,*) "CO_SUM_DOUBLE ENTERED @", this_image()

    if (present(result_image)) then

      !write(*,*) "CO_SUM_DOUBLE/RESULT_IMAGE PRESENT"

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
      !write(*,*) "CO_SUM_DOUBLE/RESULT_IMAGE <NOT> PRESENT"
      call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_max, mpi_comm_world)
    end if
    call TimerStop('co_max_double')
  end subroutine

end module
