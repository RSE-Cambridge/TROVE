module coarray_aux
  use mpi_f08
  implicit none

  interface co_sum
    module procedure :: co_sum_double
  end interface

contains

  subroutine co_sum_double(x, result_image)
    real*8, intent(inout), dimension(:,:) :: x
    integer, optional :: result_image
    integer :: i
    integer, save :: result_image_mpi[*]

    if (present(result_image)) then

      if (this_image() .eq. 1) then
        call mpi_comm_rank(mpi_comm_world, result_image_mpi)
        do i = 2, num_images()
         result_image_mpi[i] = result_image_mpi
        end do
      end if
      sync all

      call mpi_reduce(x, mpi_in_place, size(x), mpi_double_precision, mpi_sum, result_image_mpi, mpi_comm_world)
    else
      call mpi_allreduce(x, mpi_in_place, size(x), mpi_double_precision, mpi_sum, mpi_comm_world)
    end if
  end subroutine

end module
