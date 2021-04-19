subroutine to_file(input_filename, input_filename_2, i,j)
    use mod_streams
    use mod_sys
    implicit none
    character(len=70) :: input_filename , input_filename_2

    integer :: i, j, k
    real(mykind) :: rho, rhou, rhov, rhow, vwrite, uwrite, wwrite

    open(21,file=input_filename)
    open(22,file=input_filename_2)
    write(21, *) "rho, u, v, w"

    do k = 1, nz
        rho =  w(i,j,k,1)
        rhou = w(i,j,k,2)
        rhov = w(i,j,k,3)
        rhow = w(i,j,k,4)

        uwrite = rhou / rho
        vwrite = rhov / rho
        wwrite = rhow / rho

        write(21, "(E25.20, A1, E25.20, A1, E25.20, A1, E25.20)") rho, ",", uwrite, ",", vwrite, ",", wwrite

        rho =  w_gpu(i,j,k,1)
        rhou = w_gpu(i,j,k,2)
        rhov = w_gpu(i,j,k,3)
        rhow = w_gpu(i,j,k,4)

        uwrite = rhou / rho
        vwrite = rhov / rho
        wwrite = rhow / rho
        write(22, "(E25.20, A1, E25.20, A1, E25.20, A1, E25.20)") rho, ",", uwrite, ",", vwrite, ",", wwrite
    enddo

    close(21)
end subroutine to_file

subroutine write_probe_data()
    use mod_streams
    use mod_sys

    implicit none
    integer :: i, j, k, dead
    character(len=70) :: filename, filename2

    ! if we are the main thread
    if (masterproc) then
        ! write probe data 1/4th in the x direction and half in y direction
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_1_", icyc, ".csv"
        write(filename2, "(A28, I5.5, A4)") "csv_data/w_gpu_probe_data_1_", icyc, ".csv"
        i = nx *1/4
        j = ny / 2
        !print *, filename
        call to_file(filename, filename2, i, j)

        ! write probe data 1/2 in the x direction and half in y direction
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_2_", icyc, ".csv"
        write(filename2, "(A28, I5.5, A4)") "csv_data/w_gpu_probe_data_2_", icyc, ".csv"
        i = nx*2/4
        !print *, filename
        call to_file(filename, filename2, i, j)

        ! write probe data 3/4 in the x direction and half in y direction
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_3_", icyc, ".csv"
        write(filename2, "(A28, I5.5, A4)") "csv_data/w_gpu_probe_data_3_", icyc, ".csv"
        i = nx*3/4
        !print *, filename
        call to_file(filename, filename2, i, j)
    endif

end subroutine write_probe_data


subroutine init_write_telaps()
    use mod_streams
    implicit none

    if (masterproc) then
        open(995, file="timesteps.csv", action="write", status="replace")
        write(995, *) "telaps"
    endif

end subroutine init_write_telaps

subroutine write_telaps(telaps_in)
    use mod_streams
    use mod_sys
    implicit none

    real(mykind) :: telaps_in

    if (masterproc) then
        open(995, file="timesteps.csv", action="write", position="append")
        write(995, "(E15.10, A1)") telaps_in, ","
    endif
end subroutine write_telaps

! average out all of the values in a span for the vector of conservative variables w
subroutine write_span_averaged
    use mod_streams
    use mod_sys

    ! local variables
    character(len=5) :: current_cycle
    character(len=70) :: current_filename

    if (masterproc) then
        ! works for icyc < 99,999 iterations
        write(current_cycle, "(I5)") icyc

        current_filename = "span_average_"//current_cycle//"_rho.dat"
        call helper_write_span_averaged(current_filename, 1)
    endif

end subroutine write_span_averaged

! A helper to average out the values on each of the spans for different variables
! Called from the write_span_averaged subroutine with the name of the file to write to
! and the slice of the variable that we are averaging
subroutine helper_write_span_averaged(filename, slice_var)
    use mod_streams
    use mod_sys

    ! arguments
    integer:: slice_var
    character(len=70) :: filename

    ! local variables
    real(8), dimension(:,:), allocatable :: span_average
    real(8) :: current_average

    do i = 1,nx
        do j = 1, ny
            current_average = 0.0

            ! sum up all of the values on the z axis
            do k = 1, nz
                ! hopefully the compiler can inline each of these branches in the hot loop
                ! otherwise this is painfully slow
                if (slice_var == 1) then
                    ! if we are dealing with rho
                    current_average = current_average + w(i,j,k, slice_var)
                else
                    ! if we are not dealing with rho then we need to divide by it
                    current_average = current_average + (w(i,j,k, slice_var)/w(i,j,k, 1))
                end if

            enddo
            ! calculate the mean of the data
            ! TODO: hopefully this is not an overflow somewhere
            current_average = current_average / nz

            span_average(i,j) = current_average
        enddo
    enddo

    ! write everthing to files and cleanup
    open(23, file=filename)
        write(23, *) slice
    close(23)
end subroutine helper_write_span_averaged
