subroutine to_file(input_filename, i,j)
    use mod_streams
    use mod_sys
    implicit none
    character(len=50) :: input_filename

    integer :: i, j, k
    real(mykind) :: rho, rhou, rhov, rhow, vwrite, uwrite, wwrite

    open(21,file=input_filename)
    write(21, *) "rho, u, v, w"

    do k = 1, nz
        rho = w(i,j,k,1)
        rhou = w(i,j,k,2)
        rhov = w(i,j,k,3)
        rhow = w(i,j,k,4)

        uwrite = rhou / rho
        vwrite = rhov / rho
        wwrite = rhow / rho

        write(21, "(E15.10, A1, E15.10, A1, E15.10, A1, E15.10)") rho, ",", uwrite, ",", vwrite, ",", wwrite
    enddo

    close(21)
end subroutine to_file 

subroutine write_probe_data()
    use mod_streams
    use mod_sys

    implicit none
    integer :: i, j, k, dead
    character(len=50) :: filename

    ! if we are the main thread
    if (masterproc) then
        ! write probe data 1/4th in the x direction and half in y direction
        write(filename, "(A22, I5.5, A4)") "csv_data/probe_data_1_", icyc, ".csv"
        i = nx *1/4
        j = ny / 2
        !print *, filename
        call to_file(filename, i, j)

        ! write probe data 1/2 in the x direction and half in y direction
        write(filename, "(A22, I5.5, A4)") "csv_data/probe_data_2_", icyc, ".csv"
        i = nx*2/4
        !print *, filename
        call to_file(filename, i, j)

        ! write probe data 3/4 in the x direction and half in y direction
        write(filename, "(A22, I5.5, A4)") "csv_data/probe_data_3_", icyc, ".csv"
        i = nx*3/4
        !print *, filename
        call to_file(filename, i, j)
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
