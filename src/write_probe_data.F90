module mod_probe 
    real(8), dimension(:,:), allocatable :: span_average_rho, span_average_u, span_average_v, span_average_w, span_average
    character(len=5) :: nxstr, nystr
end module mod_probe

subroutine to_file(input_filename, i,j)
    use mod_streams
    use mod_sys
    implicit none
    character(len=70) :: input_filename

    integer :: i, j, k
    real(mykind) :: rho, rhou, rhov, rhow, vwrite, uwrite, wwrite

    open(21,file=input_filename)
    write(21, *) "rho, u, v, w"

    do k = 1, nz
        rho =  w(1,i,j,k)
        rhou = w(2,i,j,k)
        rhov = w(3,i,j,k)
        rhow = w(4,i,j,k)

        uwrite = rhou / rho
        vwrite = rhov / rho
        wwrite = rhow / rho

        write(21, "(E36.30, A1, E36.30, A1, E36.30, A1, E36.30)") rho, ",", uwrite, ",", vwrite, ",", wwrite
    enddo

    close(21)
end subroutine to_file

subroutine write_probe_data()
    use mod_streams
    use mod_sys

    implicit none
    integer :: i, j, k, dead
    character(len=70) :: filename

    ! The three probe locations 
    !  |____________________|____________________|___________________|____________________|
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |          X         |                   X|                   |         X          |    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |____________________|____________________|___________________|____________________|
    !
    !

    ! if we are the main thread
    if (nrank==0) then
        ! write probe data 1/4th in the x direction and half in y direction
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_1_", icyc, ".csv"
        i = nx *1/4
        j = ny / 2
        !print *, filename
        call to_file(filename, i, j)

    elseif (nrank==1) then
        ! write probe data 1/2 in the x direction and half in y direction
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_2_", icyc, ".csv"
        i = nx-1
        !print *, filename
        call to_file(filename, i, j)

    elseif (nrank ==3) then
        ! write probe data in the middle of the fourth process area, which
        ! comes out to be ~3/4ths of the entire x distance
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_3_", icyc, ".csv"
        i = nx*2/4
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

! average out all of the values in a span for the vector of conservative variables w
subroutine write_span_averaged
    use mod_streams
    use mod_sys
    use mod_probe

    ! local variables
    character(len=5) :: current_cycle
    character(len=2) :: mpi_process_number
    character(len=70) :: current_filename

    write(mpi_process_number, "(I2.2)") nrank

    ! works for icyc < 99,999 iterations
    write(current_cycle, "(I5.5)") icyc

    current_filename = "spans/span_average_"//mpi_process_number//"_"//current_cycle//"_rho.vtk"

    allocate(span_average(nx, ny))
        call helper_average_span(1)
    call move_alloc(span_average, span_average_rho)

    allocate(span_average(nx, ny))
        call helper_average_span(2)
    call move_alloc(span_average, span_average_u)

    allocate(span_average(nx, ny))
        call helper_average_span(3)
    call move_alloc(span_average, span_average_v)

    allocate(span_average(nx, ny))
        call helper_average_span(4)
    call move_alloc(span_average, span_average_w)

    call write_span_vtk(current_filename)

    ! everything else is deallocated by moving the allocation to span_average 
    ! when writing to vtk
    ! we deallocate here so that we dont error when allocating at the start of 
    ! this subroutine on the second time that it is called
    deallocate(span_average  )

end subroutine write_span_averaged

! A helper to average out the values on each of the spans for different variables
! Called from the write_span_averaged subroutine with the name of the file to write to
! and the slice of the variable that we are averaging
subroutine helper_average_span(slice_var)
    use mod_streams
    use mod_sys
    use mod_probe
    implicit none

    ! arguments
    integer:: slice_var

    ! local variables
    real(8) :: current_average
    integer:: i, j, k

    do i = 1,nx
        do j = 1, ny
            current_average = 0.0

            ! sum up all of the values on the z axis
            do k = 1, nz
                ! hopefully the compiler can inline each of these branches in the hot loop
                ! otherwise this is painfully slow
                if (slice_var == 1) then
                    ! if we are dealing with rho
                    current_average = current_average + w(slice_var,i,j,k)
                else
                    ! if we are not dealing with rho then we need to divide by it
                    current_average = current_average + (w(slice_var,i,j,k) / w(1,i,j,k))
                end if

            enddo
            ! calculate the mean of the data
            ! TODO: hopefully this is not an overflow somewhere
            current_average = current_average / nz

            span_average(i,j) = current_average
        enddo
    enddo

end subroutine helper_average_span

subroutine write_span_vtk(filename)
    use mod_streams
    use mod_probe
    implicit none

    ! the number of points in each of the directions
    character(len=1000) :: xml ! this allocation errors with nvhpc compiler
    character(len=16) :: curr_cycle
    character(len=70) :: filename

    write(nxstr, "(I4)")  nxmax
    write(nystr, "(I4)")  nymax

    open(23, file=filename, form='formatted')

    ! start by writing the general header that we use for everything
    call write_vtk_header(xml)

    call write_vtk_dataarray(xml, 1)
    call write_vtk_dataarray(xml, 2)
    call write_vtk_dataarray(xml, 3)
    call write_vtk_dataarray(xml, 4)

    call write_vtk_footer(xml)

    close(23)
end subroutine write_span_vtk

subroutine write_vtk_header(xml)
    use mod_streams
    use mod_sys
    use mod_probe
    implicit none

    character(len=16) :: curr_cycle
    character(len=1000) :: xml ! this allocation errors with nvhpc compiler
    character(len=5) :: gloabl_x_end, global_x_start
    integer:: i, j, k

    write(global_x_start, "(I5)")  (nx *nrank) +1
    write(gloabl_x_end, "(I5)")  nx * (nrank+1)

    xml = '<?xml version="1.0"?>'//new_line('a')//' &
       & <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//new_line('a')//' &
       ! start whole extent
       &  <RectilinearGrid WholeExtent="'//trim(adjustl(global_x_start))//' &
       & '//trim(adjustl(gloabl_x_end))//' 1 '//trim(adjustl(nystr))//' 1 1">'//new_line('a')//' &
       ! finished whole extent
       &   <Piece Extent="'//trim(adjustl(global_x_start))//' '//trim(adjustl(gloabl_x_end))//'& 
       & 1 '//trim(adjustl(nystr))//' 1 1">'//new_line('a')//' &
       &    <Coordinates>'//new_line('a')//' &
       &     <DataArray type="Float64" NumberOfComponents="1" Name="X" format="ascii"> '

    write(23, "(a)", advance="no") trim(adjustl(xml))

    ! this indexing ensures that we are writing only the information that we are directly
    ! in control of
    ! for a nxmax = 1000 and 4 mpi process splitting this in the x direction nx = 250
    ! nrank = 0 would cover ranges i = 1, 250
    ! nrank = 1 would cover ranges i = 251, 500
    ! etc
    do i = (nx * nrank) +1 ,nx * (nrank + 1)
        write(curr_cycle, "(E15.10)") xg(i)
        write(23, "(A1, a)", advance="no") ' ', curr_cycle
    enddo

    xml = '</DataArray>' // new_line('a') 
    xml = trim(xml) // '      <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="ascii">'
    write(23, "(a)", advance="no") trim(xml)

    do j = 1,nymax
        write(curr_cycle, "(E15.10)") yg(j)
        write(23, "(a, A15)", advance="no") ' ', curr_cycle
    enddo


    xml = '</DataArray>' // new_line('a')  

    xml = trim(xml) // '      <DataArray type="Float64" NumberOfComponents="1" name="Z" format="ascii">0.0</DataArray>'  
    xml = trim(xml) // new_line('a')
    xml = trim(xml) // '     </Coordinates>' // new_line('a')
    xml = trim(xml) // '   <PointData>' // new_line('a')
    write(23, '(a)') xml

end subroutine write_vtk_header

subroutine write_vtk_dataarray(xml, slice_var)
    use mod_streams
    use mod_probe
    implicit none

    integer :: slice_var, i, j, k
    character(len=1000) :: xml
    character(len=6) :: variable_name
    character(len=16) :: curr_cycle

    if (slice_var == 1) then
        variable_name = "rho"
        call move_alloc(span_average_rho, span_average)
    else if (slice_var == 2) then 
        variable_name = "u"
        call move_alloc(span_average_u, span_average)
    else if (slice_var == 3) then 
        variable_name = "v"
        call move_alloc(span_average_v, span_average)
    else if (slice_var == 4) then 
        variable_name = "w"
        call move_alloc(span_average_w, span_average)
    else if (slice_var == 5) then 
        variable_name = "energy"
        call move_alloc(span_average_w, span_average)
    endif

    write(23, '(a)',advance="no") '    <DataArray type="Float64" NumberOfComponents="1" Name="'&
        & //trim(variable_name)//'" format="ascii">'


    do j = 1, ny
        do i = 1, nx
            write(curr_cycle, '(E16.10)') span_average(i,j)
            !xml = trim(xml) // ' ' // curr_cycle // ' '
            write(23, "(a, A16)", advance="no") ' ', curr_cycle
        enddo
    enddo

    xml = '</DataArray>' // new_line('a') 
    write(23, '(a)') xml

end subroutine write_vtk_dataarray

subroutine write_vtk_footer(xml)
    character(len=1000) :: xml

    xml =  '    </PointData>' // new_line('a') // '   </Piece>'//new_line('a') 
    xml = trim(xml) // '  </RectilinearGrid>' //new_line('a') // '</VTKFile>'

    write(23, '(a)') trim(adjustl(xml))
end subroutine write_vtk_footer
