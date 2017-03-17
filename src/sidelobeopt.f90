!****************************************************
!*             Sidelobe Optimization
!*	      
!*  Rasmus Handberg
!*  Department of Physics and Astronomy
!*  Aarhus University
!*
!*  Input:
!*       sidelobeopt [options] [inputfile] [outputfile]
!*
!*  Options:
!*       -tsec  Treat times as seconds (default)
!*       -tday  Treat times as days
!*
!****************************************************

program sidelobeopt
    use kind_defs
	use timeseries
    implicit none
    integer :: N, i, j, k, err, Nclean
    real(dbl) :: nu_min, nu_max
    real(dbl) :: nu, numax, fmean, Amax, alpha, beta, phase
    character(200) :: buffer, inputfile = '', outputfile = ''
    logical :: bolUseWeights
	integer :: TimeType = 0 ! Default time type is seconds
   	real(dbl) :: t1, t2, ti
	integer :: nargs, iargc
	!real(dbl) :: gain = 1.0_dbl

    print *, "--------------------------------"
    print *, "             CLEAN              "
    print *, "  Iterative Sine Wave Fitting   "
    print *, "--------------------------------"

	! Extract parameters from command-line:
	nargs = iargc()
	do j = 1,nargs
		call getarg(j, buffer)
		select case (trim(buffer))
			case ("-tsec")
				TimeType = 0
			case ("-tday")
				TimeType = 1

			case default
				if (.not. buffer(1:1) == "-" .and. inputfile == "") then
					inputfile = trim(buffer)
				elseif (.not. buffer(1:1) == "-" .and. outputfile == "") then
					outputfile = trim(buffer)
				else
					print *, "UNKNOWN ARGUMENT: " // trim(buffer)
					stop
				endif
		end select
	enddo
	! If no files were given, ask for input and make output:
	if (inputfile == "") then
        print *, 'Enter input-file:'
        read '(a50)', inputfile
	endif
	if (outputfile == '') then
		j = scan(inputfile, '/', back=.true.)
		if (j /= 0) then
			outputfile = trim(inputfile(j+1:)) // '.clean'
		else
			outputfile = trim(inputfile) // '.clean'
		endif
	endif

    print *, "Reading data..."
    open(50, file=inputfile, status='old', action='read')		
        ! New count:
        call CountData(50, N, bolUseWeights)
	
        ! Allocate the data:
    	allocate(t(1:N), f(1:N), w(1:N), stat=err)
    	if (err /= 0) stop "Couldn''t allocate memory!"
        ! Import the data:
        call ImportData(50, bolUseWeights, t, f, w)
    close(50)
    
	! Convert time vector to seconds, based on the selected input:
	if (TimeType == 1) then
		t = t*86400_dbl
	endif

    ! Print information about series:
	call print_info(t)
    

	! Subtract mean:
	call cpu_time(t1)


	smooth(w, ???)

	
	! Save in file:
	open(52, file=outputfile, status='replace', action='write')
	do i = 1,N
		! TODO: Stuid to do this in the loop:
		if (TimeType == 1) then
			ti = t(i)/86400_dbl
		else
			ti = t(i)
		endif
	
		write(52, '(3ES22.12E3)') ti, f(i), 1.0_dbl/sqrt(w(i))
	enddo
	close(52)
	
	! Clean up memory:
	deallocate(t, f, w)
	
	call cpu_time(t2)
	print *, "---------------------------------"
    call print_time(t2-t1)
end program sidelobeopt
