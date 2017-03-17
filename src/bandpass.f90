!****************************************************
!   BANDPASS
!   Bandpass filtering using Weighted Least Squares
!
!   Input:
!       bandpass [options] [inputfile] [outputfile]
!
!   Options:
!       -tday      Treat times as days (default).
!       -tsec      Treat times as seconds.
!       -version   Print program version information.
!
!   Written by
!     Rasmus Handberg
!     Department of Physics and Astronomy
!     Aarhus University
!****************************************************

program bandpas
	use kind_defs
	use timeseries
	implicit none
	real(dbl), allocatable :: bp(:)
	real(dbl) :: dF, nu_max, nu1, nu2, nu_min
	real(dbl) :: t1, t2
	integer :: N, err, j
	logical :: bolUseWeights
	character(200) :: buffer, inputfile, outputfile
   	integer :: TimeType = 1 ! Default time type is days
	integer :: nargs, iargc	
	
	! Extract parameters from command-line:
	nargs = iargc()
	do j = 1,nargs
		call getarg(j, buffer)
		select case (trim(buffer))
			case ("-tsec")
				TimeType = 0
			case ("-tday")
				TimeType = 1
			case ("-version")
				print *, "Current version:"
				print *, "$Revision: 119 $"
				print *, "$Date: 2016-02-09 09:57:56 +0100 (ti, 09 feb 2016) $"
				stop

			case default
				buffer = trim(buffer)
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
	
	print *, "--------------------------------"
    print *, "         BANDPASS-filter        "
    print *, "--------------------------------"
	
	! If no files were given, ask for input and make output:
	if (inputfile == "") then
        print *, 'Enter input-file:'
        read '(a50)', inputfile
	endif
	if (outputfile == '') then
		outputfile = 'bandpass.dat'
	endif

    print *, "Reading data..."
    open(50, file=inputfile, status='old', action='read')
        ! Count number of datapoints:
        call CountData(50, N, bolUseWeights)
        ! Allocate the data:d
    	allocate(t(1:N), f(1:N), w(1:N), bp(1:N), stat=err)
    	if (err /= 0) stop "Couldn't allocate memory!"
	    ! Ask to use the statistical weights:
		if (bolUseWeights) then
			print *, "Do you want to use the statistical weights? (1=yes,0=no)"
			read(*, '(i1)') j
			if (j == 0) bolUseWeights = .false.
		endif
        ! Import the data:
        call ImportData(50, bolUseWeights, t, f, w)
    close(50)

	! Convert time vector to seconds, based on the selected input:
	if (TimeType == 1) then
		t = t*86400_dbl
	endif

    ! Print information about series:
	call print_info(t)

	! Ask for frequencies:
	print *, ''
	print *, "min-frequency: (microHz)"
	read *, nu_min
	print *, "max-frequency: (microHz)"
	read *, nu_max
	print *, "Frequency resolution: (microHz)"
	read *, dF
	print *, "lower frequency: (microHz)"
	read *, nu1
	print *, "higher frequency: (microHz)"
	read *, nu2
	print *, ''

	! Run the bandpass-filter:
	call cpu_time(t1)
	bp = bandpass(t, f, w, nu1, nu2, nu_min, nu_max, dF)
	
	! Convert time vector back to the original units, based on the selected input:
	if (TimeType == 1) then
		t = t/86400_dbl
	endif
	
	! Write results to file:
	open(51, file=outputfile, status='replace', action='write')
	do j = 1,N
		write(51, '(3ES26.16E3)') t(j), bp(j), 1.0_dbl/sqrt(w(j))
	enddo
	close(51)
	
	! Clean up memory:
	deallocate(t, f, w, bp)
	
	call cpu_time(t2)
	print *, "--------------------------------"
    call print_time(t2-t1)	
end program bandpas
