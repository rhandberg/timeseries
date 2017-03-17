!****************************************************
!*             CLEAN
!*  Iterative Sine Wawe Fitting
!*
!*  Input:
!*       clean [options] [inputfile] [outputfile]
!*
!*  Options:
!*       -tday      Treat times as days (default)
!*       -tsec      Treat times as seconds
!*       -version   Print program version information
!*
!*  Written by
!*    Rasmus Handberg
!*    Department of Physics and Astronomy
!*    Aarhus University
!****************************************************

program clean
    use kind_defs
	use timeseries
	use FindMax
    implicit none
    integer :: N, i, j, k, err, Nclean
    real(dbl) :: nu_min, nu_max
    real(dbl) :: nu, numax, fmean, Amax, alpha, beta, phase
    character(200) :: buffer, inputfile = '', outputfile = ''
    logical :: bolUseWeights
	integer :: TimeType = 1 ! Default time type is days
   	real(dbl) :: t1, t2, ti
	integer :: nargs, iargc

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
			case ("-version")
				print *, "Current version:"
				print *, "$Revision: 119 $"
				print *, "$Date: 2016-02-09 09:57:56 +0100 (ti, 09 feb 2016) $"
				stop

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
			outputfile = trim(inputfile(j+1:))
		else
			outputfile = trim(inputfile)
		endif
		j = scan(outputfile, '.', back=.true.)
		if (j > 1) then
			outputfile = trim(outputfile(1:j-1)) // '.clean'
		else
			outputfile = trim(outputfile) // '.clean'
		endif
	endif

    print *, "Reading data..."
    open(50, file=inputfile, status='old', action='read')		
        ! New count:
        call CountData(50, N, bolUseWeights, 0)
	
        ! Allocate the data:
    	allocate(t(1:N), f(1:N), w(1:N), stat=err)
    	if (err /= 0) stop "Couldn''t allocate memory!"
    	
    	! Ask to use the statistical weights:
		if (bolUseWeights) then
			print *, "Do you want to use the statistical weights? (1=yes,0=no)"
			read(*, '(i1)') j
			if (j == 0) bolUseWeights = .false.
		endif
		if (.not. quiet) then
			if (bolUseWeights) then
				print *, "Using weights"
			else
				print *, "Not using weights"
			endif
		endif	
    	
        ! Import the data:
        call ImportData(50, bolUseWeights, t, f, w, 0)
    close(50)
    
	! Convert time vector to seconds, based on the selected input:
	if (TimeType == 1) then
		t = t*86400_dbl
	endif

    ! Print information about series:
	call print_info(t)
    
    ! Select frequencies:
    print *, ''
    print *, 'min-frequency (microHz):'
    read *, nu_min
    print *, 'max-frequency (microHz):'
    read *, nu_max
    
    ! Number of points to clean:
    print *, ""
    print *, "Number of frequencies to clean:"
    read *, Nclean
	!print *, "Loop gain: (0-1)"
	!read *, gain
    print *, ""
    print *, "---------------------------------"
    !if (gain > 1.0 .or. gain <= 0.0) stop "Gain loop must be between 0 and 1."

	! Subtract mean:
	call cpu_time(t1)
	fmean = mean(f,w)
	f = f - fmean

    ! Open file to store frequencies:
    open(51, file='clean.log', action='write', position='append')
	! Write header to clean.dat file:
	write(51, '(a)') '#----------------------------------'
	write(51, '(a)') '# CLEAN'
	write(51, '(a)') '# ' // trim(inputfile)
	write(51, "('# min: ',f10.2,' microHz')") nu_min
	write(51, "('# max: ',f10.2,' microHz')") nu_max
	write(51, "('# Nclean: ',i7)") Nclean
	write(51, '(a)') '# Number, Power, Frequencies (microHz), Phase'
	write(51, '(a)') '#----------------------------------'
    do k = 1,Nclean
    	! Find maximum peak in spectrum:
    	call FindSpectrumMax(t, f, w, nu_min, nu_max, alpha, beta, numax)
		Amax = alpha**2 + beta**2
        
		! Write the progress:
		phase = atan2(alpha,beta)
		print '(i5,3f14.4)', k, Amax, numax, phase
		write(51,'(i5,3f14.4)') k, Amax, numax, phase
        
		! Clean the series with Iterativ Sine Fitting:
		f = f - alpha*sin(twopi_micro*numax*t) - beta*cos(twopi_micro*numax*t)

		! Save in file:
		open(52, file=outputfile, status='replace', action='write')
		do i = 1,N
			! TODO: Stuid to do this in the loop:
			if (TimeType == 1) then
				ti = t(i)/86400_dbl
			else
				ti = t(i)
			endif
		
			write(52, '(3ES22.12E3)') ti, f(i)+fmean, 1.0_dbl/sqrt(w(i))
		enddo
		close(52)
    enddo
	close(51)
	
	! Clean up memory:
	deallocate(t, f, w)
	
	call cpu_time(t2)
	print *, "---------------------------------"
    call print_time(t2-t1)
end program clean
