!****************************************************
!*             LCLEAN
!*  Iterative Sine Wawe Fitting
!*
!*  Description:
!*    Cleans a timeseries until given limit.
!*
!*  Input:
!*    lclean [options] [inputfile] [outputfile]
!*
!*  Options:
!*    -p         Power spectrum (default)
!*    -a         Amplitude spectrum
!*    -pd        Power density spectrum
!*    -tday      Treat times as days (default)
!*    -tsec      Treat times as seconds
!*    -version   Print program version information
!*
!*  Written by
!*    Rasmus Handberg
!*    Department of Physics and Astronomy
!*    Aarhus University
!****************************************************

program lclean
    use kind_defs
	use timeseries
	use FindMax
    implicit none
	integer, parameter :: Nmax = 1000
    integer :: N, i, j, k, err
    real(dbl) :: nu_min, nu_max, Pint
    real(dbl) :: nu, numax, fmean, Amax, alpha, beta, limit, phase
    character(200) :: buffer, inputfile = "", outputfile = ""
    logical :: bolUseWeights
   	real(dbl) :: t1, t2, ti
   	integer :: SpectrumType = 1 ! Default spectrum type is power
	character(50) :: StrLimitInput = 'Power'
	integer :: TimeType = 1 ! Default time type is days
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
			case ("-a")
				SpectrumType = 0
				StrLimitInput = "Amplitude"
			case ("-p")
				SpectrumType = 1
				StrLimitInput = "Power"
			case ("-pd")
				SpectrumType = 2
				StrLimitInput = "Power density"
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
	if (inputfile == '') then
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
    print *, ""
    print *, 'min-frequency (microHz):'
    read *, nu_min
    print *, 'max-frequency (microHz):'
    read *, nu_max
	if (SpectrumType == 2) then
		print *, 'Integral of Spectral Window function (1/microHz):'
		read *, Pint
	endif
    
    ! Number of points to clean:
    print *, ""
    print *, "Limit to clean frequencies above: (" // trim(StrLimitInput) // ")"
    read *, limit
    print *, ""
    print *, "---------------------------------"
	
	! Subtract mean:
	call cpu_time(t1)
	fmean = mean(f,w)
	f = f - fmean

	! Open file to store frequencies:
	k = 1
	open(51, file='clean.log', action='write', position='append')
	! Write header to clean.dat file:
	write(51, '(a)') '#----------------------------------'
	write(51, '(a)') '# LCLEAN'
	write(51, '(a)') '# ' // trim(inputfile)
	write(51, "('# min: ',f10.2,' microHz')") nu_min
	write(51, "('# max: ',f10.2,' microHz')") nu_max
	write(51, "('# limit: ',f10.2,' (',a,')')") limit, trim(StrLimitInput)
	write(51, '(a)') '# Number, ' // trim(StrLimitInput) // ', Frequencies (microHz), Phase'
	write(51, '(a)') '#----------------------------------'
	do
		! Find maximum peak in spectrum:
		call FindSpectrumMax(t, f, w, nu_min, nu_max, alpha, beta, numax)
		Amax = alpha**2 + beta**2
		
		! Convert the limit to power, depending on the input:
		if (SpectrumType == 0) then
			Amax = sqrt(Amax)
		elseif (SpectrumType == 2) then
			Amax = Amax/Pint
		endif
		
		write(*, '(a1)', advance='NO') char(13)
		if (Amax < limit .or. k > Nmax) exit

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
		
			write(52, '(3ES26.16E3)') ti, f(i)+fmean, 1.0_dbl/sqrt(w(i))
		enddo
		close(52)
		
		k = k + 1
	enddo
	close(51)
	
	! Clean up memory:
	deallocate(t, f, w)
	
	call cpu_time(t2)
	print *, "---------------------------------"
	call print_time(t2-t1)
end program lclean
