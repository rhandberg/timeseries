!****************************************************
!*                 SCLEAN
!*  Simultaneous Iterative Sine Wave Fitting
!*
!*  Description:
!*    Cleans a timeseries for coherent signals,
!*    taking into account the influences of the
!*    frequencies on each other.
!*
!*  Input:
!*    sclean [options] [inputfile] [outputfile]
!*
!*  Options:
!*    -tday      Treat times as days (default).
!*    -tsec      Treat times as seconds.
!*    -version   Print program version information.
!*
!*  Written by
!*    Rasmus Handberg
!*    Department of Physics and Astronomy
!*    Aarhus University
!****************************************************

program sclean
	use kind_defs
	use timeseries
	use FindMax
	implicit none
	character(200) :: buffer, inputfile = '', outputfile = ''
	logical :: bolUseWeights
	integer :: i, j, err, N, Nclean, k, Nspectra, total_spectra
	real(dbl), allocatable :: alpha(:), beta(:), nu(:)
	real(dbl) :: nu_min, nu_max, fmean
	integer :: TimeType = 1 ! Default time type is days
   	real(dbl) :: t1, t2
   	real(dbl) :: nu_low, nu_high, dF, Pmax, phase, nyquist
	integer :: nargs, iargc

	print *, "--------------------------------------------"
	print *, "                  SCLEAN                    "
	print *, "  Simultaneous Iterative Sine Wave Fitting  "
	print *, "--------------------------------------------"

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
	
        ! Import the data:
        call ImportData(50, bolUseWeights, t, f, w, 0)
    close(50)
    
	! Convert time vector to seconds, based on the selected input:
	if (TimeType == 1) then
		t = t*86400_dbl
	endif
	
    ! Print information about series:
	! Get fundamental frequency resolution in case of equividistant points.
	! A good measure of the resolution in the spectrum.
    call print_info(t, nyquist, dF)

    ! Select frequencies:
    print *, ""
    print *, 'min-frequency (microHz):'
    read *, nu_min
    print *, 'max-frequency (microHz):'
    read *, nu_max
    
    ! Number of points to clean:
    print *, ""
    print *, "Number of frequencies to clean:"
    read *, Nclean
    print *, ""

	! Subtract mean:
	call cpu_time(t1)
	fmean = mean(f,w)
	f = f - fmean

	! Allocate space for the results:
   	allocate(alpha(1:Nclean), beta(1:Nclean), nu(1:Nclean), stat=err)
    if (err /= 0) stop "Couldn''t allocate memory!"

	! Calculate the total number of spectra needed to be calculated:
	total_spectra = Nclean-1 + Nclean*(Nclean+1)/2.0

	print '("  A total of ",i3," spectra needed.")', total_spectra
	print *, "--------------------------------------------"

	Nspectra = 0
	do i = 1,Nclean
		! Find the maximum in the spectrum and subtract it:
		call FindSpectrumMax(t, f, w, nu_min, nu_max, alpha(i), beta(i), nu(i))
		f = f - alpha(i)*sin(twopi_micro*nu(i)*t) - beta(i)*cos(twopi_micro*nu(i)*t)
		Nspectra = Nspectra+1
		print '("  Spectrum ",i3,"/",i3," done.")', Nspectra, total_spectra

		if (i /= 1) then
			do j = 1,i
				! Add the oscillation:
				f = f + alpha(j)*sin(twopi_micro*nu(j)*t) + beta(j)*cos(twopi_micro*nu(j)*t)

				! Narrow search (nu_min,nu_max) to area arround nu(j)
				nu_low  = max(nu_min, nu(j)-5*dF)
				nu_high = min(nu_max, nu(j)+5*dF)
				
				! Find the maximum:
				call FindSpectrumMax(t, f, w, nu_low, nu_high, alpha(j), beta(j), nu(j))
				Nspectra = Nspectra+1
				print '("  Spectrum ",i3,"/",i3," done.")', Nspectra, total_spectra

				! Remove the oscillation again:
				f = f - alpha(j)*sin(twopi_micro*nu(j)*t) - beta(j)*cos(twopi_micro*nu(j)*t)
			enddo
		endif
	enddo
	print *, '--------------------------------------------'

	! Add the mean again and convert the time back to original format:
	f = f + fmean
	if (TimeType == 1) then
		t = t/86400_dbl
	endif

	! Write the found frequencies and amplitudes to file:
	open(51, file='clean.log', action='write', position='append')
	! Write header to clean.dat file:
	write(51, '(a)') '#----------------------------------'
	write(51, '(a)') '# SCLEAN'
	write(51, '(a)') '# ' // trim(inputfile)
	write(51, "('# min: ',f10.2,' microHz')") nu_min
	write(51, "('# max: ',f10.2,' microHz')") nu_max
	write(51, "('# Nclean: ',i7)") Nclean
	write(51, '(a)') '# Number, Power, Frequencies (microHz), Phase'
	write(51, '(a)') '#----------------------------------'
	do i = 1,Nclean
		Pmax = alpha(i)**2 + beta(i)**2
		phase = atan2(alpha(i),beta(i))
		write(51,'(i5,3ES26.16E3)') i, Pmax, nu(i), phase
		print '(i5,3f14.4)', i, Pmax, nu(i), phase
	enddo
	close(51)

	! Save in file:
	open(52, file=outputfile, status='replace', action='write')
	do i = 1,N
		write(52, '(3ES26.16E3)') t(i), f(i), 1.0_dbl/sqrt(w(i))
	enddo
	close(52)

	! Clean up the memory:
	deallocate(t, f, w, alpha, beta, nu)

	call cpu_time(t2)
	print *, '--------------------------------------------'
    call print_time(t2-t1)
end program
