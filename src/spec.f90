!****************************************************
!   Weighted Least Squares Spectrum
!
!   Input:
!       spec [options] [inputfile] [outputfile]
!
!   Options:
!       -p         Power spectrum (default)
!       -a         Amplitude spectrum
!       -pd        Power density spectrum
!       -tday      Treat times as days (default)
!       -tsec      Treat times as seconds
!       -kplrraw   Input file is a Kepler data file (use raw columns)
!       -kplrpdc   Input file is a Kepler data file (use PDC columns)
!       -kplrwg    Input file is a Kepler data file (use WG columns)
!       -auto      Automatic. No interaction required (experimental)
!       --ofac=number
!       -quiet     Print nothing to screen
!       -version   Print program version information
!
!   Written by
!     Rasmus Handberg
!     Department of Physics and Astronomy
!     Aarhus University
!****************************************************

program spec
	use kind_defs
	use timeseries
	implicit none
	integer :: N, j, Nspec, err
	real(dbl) :: nu_min, nu_max, dF, Pint
	real(dbl) :: nu, Tot, numax, Amax = 0
	real(dbl) :: alpha, beta, priv_max, priv_numax
	real(dbl), allocatable :: A(:), delta(:)
	real(dbl) :: ofac = 2 ! Oversampling factor when using auto
	character(200) :: buffer, inputfile = '', outputfile = '', StrSpectrumType
	logical :: bolUseWeights, UseFundamental = .false., RunAuto = .false.
	integer :: PPMNormalisation = 0 ! Default normalizaton to PPM is OFF
   	integer :: SpectrumType = 1 ! Default spectrum type is power
   	integer :: TimeType = 1 ! Default time type is days
	integer :: inputtype = 0
	real(dbl) :: t1, t2
	integer :: nargs, iargc

	! Extract parameters from command-line:
	nargs = iargc()
	do j = 1,nargs
		call getarg(j, buffer)
		select case (trim(buffer))
			case ("-a")
				SpectrumType = 0
			case ("-p")
				SpectrumType = 1
			case ("-pd")
				SpectrumType = 2
			case ("-tsec")
				TimeType = 0
			case ("-tday")
				TimeType = 1
			case ("-auto")
				RunAuto = .true.
			case ("-kplrraw")
				TimeType = 1
				PPMNormalisation = 1
				inputtype = 1
			case ("-kplrpdc")
				TimeType = 1
				PPMNormalisation = 1
				inputtype = 2
			case ("-kplrwg")
				TimeType = 1
				PPMNormalisation = 1
				inputtype = 3
			case ("-quiet")
				quiet = .true.
			case ("-version")
				print *, "Current version:"
				print *, "$Revision: 119 $"
				print *, "$Date: 2016-02-09 09:57:56 +0100 (ti, 09 feb 2016) $"
				stop

			case default
				buffer = trim(buffer)
				if (scan(buffer, "--ofac=") == 1) then
					buffer = buffer(8:)
					read(buffer, *) ofac
				else
					if (.not. buffer(1:1) == "-" .and. inputfile == "") then
						inputfile = trim(buffer)
					elseif (.not. buffer(1:1) == "-" .and. outputfile == "") then
						outputfile = trim(buffer)
					else
						print *, "UNKNOWN ARGUMENT: " // trim(buffer)
						stop
					endif
				endif
		end select
	enddo
	
	if (.not. quiet) then
		print *, "---------------------------------"
		print *, " Weighted Least Squares Spectrum "
		print *, "---------------------------------"
	endif
	
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
			outputfile = trim(outputfile(1:j-1)) // '.spec'
		else
			outputfile = trim(outputfile) // '.spec'
		endif
	endif
	
	! Read the data from the input-file:
	if (.not. quiet) then
		print *, "Reading data..."
	endif
	open(50, file=inputfile, status='old', action='read')
		! Count number of datapoints:
		call CountData(50, N, bolUseWeights, inputtype)
	
		! Allocate the data:
		allocate(t(1:N), f(1:N), w(1:N), stat=err)
		if (err /= 0) stop "Couldn't allocate memory!"
	
		! Ask to use the statistical weights:
		if (bolUseWeights .and. .not. RunAuto) then
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
		call ImportData(50, bolUseWeights, t, f, w, inputtype)
	close(50)
    
	! Convert time vector to seconds, based on the selected input:
	if (TimeType == 1) then
		t = t*86400_dbl
	endif
    
	! Apply desired scaling to PPM:
    if (PPMNormalisation == 1) then
    	f = 1e6_dbl*(f/median(f) - 1.0_dbl)
    elseif (PPMNormalisation == 2) then
    	f = 1e6_dbl*(f - 1.0_dbl)
    elseif (PPMNormalisation == 3) then
    	f = 1e6_dbl*f
    endif
    
	! Print information about series:
	call print_info(t, nu_max, dF)
    
	! Select frequencies:
	if (RunAuto) then
		nu_min = 0.0_dbl
		Pint = dF
		dF = dF/ofac
		if (.not. quiet) then
			print *, 'WARNING: Fundamental resolution not fully supported!'
			print *, 'Oversampling by factor of ' , ofac , ': ', dF
		endif
	else
		print *, ''
		print *, 'min-frequency (microHz):'
		read *, nu_min
		print *, 'max-frequency (microHz):'
		read *, nu_max
		print *, 'DeltaF (microHz):'
		read *, dF
		if (SpectrumType == 2) then
			print *, 'Integral of Spectral Window function (1/microHz):'
			read *, Pint
		endif
	endif

	! Calculate number of points in spectrum:
	call cpu_time(t1)
	if (nu_min == 0) nu_min = nu_min + dF ! Don't calculate nu=0
	Nspec = ceiling((nu_max-nu_min)/dF)+1
	if (Nspec <= 0 .or. nu_min < 0 .or. nu_max < 0) stop "Wrong frequencies!"
	
	if (.not. quiet) then
		print *, ""
		print *, "Points in spectrum:", Nspec
		print *, "--------------------------------"
	endif

	! Allocate memory for results:
	allocate(A(1:Nspec), delta(1:Nspec), stat=err)
	if (err /= 0) stop "Couldn't allocate memory!"

	! Subtract the mean from the data:
	f = f - mean(f,w) ! Weighted mean

	! Calcultate spectrum:
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(alpha,beta,nu,priv_max,priv_numax)
	priv_max = 0.0_dbl
	!$OMP DO SCHEDULE(STATIC)
	do j = 1,Nspec
		nu = nu_min+dF*(j-1)
      
		! Calculate alpha and beta:
		call CalcAlphaBeta(t,f,w,nu,alpha,beta)
		
		! Calculate power and phase:
		A(j) = alpha**2 + beta**2
		delta(j) = atan2(alpha,beta)
        
		! Save the local maximum in the spectrum:
		if (A(j) > priv_max) then
			priv_max = A(j)
			priv_numax = nu
		endif
		
		! Update progressbar:
		call Progress2(Nspec)
	enddo
	!$OMP END DO NOWAIT
	!$OMP FLUSH (Amax)
	if (priv_max > Amax) then
		!$OMP CRITICAL
			if (priv_max > Amax) then
				Amax = priv_max
				numax = priv_numax
			endif
		!$OMP END CRITICAL
	endif
	!$OMP END PARALLEL
    
	! Convert to the desired spectrum type
	! and construct spectrum type string:
	if (SpectrumType == 1) then
		StrSpectrumType = 'Power spectrum'
	elseif (SpectrumType == 2) then
		StrSpectrumType = 'Power density spectrum'
		A = A/Pint
		Amax = Amax/Pint
	else
		StrSpectrumType = 'Amplitude spectrum'
		A = sqrt(A)
		Amax = sqrt(Amax)
	endif
	if (bolUseWeights) then
		StrSpectrumType = 'Weighted ' // StrSpectrumType
	endif

	if (.not. quiet) then
		print *, ''
		print *, "--------------------------------"
		print '(" Amax  = ",f10.3)', Amax
		print '(" Numax = ",f10.3," microHz")', numax
    endif
	
	! Write spectrum to file:
	open(51, file=outputfile, status='replace', action='write')
		write(51, '(a)') '#-----------------------------------------'
		write(51, '(a)') '# SPEC'
		write(51, '(a)') '# ' // trim(inputfile)
		write(51, '(a)') '# Type: ' // trim(StrSpectrumType)
		if (SpectrumType == 2) then
			write(51, '(a,f6.4,a)') '# Pint: ', Pint, ' (1/microHz)'
		endif
		write(51, '(a,f6.4,a)') '# Resolution: ', dF, ' microHz'
		write(51, '(a)') '# Frequencies in microHz, Spectrum, Phases'
		write(51, '(a)') '#-----------------------------------------'
		do j = 1,Nspec
			nu = nu_min+dF*(j-1)
			write(51, '(3ES26.16E3)') nu, A(j), delta(j)
		enddo
	close(51)

	! Clean up memory:
	deallocate(A, delta, t, f, w)
    
    if (.not. quiet) then
		call cpu_time(t2)
		print *, ''
		call print_time(t2-t1)
	endif
end program spec
