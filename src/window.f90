!*************************************************
!   WINDOW
!   Calculate the Spectral window-function
!   of a timeseries.
!
!   Input:
!       window [options] [inputfile] [outputfile]
!
!   Options:
!       -tday      Treat times as days (default)
!       -tsec      Treat times as seconds
!       -kplrraw   Input file is a Kepler data file (use raw columns)
!       -kplrpdc   Input file is a Kepler data file (use PDC columns)
!       -kplrwg    Input file is a Kepler data file (use WG columns)
!       -quiet     Print nothing to screen
!       -version   Print program version information
!
!   Written by
!     Rasmus Handberg
!     Department of Physics and Astronomy
!     Aarhus University
!*************************************************

program window
	use kind_defs
	use numinteg
	use timeseries
	implicit none
	real(dbl) :: nuw, dnu, dF
	real(dbl) :: nu_min, nu_max, nu
	integer :: N, j, Nspec, err
	logical :: bolUseWeights
	character(200) :: buffer, inputfile = "", outputfile = ""
	integer :: TimeType = 1 ! Default time type is days
	integer :: inputtype = 0
	real(dbl) :: t1, t2, Psum, Pint, omegaw
	real(dbl), allocatable :: win(:), fs(:), fc(:)
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
			case ("-kplrraw")
				TimeType = 1
				inputtype = 1
			case ("-kplrpdc")
				TimeType = 1
				inputtype = 2
			case ("-kplrwg")
				TimeType = 1
				inputtype = 3
			case ("-quiet")
				quiet = .true.
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
	
	if (.not. quiet) then
		print *, "--------------------------------"
		print *, "            WINDOW              "
		print *, "   Calculate window-function    "
		print *, "--------------------------------"
	endif
	
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
			outputfile = trim(outputfile(1:j-1)) // '.window'
		else
			outputfile = trim(outputfile) // '.window'
		endif
	endif

	! Import data:
	if (.not. quiet) then
		print *, "Reading data..."
    endif
    open(50, file=inputfile, status='old', action='read')
        ! New count:
        call CountData(50, N, bolUseWeights, inputtype)
        ! Allocate the data:
    	allocate(t(1:N), f(1:N), w(1:N), fs(1:N), fc(1:N), stat=err)
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
        call ImportData(50, bolUseWeights, t, f, w, inputtype)
    close(50)
    
    ! Convert time vector to seconds, based on the selected input:
    if (TimeType == 1) then
    	t = t*86400_dbl
    endif
    
    ! Print information about series:
	call print_info(t)
    
	! Frequencies to calculate the window for:
	print *, ""
	print *, "Windowfunction Frequency (microHz):"
	read *, nuw
	print *, "Frequency interval (microHz):"
	read *, dnu
	print *, "Frequency resolution (microHz):"
	read *, dF
	print *, ""

	! A bit of frequency-calculations:
	call cpu_time(t1)
	nu_min = nuw - dnu
	Nspec = ceiling(dnu/dF)
	if (Nspec <= 0 .or. dnu <= 0 .or. nu_min < 0) stop "Wrong frequencies!"

	! Allocate space for windowfunction:
	allocate(win(1:(2*Nspec+1)), stat=err)
	if (err /= 0) stop "Couldn't allocate memory!"

	! Calculate fake datasets:
	omegaw = pi*nuw*2e-6_dbl
	fs(:) = w(:)*sin(omegaw*t(:)) ! Sine-function "dataset"
	fc(:) = w(:)*cos(omegaw*t(:)) ! Cosine-function "dataset"

	! Write the windowfunction to the outputfile:
	if (.not. quiet) then
		print *, "--------------------------------"
	endif
	!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nu,j) SCHEDULE(STATIC)
	do j = -Nspec,Nspec
		nu = nuw+j*dF

		win(Nspec+1+j) = windowfunction(t, w, nuw, nu, fs, fc)

		call Progress2(2*Nspec+1)
	enddo
	!$OMP END PARALLEL DO

	! Calculate sum and integral of window:
	if (.not. quiet) then
		Psum = sum(win)
		Pint = simpson(dF, win)
		print *, ''
		print *, "--------------------------------"
		print *, 'Psum = ', Psum
		print *, 'Pint = ', Pint
    endif

	! Write spectrum to file:
    open(51, file=outputfile, status='replace', action='write')
	write(51, '(a)') '#-----------------------------------------'
	write(51, '(a)') '# WINDOW'
	write(51, '(a)') '# ' // trim(inputfile)
	write(51, '(a,f8.2,a)') '# Window frequency: ', nuw, ' microHz'
	write(51, '(a,f6.4,a)') '# Resolution: ', dF, ' microHz'
	write(51, '(a)') '# Frequencies in microHz, Window (power)'
	write(51, '(a)') '#-----------------------------------------'
    do j = -Nspec,Nspec
    	nu = j*dF
        write(51, '(f17.10,ES26.16E3)') nu, win(Nspec+1+j)
    enddo
    close(51)

	! Clean the memory:
	deallocate(t,f,w,win,fs,fc)
	
	call cpu_time(t2)
    call print_time(t2-t1)
end program window
