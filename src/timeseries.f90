!****************************************
!* Analysis of Timeseries-data
!* 
!* Rasmus Handberg
!* Department of Physics and Astronomy
!* University of Aarhus
!****************************************

module kind_defs
	implicit none
	integer, parameter :: &
		short = SELECTED_INT_KIND(5),  &
		long  = SELECTED_INT_KIND(10), &
		sgl   = SELECTED_REAL_KIND(6), &
		dbl   = SELECTED_REAL_KIND(14), &
		quad  = SELECTED_REAL_KIND(30)
	      
	real(dbl), parameter :: pi = 3.14159265358979323846_dbl
	real(dbl), parameter :: twopi = 6.2831853071795864769252867665590057683943387987502116419_dbl
	real(dbl), parameter :: twopi_micro = 6.2831853071795864769252867665590057683943387987502116419e-6_dbl
end module kind_defs

module timeseries
	use kind_defs
	implicit none
	real(dbl), allocatable :: t(:), f(:), w(:)
	real(dbl) :: sumWeights
	logical :: quiet = .false.
	
	interface
		subroutine CalcAlphaBeta(t, f, w, nu, alpha, beta)
			use kind_defs
			implicit none
			real(dbl), intent(in) :: t(:), f(:), w(:), nu
			real(dbl), intent(out) :: alpha, beta
		end subroutine

		subroutine CountItems(fid, N, bolUseWeights, inputtype)
			use kind_defs
			implicit none
			integer, intent(in) :: fid
			integer, intent(out) :: N
			logical, intent(out) :: bolUseWeights
			integer, intent(in), optional :: inputtype
		end subroutine

		subroutine ImportData(fid, bolUseWeights, t, f, w, inputtype)
			use kind_defs
			implicit none
			integer, intent(in) :: fid
			logical, intent(in) :: bolUseWeights
			real(dbl), dimension(:), intent(out) :: t, f, w
			integer, intent(in), optional :: inputtype
		end subroutine
		
		subroutine print_info(t, nyquist, dnu)
			use kind_defs
			real(dbl), intent(in) :: t(:)
			real(dbl), intent(out), optional :: nyquist
			real(dbl), intent(out), optional :: dnu
		end subroutine
		
		subroutine print_time(dt)
			use kind_defs
			real(dbl), intent(in) :: dt
		end subroutine
		
		function mean(f, w)
			use kind_defs
			implicit none
			real(dbl), intent(in) :: f(:)
			real(dbl), intent(in), optional :: w(:)
			real(dbl) :: mean
		end function mean
		
		function median(f)
			use kind_defs
			implicit none
			real(dbl), intent(in) :: f(:)
			real(dbl) :: median
		end function median
		
		function windowfunction(t, w, nuw, nu, fs, fc)
			use kind_defs
			real(dbl), intent(in) :: t(:), w(:), nuw, nu, fs(:), fc(:)
			real(dbl) :: windowfunction
		end function windowfunction
		
		function bandpass(t,f,w,nu1,nu2,nu_min,nu_max,dF)
			use kind_defs
			real(dbl), intent(in) :: t(:), f(:), w(:), nu1, nu2, dF, nu_max, nu_min
			real(dbl) :: bandpass(size(t))
		end function bandpass
		
		function lowpass(t,f,w,nu_lowpass,nu_min,nu_max,dF)
			use kind_defs
			real(dbl), intent(in) :: t(:), f(:), w(:), nu_lowpass, dF, nu_max, nu_min
			real(dbl) :: lowpass(size(t))
		end function lowpass
		
		function highpass(t,f,w,nu_highpass,nu_min,nu_max,dF)
			use kind_defs
			real(dbl), intent(in) :: t(:), f(:), w(:), nu_highpass, dF, nu_max, nu_min
			real(dbl) :: highpass(size(t))
		end function highpass
		
		function bandstop(t,f,w,nu1,nu2,nu_min,nu_max,dF)
			use kind_defs
			real(dbl), intent(in) :: t(:), f(:), w(:), nu1, nu2, dF, nu_max, nu_min
			real(dbl) :: bandstop(size(t))
		end function bandstop
	end interface
end module

! Count number of datapoints in data-file.
subroutine CountData(fid, N, bolUseWeights, inputtype)
	use kind_defs
	use strings, only : Count_Items, strlower
	implicit none
	integer, intent(in) :: fid
	integer, intent(out) :: N
	logical, intent(out) :: bolUseWeights
	integer, intent(in), optional :: inputtype
	logical :: AllWeightsOne
	real(dbl) :: wi
	integer :: ios, inptype = 0, Ncolumns, k
	character(200) :: line, tb, fb, wb, dummy1, dummy2, dummy3, dummy4

	if (present(inputtype)) then
		inptype = inputtype
	endif
	
	k = 1
	N = 0
	bolUseWeights = .true.
	AllWeightsOne = .true.
	! Count number of datapoints:
	do
		read(fid, '(a)', iostat=ios) line
		if (ios /= 0) exit
		if (len_trim(line) == 0) cycle ! Skip blank lines
		line = adjustl(line)
		if (line(1:1) == "#") cycle ! Skip lines starting with #
		
		! Count the number of columns in the file:
		Ncolumns = Count_Items(line)
		
		! Read time, data and weight from files, based on the specified format
		if (inputtype == 0 .and. Ncolumns == 2) then     ! Default files (3 columns)
			read(line, *) tb, fb
			wb = ""
		elseif (inputtype == 0 .and. Ncolumns >= 3) then     ! Default files (3 columns)
			read(line, *) tb, fb, wb
		elseif (inputtype == 1 .and. Ncolumns >= 3) then ! Kepler Raw data
			read(line, *) tb, fb, wb
		elseif (inputtype == 2 .and. Ncolumns >= 5) then ! Kepler PDC Corrected data
			read(line, *) tb, dummy1, dummy2, fb, wb
		elseif (inputtype == 3 .and. Ncolumns >= 7) then ! Kepler WG Corrected files
			read(line, *) tb, dummy1, dummy2, dummy3, dummy4, fb, wb
		else
			print *, "('*** ERROR: Unable to parse line ',i4,' in input file. ***')", k 
			print *, ""
			exit
		endif
		
		! Skip lines where the input is wrong:
		if (strlower(trim(tb)) == "-inf" .or. strlower(trim(tb)) == "inf" .or. strlower(trim(tb)) == "nan") cycle
		if (strlower(trim(fb)) == "-inf" .or. strlower(trim(fb)) == "inf" .or. strlower(trim(fb)) == "nan") cycle
		if (strlower(trim(wb)) == "-inf" .or. strlower(trim(wb)) == "inf" .or. strlower(trim(wb)) == "nan") cycle
		
		! Find out if the file contains weights:
		if (bolUseWeights) then
			if (trim(wb) == "") then
				bolUseWeights = .false.
			elseif (AllWeightsOne) then
				read(wb, *) wi
				if (wi /= 1.0_dbl) AllWeightsOne = .false.
			endif
		endif
		
		! Count the number of good datapoints:
		N = N + 1
		k = k + 1
	enddo
	! If all the weights were set to one, turn off the weights:
	if (AllWeightsOne) bolUseWeights = .false.
	rewind(fid)
end subroutine CountData

! Import data from file to arrays.
subroutine ImportData(fid, bolUseWeights, t, f, w, inputtype)
	use kind_defs
	use strings, only : Count_Items, strlower
	use timeseries, only: sumWeights
	implicit none
	integer, intent(in) :: fid
	logical, intent(in) :: bolUseWeights
	real(dbl), intent(out) :: t(:), f(:), w(:)
	integer, intent(in), optional :: inputtype
	integer :: ios, i = 1, inptype = 0, Ncolumns
	character(200) :: line, tb, fb, wb, dummy1, dummy2, dummy3, dummy4

	if (present(inputtype)) then
		inptype = inputtype
	endif

	do
		! Read the line from the file:
		read(fid, '(a)', iostat=ios) line
		if (ios /= 0) exit
		if (len_trim(line) == 0) cycle ! Skip blank lines
		line = adjustl(line)
		if (line(1:1) == "#") cycle ! Skip lines starting with #

		! Count the number of columns in the file:
		Ncolumns = Count_Items(line)

		! Read time, data and weight from files, based on the specified format
		if (inputtype == 0 .and. Ncolumns == 2) then
			read(line, *) tb, fb
			wb = ""
		elseif (inputtype == 0 .and. Ncolumns >= 3) then     ! Default files (3 columns)
			read(line, *) tb, fb, wb
		elseif (inputtype == 1 .and. Ncolumns >= 3) then ! Kepler Raw data
			read(line, *) tb, fb, wb
		elseif (inputtype == 2 .and. Ncolumns >= 5) then ! Kepler PDC Corrected data
			read(line, *) tb, dummy1, dummy2, fb, wb
		elseif (inputtype == 3 .and. Ncolumns >= 7) then ! Kepler WG Corrected files
			read(line, *) tb, dummy1, dummy2, dummy3, dummy4, fb, wb
		else
			exit
		endif
		
		! Skip lines where the input is wrong:
		if (strlower(trim(tb)) == "-inf" .or. strlower(trim(tb)) == "inf" .or. strlower(trim(tb)) == "nan") cycle
		if (strlower(trim(fb)) == "-inf" .or. strlower(trim(fb)) == "inf" .or. strlower(trim(fb)) == "nan") cycle
		if (strlower(trim(wb)) == "-inf" .or. strlower(trim(wb)) == "inf" .or. strlower(trim(wb)) == "nan") cycle
		
		! Read the line into the variables:
		read(tb, *) t(i)
		read(fb, *) f(i)
		if (bolUseWeights) then
			read(wb, *) w(i)
		endif
		
		i = i + 1
	enddo
	
	! Statistical Weights:
	! Save the sum of the weights,
	! which is used in the further calculations
	! (see CalcAlphaBeta)
	if (bolUseWeights) then
		if (any(w == 0)) then
			print *, "SIGMA zero values detected. SIGMA can't be zero."
			stop
		endif
		w = 1.0_dbl/(w**2)
		sumWeights = sum(w)
	else
		w = 1.0_dbl
		sumWeights = i-1
	endif
end subroutine ImportData

!subroutine sincos(x, sx, cx)
!	use kind_defs
!	implicit none
!	real(dbl), intent(in) :: x
!	real(dbl), intent(out) :: sx, cx
!	
!#if defined _ACML
!	call fastsincos(x, sx, cx)
!#else
!	sx = sin(x)
!	cx = cos(x)
!#fi
!end subroutine

! Calculates the alpha and beta-values for a given frequency.
! Input:	t - vector with times in seconds
!			f - vector with data
!			w - vector with weights
!			nu - the frequency in microHz
! Output:	alpha
!			beta
subroutine CalcAlphaBeta(t, f, w, nu, alpha, beta)
	use kind_defs
	use timeseries, only: sumWeights
	implicit none
	real(dbl), intent(in) :: t(:), f(:), w(:), nu
	real(dbl), intent(out) :: alpha, beta
	real(dbl) :: omega, s, c, cs, ss, cc, D, sx, cx, tx
	integer :: i, N

	N = size(t)
	omega = nu*twopi_micro
	
	! Calcultae sums:
	s = 0.0_dbl; c = 0.0_dbl;
	cs = 0.0_dbl; cc = 0.0_dbl;
	do i = 1,N
		sx = sin(omega*t(i))
		cx = cos(omega*t(i))
		!call sincos(omega*t(i), sx, cx)
		s = s + f(i) * sx * w(i)
		c = c + f(i) * cx * w(i)
		cs = cs + sx*cx * w(i)
		cc = cc + cx*cx * w(i)
	enddo
	! Calculate ss on basis of cc and the sum
	! of the weights, which was calculated in
	! ImportData:
	ss = sumWeights - cc
	
	! Calculate amplitude and phase:
	D = ss*cc - cs*cs
	alpha = (s * cc - c * cs)/D
	beta  = (c * ss - s * cs)/D
end subroutine CalcAlphaBeta

! Weighted Mean
function mean(f, w)
	use kind_defs
	implicit none
	real(dbl), intent(in) :: f(:)
	real(dbl), intent(in), optional :: w(:)
	real(dbl) :: mean

	if (present(w)) then
		mean = dot_product(f,w)/sum(w)
	else
		mean = sum(f)/size(f)
	endif
end function mean

! Median
function median(x)
	use kind_defs
	implicit none
	real(dbl), intent(in) :: x(:)
	real(dbl) :: median
	integer :: N, err
	integer, dimension(size(x)) :: iperm
	N = size(x)

	! Sort the copy
	call dpsort(x, N, iperm, 1, err)
	
	! Compute the median:
	if (mod(N,2) == 0) then
		median = 0.5_dbl * ( x(iperm(N/2)) + x(iperm(N/2+1)) )
	else
		median = x(iperm(N/2+1))
	end if
end function median

! Display progress for process.
subroutine Progress2(N)
	use timeseries, only: quiet
	implicit none
	integer, intent(in) :: N
	character, parameter :: CR = achar(13)
	integer, save :: j
	real, save :: nwr
	integer :: pro
	
	if (quiet) then
		return
	endif
	
	j = j+1
	if (j == 1) nwr = 0.0
	if (j >= nwr .or. j >= N) then
		! Calculate the progress:
		if (j >= N) then
			pro = 100
			j = 0
		else
			pro = nint(100*real(j)/N)
		endif
		
		! Write to screen:
		write(*, '(a1," ",i4,"%")', ADVANCE='NO') CR, pro
		
		! Set next write-point:
		if (N < 100) then
			nwr = j+1
		else
			nwr = nwr + (N*0.01)
		endif
	endif
end subroutine Progress2

! Print information about timeseries:
subroutine print_info(t, nyquist, dnu)
	use kind_defs
	use timeseries, only: median, quiet
	implicit none
	real(dbl), intent(in) :: t(:)
	real(dbl), intent(out), optional :: nyquist
	real(dbl), intent(out), optional :: dnu
	real(dbl) :: dF, max_nu
	integer, dimension(size(t)) :: iperm
	integer :: N, err
	
	N = size(t)
	call dpsort(t, N, iperm, 1, err)
	dF = 1.0e6_dbl/( t(iperm(N)) - t(iperm(1)) )
	max_nu = 0.5e6_dbl/median( t(iperm(2:N)) - t(iperm(1:N-1)) )
	
	! Write the results to screen:
	if (.not. quiet) then
		print '(" #points: ", i12)', N
		print '(" Nyquist: ", f12.6, " microHz")', max_nu
		print '(" DeltaNu: ", f12.6, " microHz")', dF
	endif
	
	! If asked for it, return the parameters:
	if (present(nyquist)) then
		nyquist = max_nu
	endif
	if (present(dnu)) then
		dnu = dF
	endif
end subroutine print_info

! Print elapsed time:
subroutine print_time(dt)
	use timeseries, only: quiet
	use kind_defs
	implicit none
	real(dbl), intent(in) :: dt
	
	if (quiet) then
		return
	endif
	
	if (dt > 3600) then
		print '(" Elapsed time: ",f7.1," hours.")', dt/3600.0
	elseif (dt > 60) then
		print '(" Elapsed time: ",f7.1," minutes.")', dt/60.0
	else
		print '(" Elapsed time: ",f7.1," seconds.")', dt
	endif
end subroutine print_time

! Calculate the Windowfunction for a given frequency:
function windowfunction(t, w, nuw, nu, fs, fc)
	use kind_defs
	use timeseries, only: sumWeights
	implicit none
	real(dbl), intent(in) :: t(:), w(:), nuw, nu
	real(dbl), intent(in) :: fs(:), fc(:)
	real(dbl) :: windowfunction
	real(dbl) :: omega, s_s, c_s, s_c, c_c
	real(dbl) :: cs, ss, cc, D, sx, cx
	real(dbl) :: alphas, betas, alphac, betac
	integer :: i, N

	! No need to calculate the central frequency:
	if (nu == nuw) then
		windowfunction = 1.0_dbl
		return
	endif

	N = size(t)
	omega = nu*twopi_micro

	! Calcultae sums:
	s_s = 0.0_dbl; c_s = 0.0_dbl; s_c = 0.0_dbl;
	c_c = 0.0_dbl; cs = 0.0_dbl; cc = 0.0_dbl;
	do i = 1,N
		sx = sin(omega*t(i))
		cx = cos(omega*t(i))
		!call fastsincos(omega*t(i), sx, cx)	
		s_s = s_s + fs(i) * sx
		c_s = c_s + fs(i) * cx
		s_c = s_c + fc(i) * sx
		c_c = c_c + fc(i) * cx
		cs = cs + sx*cx * w(i)
		cc = cc + cx*cx * w(i)
	enddo
	! Calculate ss on basis of cc and the sum
	! of the weights, which was calculated in
	! ImportData:
	ss = sumWeights - cc

	! Calculate alpha and beta for both:
	D = 1.0_dbl/(ss*cc - cs**2)
	alphas = (s_s * cc - c_s * cs)*D
	betas  = (c_s * ss - s_s * cs)*D
	alphac = (s_c * cc - c_c * cs)*D
	betac  = (c_c * ss - s_c * cs)*D
	
	! The window function is then the average:
	windowfunction = 0.5_dbl*(alphas**2+betas**2 + alphac**2+betac**2)
end function windowfunction

! Perform a bandpass-filter:
function bandpass(t, f_in, w, nu1_in, nu2, nu_min_in, nu_max, dF) result (bp)
	use kind_defs
	use timeseries, only : windowfunction, CalcAlphaBeta, mean
	implicit none
	real(dbl), intent(in) :: t(:), f_in(:), w(:), nu1_in, nu2, dF, nu_min_in, nu_max
	real(dbl), dimension(size(t)) :: bp, fs, fc
	real(dbl), allocatable :: alpha(:), beta(:), nuvec(:)
	real(dbl) :: nu, nuw, omega, s, sW, fmean, nu1, f(size(f_in)), nu_min, dnuwin, priv_sW, omegaw
	integer :: N, j, i, err, Nspec, Nband
	
	! Do so we can alter the input-variables:
	f = f_in
	nu1 = nu1_in
	nu_min = nu_min_in

	! The Algoritm is not defined for zero:
	if (nu1 <= 0) nu1 = dF
	if (nu_min <= 0) nu_min = dF
	
	nuw = (nu1+nu2)/2.0_dbl ! The frequency used for the window-function
	N = size(t)
	Nspec = ceiling((nu_max-nu_min)/dF)+1
	Nband = ceiling((nu2-nu_min)/dF) - ceiling((nu1-nu_min)/dF) + 1

	! Subtract mean from data:
	fmean = mean(f,w)
	f = f - fmean

	! Calculate windowfunction datasets:
	omegaw = nuw*twopi_micro
	fs(:) = w(:)*sin(omegaw*t(:)) ! Sine-function "dataset"
	fc(:) = w(:)*cos(omegaw*t(:)) ! Cosine-function "dataset"

	! Calculate the sum of the window-function:
	print *, 'Calcultaing window-function...'
	sW = 0.0_dbl
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(priv_sW,nu,j)
	priv_sW = 0.0_dbl
	!$OMP DO SCHEDULE(STATIC)
	do j = 1,Nspec
		nu = (j-1)*dF+nu_min
		priv_sW = priv_sW + windowfunction(t,w,nuw,nu,fs,fc)
		call Progress2(Nspec)
	enddo
	!$OMP END DO NOWAIT
	!$OMP CRITICAL
		sW = sW + priv_sW
	!$OMP END CRITICAL
	!$OMP END PARALLEL
	print *, 'Done.'
	
	! Allocate memory for alpha and betas:
	allocate(alpha(1:Nband), beta(1:Nband), nuvec(1:Nband), stat=err)
	if (err /= 0) stop "Couldn''t allocate memory!"
	
	! Calculate spectrum from 0 to cutoff:
	print *, 'Calculating spectrum...'
	!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
	do i = 1,Nband
		nuvec(i) = (i-1)*dF+nu1
		call CalcAlphaBeta(t,f,w,nuvec(i),alpha(i),beta(i))
		call Progress2(Nband)
	enddo
	!$OMP END PARALLEL DO
	print *, 'Done.'
	
	! Change to angular frequency:
	nuvec = nuvec*pi*2e-6_dbl

	! Calculate the new timeseries:
	print *, 'Calculating new timeseries...'
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,s)
	s = 0.0_dbl
	!$OMP DO SCHEDULE(STATIC)
	do j = 1,N
		! Calculate sum:
		s = sum( alpha*sin(nuvec*t(j)) + beta*cos(nuvec*t(j)) )
		! The bandpass'ed timeseries:
		bp(j) = s/sW
		! Write progress:
		call Progress2(N)
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	print *, 'Done.'
	
	! Add the mean again:
	bp = bp + fmean
	
	! Remove the memory for alpha and beta:
	deallocate(alpha, beta, nuvec)
end function bandpass

! Perform a lowpass-filter:
function lowpass(t, f, w, nu_lowpass, nu_min, nu_max, dF)
	use kind_defs
	use timeseries, only : bandpass
	implicit none
	real(dbl), intent(in) :: t(:), f(:), w(:), nu_lowpass, dF, nu_min, nu_max
	real(dbl) :: lowpass(size(t))

	lowpass = bandpass(t, f, w, 0.0_dbl, nu_lowpass, nu_min, nu_max, dF)
end function lowpass

! Perform a highpass-filter:
function highpass(t, f, w, nu_highpass, nu_min, nu_max, dF)
	use kind_defs
	use timeseries, only : bandpass
	implicit none
	real(dbl), intent(in) :: t(:), f(:), w(:), nu_highpass, dF, nu_min, nu_max
	real(dbl) :: highpass(size(t))
	
	highpass = f - bandpass(t, f, w, 0.0_dbl, nu_highpass, nu_min, nu_max, dF)
end function highpass

! Perform a bandstop-filter:
function bandstop(t, f, w, nu1, nu2, nu_min, nu_max, dF)
	use kind_defs
	use timeseries, only : bandpass
	implicit none
	real(dbl), intent(in) :: t(:), f(:), w(:), nu1, nu2, dF, nu_min, nu_max
	real(dbl) :: bandstop(size(t))
	
	bandstop = f - bandpass(t, f, w, nu1, nu2, nu_min, nu_max, dF)
end function
