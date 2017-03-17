module FindMax
	interface
		subroutine FindSpectrumMax(t, f, w, nu_min, nu_max, alpha, beta, numax)
			use kind_defs
			implicit none
			real(dbl), intent(in) :: t(:), f(:), w(:), nu_min, nu_max
			real(dbl), intent(out) :: alpha, beta, numax
		end subroutine
	end interface
end module

subroutine FindSpectrumMax(t, f, w, nu_min_in, nu_max, alpha, beta, numax)
	use kind_defs
	use timeseries, only : CalcAlphaBeta
	implicit none
	real(dbl), intent(in) :: t(:), f(:), w(:), nu_min_in, nu_max
	real(dbl), intent(out) :: alpha, beta, numax
	real(dbl), external :: Amp
	real(dbl), parameter :: eps = 1e-5_dbl
	real(dbl) :: dF, Amax, nu, A, nu_min, nu1, nu2
	real(dbl) :: priv_max, priv_numax
	real(dbl) :: fmin
	integer :: j, Nspec

	! Set the crude spacing to 1/(4*T):
	! FIXME: Serious wast of time!
	! We are doing *2 more points!
	nu_min = nu_min_in
	dF = 0.25e6_dbl/(maxval(t)-minval(t))
	if (nu_min == 0) nu_min = dF ! Don't calculate nu=0
	Nspec = ceiling((nu_max-nu_min)/dF)+1
	if (Nspec <= 0 .or. nu_min < 0 .or. nu_max < 0) stop "Wrong frequencies!"
	
	! Calcultate spectrum:
	Amax = 0.0_dbl
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,alpha,beta,nu,A,priv_max,priv_numax)
	priv_max = 0.0_dbl
	!$OMP DO SCHEDULE(STATIC)
	do j = 1,Nspec
		nu = nu_min+dF*(j-1)
	
		! Calculate amplitude:
		call CalcAlphaBeta(t,f,w,nu,alpha,beta)
		A = alpha**2 + beta**2
        
		! Save the maximum amplitude in this part:
		if (A > priv_max) then
			priv_max = A
			priv_numax = nu
		endif
		
		! Write the progress:
		call Progress2(Nspec)
	enddo
	!$OMP END DO NOWAIT
	! Compare the maximums found by different threads:
	!$OMP FLUSH
	if (priv_max > Amax) then
		!$OMP CRITICAL
			! Save only the very highest peak:
			if (priv_max > Amax) then
				Amax = priv_max
				numax = priv_numax
			endif
		!$OMP END CRITICAL
	endif
	!$OMP END PARALLEL
	
	! Be careful -- If using old version of the PGI-compiler, you should at max use -O1 as optimization.
	! Otherwise, the above will fail, as Amax is not used later on. OR the problem can be 
	! worked around by using Amax for something:
	alpha = Amax - Amax
	
	! Minimize to find the true maximum:
	nu1 = max(nu_min, numax-dF)
	nu2 = min(nu_max, numax+dF)
	numax = fmin(nu1, nu2, Amp, eps)
	
	! Calculate the max alpha and beta corresponding to the found frequency:
	call CalcAlphaBeta(t,f,w,numax,alpha,beta)
end subroutine FindSpectrumMax

! Function to use to minimize in the CLEAN-algorithm:
function Amp(nu)
	use kind_defs
	use timeseries, only : t, f, w, CalcAlphaBeta
	implicit none
	real(dbl), intent(in) :: nu
	real(dbl) :: Amp
	real(dbl) :: alpha, beta
	
	call CalcAlphaBeta(t,f,w,nu,alpha,beta)
	Amp = -alpha**2 - beta**2
end function Amp
