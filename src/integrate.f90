module numinteg
	interface	
		function simpson(H, FI) result (S)
			use kind_defs
			implicit none
			real(dbl), intent(in) :: H
			real(dbl), intent(in) :: FI(:)
			real(dbl) :: S
		end function
	end interface
end module

! Subroutine for integration over f(x) with the Simpson rule.
! FI: integrand f(x); H: interval; S: integral. Copyright (c) Tao Pang 1997.
function simpson(H, FI) result (S)
	use kind_defs
	implicit none
	real(dbl), intent(in) :: H
	real(dbl), intent(in) :: FI(:)
	real(dbl) :: S
	integer :: i, N
	real(dbl) :: S0,S1,S2

	N = size(FI)
	S  = 0.0_dbl
	S0 = 0.0_dbl
	S1 = 0.0_dbl
	S2 = 0.0_dbl
	do i = 2,N-1,2
		S1 = S1 + FI(i-1)
		S0 = S0 + FI(i)
		S2 = S2 + FI(i+1)
	enddo
	S = H*(S1+4.0_dbl*S0+S2)/3.0_dbl
	
	! If N is even, add the last slice separately
	if (mod(N,2) .eq. 0) then
		S = S + H*(5.0_dbl*FI(N)+8.0_dbl*FI(N-1)-FI(N-2))/12.0_dbl
	endif
end function simpson
