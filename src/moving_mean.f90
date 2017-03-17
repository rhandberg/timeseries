
program moving
	use kind_defs
	implicit none
	interface
		function moving_average(X, F) result (Y)
			use kind_defs
			implicit none
			real(dbl), intent(in) :: X(:)
			integer, intent(in) :: F
			real(dbl) :: Y(size(X))
		end function
	end interface
	integer, parameter :: N = 5000
	real(dbl) :: X(N), Y(N), ran
	integer :: i
	
	X(1) = 0
	do i = 2,N
		!call random_seed()
		call random_number(ran)
		X(i) = X(i-1) + (ran-0.5)*10
	enddo

	Y = moving_average(X, 100)
	
	open(50, file='mm.dat', status='replace', action='write')
	do i = 1,N
		write(50, '(2f20.8)') X(i), Y(i)
	enddo
	close(50)
end program

function moving_average(X, F) result (Y)
	use kind_defs
	implicit none
	real(dbl), intent(in) :: X(:)
	integer, intent(in) :: F
	real(dbl) :: Y(size(X))
	integer :: N, Wwidth, i
	real(dbl) :: cumsum

	N = size(X)
	Wwidth = 2*F + 1                     ! filter width

	! Boxcar window of length 2F+1 via recursive moving average (really fast)
	! Moving average method, except the ends:
	if (F == 0) then
		Y = X
	else
		Y(F+1) = sum(X(1:Wwidth))
		do i = F+2,N-F
			Y(i) = Y(i-1) + X(i+F) - X(i-F-1)   ! recursive moving average
		enddo
		Y = Y/Wwidth
	endif


	! Smooths the ends:
	! First end:
	!Yini = cumsum(X(1:Wwidth-2))
	!Yini = Yini(1:2:end)./(1:2:Wwidth-2)
	!Y(1:F) = Yini
	
	! Last end:
	!Yfin = cumsum(X(N:-1:N-Wwidth+3));
	!Yfin = Yfin(end:-2:1)./(Wwidth-2:-2:1)
	!Y(N-F+1:N) = Yfin
	
	! My FORTRAN-version:
	Y(1) = X(1)
	Y(2) = (X(1)+X(2))/2
	cumsum = X(1)+X(2)
	do i = 3,Wwidth-2,2
		cumsum = cumsum + X(i-1) + X(i)
		Y((i+1)/2) = cumsum / i
	enddo

	! Last end:
	Y(N) = X(N)
	Y(N-1) = (X(N)+X(N-1))/2
	cumsum = X(N)+X(N-1)
	do i = N,N-Wwidth+3,-2
		cumsum = cumsum + X(i-1) + X(i)
		Y((i+1)/2) = cumsum / i
	enddo
	
end function

