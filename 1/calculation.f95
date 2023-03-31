module calculation
	implicit none
	
contains

	subroutine tridiagonal_matrix_algorithm(n, gama, beta, alpha, omega, solution, err)
		integer :: n, i
		real, intent(in) :: gama(0:n-1), beta(0:n-1), alpha(0:n-1), omega(0:n-1)
		real, intent(out) :: solution(0:n-1)
		logical, intent(inout) :: err
		real  :: x(0:n-2), y(0:n-2), d

		! считаем определитель матрицы
		d = determinant(n, gama(1:), beta, alpha(:n-2))

		! если определитель равен 0, досрочно выходим и запоминаем это
		if (abs(d) <= 1e-8) then
			err = .true.
			return
		end if

		! обратный ход прогонки
		x(n-2) = -gama(n-1) / beta(n-1)
		y(n-2) = omega(n-1) / beta(n-1)

		do i = n-2, 1, -1
			x(i-1) = -gama(i) / (alpha(i) * x(i) + beta(i))
			y(i-1) = (omega(i) - alpha(i) * y(i)) / (alpha(i) * x(i) + beta(i))
		end do

		! прямой ход прогонки
		solution(0) = (omega(0) - alpha(0) * y(0)) / (alpha(0) * x(0) + beta(0))
		
		do i = 0, n-2
			solution(i+1) = x(i) * solution(i) + y(i)
		end do
	end subroutine tridiagonal_matrix_algorithm


	function determinant(n, c, a, b) result(d)
		integer :: n, i
		real, intent(in) :: c(n-1), a(n), b(n-1)
		real :: f(-1:n), d

		f(-1) = 0
		f(0) = 1
		f(1) = a(1)

		do i = 2, n
			f(i) = a(i) * f(i-1) - c(i-1) * b(i-1) * f(i-2)
		end do

		d = f(n)		
	end function determinant

	
end module calculation