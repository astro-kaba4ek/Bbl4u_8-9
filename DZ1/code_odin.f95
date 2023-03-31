program odin
	use calculation
	implicit none

	integer :: n, i, j
	real, allocatable :: a(:,:), gama(:), beta(:), alpha(:), omega(:), u(:)
	logical :: err
	character(10) :: arg 
	
	! парсим аргументов => получаем имя входного файла
	call get_command_argument(1, arg)

	open(1, file=trim(arg))
	read(1,*) n 

	allocate(a(n,3), gama(n), beta(n), alpha(n), omega(n), u(n))

	! считываем файл
	do i = 1, n 
		read(1,*) (a(i,j), j=1, 3)
	end do 
	read(1,*) (omega(i), i=1, n)
	
	close(1)

	
	! разбиваем на отдельные диагонали
	gama = a(:,1)
	beta = a(:,2)
	alpha = a(:,3)

	err = .false.
	! метод прогонки
	call tridiagonal_matrix_algorithm(n, gama, beta, alpha, omega, u, err)


	open(2, file="sol_"//trim(arg))

	! вывод в файл и на экран
	if (err) then
		write(*,*) "Определитель матрицы равен нулю. Решений бесконечно много. &
		Метод прогонки применять нельзя."
		write(2,*) "Определитель матрицы равен нулю. Решений бесконечно много. & 
		Метод прогонки применять нельзя."
	else
		do i = 1, n 
			write(*, "(1x, F10.3)") u(i)
			write(2, "(1x, F10.3)") u(i)
		end do
	end if

	close(2)


	deallocate(a, gama, beta, alpha, omega, u)

end program odin