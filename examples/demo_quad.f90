module my_quad
implicit none
integer, parameter :: dp = 8

contains
	real(dp) function f(x)
		implicit none
		real(dp) x
		f = exp(- x * x)
	end function f

	real(dp) function quad_d_lsquare(left, right, N)
		implicit none
		
		real(dp) left, right
		integer N
		
		integer i
		real(dp) xm, xp
		real(dp) dx, x0, tmp
		
		if (N .le. 0) then
			tmp = 0.0_dp
			quad_d_lsquare = tmp / tmp
			return
		end if
		
		dx = (right - left) / N
		x0 = left
		
		quad_d_lsquare = 0
		do i = 1, N
			xm = x0 + (i-1) * dx
			xp = x0 + i * dx
			quad_d_lsquare = quad_d_lsquare + f(xm) * (xp - xm)
		end do
	end function quad_d_lsquare

	real(dp) function quad_d_rsquare(left, right, N)
		implicit none
		
		real(dp) left, right
		integer N
		
		integer i
		real(dp) xm, xp
		real(dp) dx, x0, tmp
		
		if (N .le. 0) then
			tmp = 0.0_dp
			quad_d_rsquare = tmp / tmp
			return
		end if
		
		dx = (right - left) / N
		x0 = left
		
		quad_d_rsquare = 0
		do i = 1, N
			xm = x0 + (i-1) * dx
			xp = x0 + i * dx
			quad_d_rsquare = quad_d_rsquare + f(xp) * (xp - xm)
		end do
	end function quad_d_rsquare

	real(dp) function quad_d_trapz(left, right, N)
		implicit none
		
		real(dp) left, right
		integer N
		
		integer i
		real(dp) xm, xp
		real(dp) dx, x0, tmp
		
		if (N .le. 0) then
			tmp = 0.0_dp
			quad_d_trapz = tmp / tmp
			return
		end if
		
		dx = (right - left) / N
		x0 = left
		
		quad_d_trapz = 0
		do i = 1, N
			xm = x0 + (i-1) * dx
			xp = x0 + i * dx
			quad_d_trapz = quad_d_trapz + 0.5_dp * (f(xm) + f(xp)) * (xp - xm)
		end do
	end function quad_d_trapz
	
	real(dp) function quad_d_trapz_kah(left, right, N)
		implicit none
		
		real(dp) left, right
		integer N
		
		integer i
		real(dp) xm, xp
		real(dp) dx, x0, tmp
		real(dp) c, t, y
		
		if (N .le. 0) then
			tmp = 0.0_dp
			quad_d_trapz_kah = tmp / tmp
			return
		end if
		
		dx = (right - left) / N
		x0 = left
		
		quad_d_trapz_kah = 0
		c = 0.0

		do i = 1, N
			xm = x0 + (i-1) * dx
			xp = x0 + i * dx
			y = 0.5_dp * (f(xm) + f(xp)) * (xp - xm) - c
			t = quad_d_trapz_kah + y
			c = (t - quad_d_trapz_kah) - y
			quad_d_trapz_kah = t
		end do
	end function quad_d_trapz_kah
	
	real(dp) function quad_simps_kah(left, right, N)
		implicit none
		
		real(dp) left, right
		integer N
		
		integer i
		real(dp) xm, xp
		real(dp) dx, x0, tmp
		real(dp) c, t, y
		
		if (N .le. 0) then
			tmp = 0.0_dp
			quad_simps_kah = tmp / tmp
			return
		end if
		
		dx = (right - left) / N
		x0 = left
		
		quad_simps_kah = 0
		c = 0.0

		do i = 1, N
			xm = x0 + (i-1) * dx
			xp = x0 + i * dx
			y = (f(xm) + 4* f(0.5_dp * (xp + xm)) + f(xp)) * (xp - xm)/6 - c
			t = quad_simps_kah + y
			c = (t - quad_simps_kah) - y
			quad_simps_kah = t
		end do
	end function quad_simps_kah
end module my_quad


program integralchik_pryamoygolnik
use my_quad
implicit none
	
	integer N, i
	real(dp) left, right, quad
	real(8), parameter :: precise = sqrt(4.0_8 * atan(1.0_8))
	real(dp), parameter :: dpow = 0.5;
	integer, parameter :: max_pow = int(20 / dpow)
	real(dp), parameter :: base = 2.0_dp
	real(dp) pwr
	
	left = 0.0_dp
	right = 15.0_dp
	write (*,*) "precise", precise
	
	write (*,*) "left squares 0...10"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = 2 * quad_d_lsquare (left, right, N)		
		write (*,*) quad, N, abs(quad - precise)
	end do
	
	write (*,*) char(10), "left squares 10...0"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = -2 * quad_d_rsquare (right, left, N)		
		write (*,*) quad, N, abs(quad - precise)
	end do
	
	write (*,*) char(10), "right squares 0...10"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = 2 * quad_d_rsquare (left, right, N)		
		write (*,*) quad, N, abs(quad - precise)
	end do
	
	write (*,*) char(10), "right squares 10...0"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = -2 * quad_d_lsquare (right, left, N)		
		write (*,*) quad, N, abs(quad - precise)
	end do
	
	write (*,*) char(10), "trapz 0...10"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = 2 * quad_d_trapz (left, right, N)		
		write (*,*) quad, pwr, abs(quad - precise)
	end do
	
	write (*,*) char(10), "trapz 10...0"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = -2 * quad_d_trapz (right, left, N)		
		write (*,*) quad, pwr, abs(quad - precise)
	end do

	write (*,*) char(10), "trapz_kah 0...10"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = 2 * quad_d_trapz_kah (left, right, N)		
		write (*,*) quad, pwr, abs(quad - precise)
	end do
	
	write (*,*) char(10), "trapz_kah 10...0"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = -2 * quad_d_trapz_kah (right, left, N)		
		write (*,*) quad, pwr, abs(quad - precise)
	end do
	
	write (*,*) char(10), "simps_kah 0...10"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = 2 * quad_simps_kah (left, right, N)		
		write (*,*) quad, pwr, abs(quad - precise)
	end do
	
	write (*,*) char(10), "simps_kah 10...0"
	do i = 1, max_pow
		pwr = (i*dpow)
		N = int(base**pwr)
		quad = -2 * quad_simps_kah (right, left, N)		
		write (*,*) quad, pwr, abs(quad - precise)
	end do
end program 
