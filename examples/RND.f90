program random_number_generator
implicit none
integer,parameter:: xSize = 20
real    x(xSize)          !initialize array for random numbers
integer x_int(xSize)      !initializr array for random integers

!try running this program without lines 8-22
integer :: n, un, istat 
integer, allocatable ::seed(:)

call random_seed(size = n)
allocate(seed(n))
open(newunit=un, file="/dev/urandom", access="stream", &
          form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
    read(un) seed
    close(un)
else
    write (*,*) "can't read urandom"
    stop 1
end if
call random_seed(put=seed)

call random_number(x)   !initialize x with random numbers 0<=x(i)<1
x_int = int(10*x)         !initialize x_int with random integers {0,..,9}
write(*,'(A)')     "real_rnd="
write(*,'(5f6.3)') x
write(*,'(A)')     "integer_rnd="
write(*,'(5i3)')   x_int
end program
