type C_structure      !declare type
  real x
  real y
end typeC_structure

program using_types
implicit none
type (C_structure) struct

struct%x=10.1
struct%y=0.1
write(*,*) struct
end program