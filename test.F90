program test

use CE1_CONSTS
use CE1_UTILS

 IMPLICIT NONE



real :: r1,r2,d

write(*,*) "INPUT:  R1 r2 d "

read(*,*) r1,r2,d

write(*,*) ovrcircle(r1,r2,d)


end program test
