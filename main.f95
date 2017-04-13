! A fortran95 program for G95
! By WQY
program main
  implicit none
  integer re_i
  integer(4) n
  real(16) n1div2
  real(16) Rir,Ri,kD
  real(16) k0q,k1q,k2q
  real(16) Kq,Iq,dKq,dIq
  real(16), parameter :: PI = acos(-1.0q0)
  !---------------------------------------------------------
  n = 6                   ! The whole part of the Bessel function number
  n1div2 = real(n+0.5,16) ! number of Bessel function
  !n1div2 = qextd(n+0.5)
  kD = 1.0q2
  Rir = 1.0q-3
  Ri = Rir * kD
  !---------------------------------------------------------
  k0q = sqrt(PI/2.0q0/Ri) * exp(-Ri)
  k1q = sqrt(PI/2.0q0/Ri)* (1.0q0 + 1.0q0 / Ri) * exp(-Ri)

  write(*,*) "Bessel functions calculation greetings!"
  call calc_bessel_q_precision(n1div2,ri,Iq,Kq,dIq,dKq) !use this subroutine


!  write(*,*) 'I know PI number:'
!  write(*,*) PI
  !write(*,*) n1div2
  !re_i = system("pause")

end
