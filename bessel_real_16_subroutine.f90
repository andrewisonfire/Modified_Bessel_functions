! using some formula from DOI: 10.2307/2006491
subroutine calc_bessel_q_precision(num,arg,Iq,Kq,dIq,dKq)
    implicit none
    integer(4), parameter :: nnmax = 100
    integer(4) nn
    real(16) arrKq(0:nnmax),arrIq(0:nnmax)
    real(16) arrdKq(0:nnmax),arrdIq(0:nnmax)
    real(16) tmp,nutmp
    !---------------------------------------------------------
    integer(4) nmax,kmax,delta,itmp
    integer(4) i,k,l,n,m
    real(16) num,arg,Iq,Kq,dIq,dKq
    real(16) GaussFrac0
    real(16), parameter :: PI = acos(-1.0q0)
!    real(16) k0q,k1q,k2q
!    real(16) i0q,i1q,IF1
    real(16) ri
    real(16) resKnu


    !---------------------------------------------------------
    Ri = arg                                                    ! argument of Bessel function
    nmax = int(num,4)                                           ! The whole part of the Bessel function number
    itmp = nmax - 1
    delta =  1                                                  ! A few steps forward to get a more accurate value of
                                                                ! the Infeld function of the required number

    kmax = 200000000                                                   ! The number of terms of the infinite Gauss fraction
    !---------------------------------------------------------
    arrKq(0) = sqrt(PI/2.0q0/Ri) * exp(-Ri)                         ! array index + 1/2 = function index
    arrKq(1) = sqrt(PI/2.0q0/Ri)* (1.0q0 + 1.0q0 / Ri) * exp(-Ri)   ! The initial conditions for the recurrence formula
    !---------------------------------------------------------
    !---------------------------------------------------------
    do n = 2,nnmax                                            ! MacDonald function index from 1/2 to nnmax + 1/2
        tmp = 2.0q0 * real(n - 0.5,16)/ri                     ! Successive iterations of computing the functions of MacDonald
        arrKq(n) = arrKq(n-2) + tmp * arrKq(n-1)
    end do
    !---------------------------------------------------------
    ! Calculation of derivatives using Vatson GN "A treatise on the theory of Bessel functions" page 93 formulas
    do n = 1,nnmax-1
        arrdKq(n) = -0.5q0 * (arrKq(n-1)+arrKq(n+1))
    end do
    arrdKq(0) = arrKq(0)/2.0q0/Ri - arrKq(1)
    arrdKq(nnmax) = -( arrKq(nnmax)/2.0q0/Ri + arrKq(nnmax-1) )
    !---------------------------------------------------------
    nutmp = real(nnmax+0.5,16)                                  ! calculating nnmax + 0.5 and nnmax - 0.5 Infeld function
    tmp = arrKq(nnmax) + arrKq(nnmax-1)*GaussFrac0(nutmp,ri,kmax)
    arrIq(nnmax-1) = 1.0q0/(tmp*ri)
    arrIq(nnmax) = arrIq(nnmax-1)*GaussFrac0(nutmp,ri,kmax)     ! using Gauss continued fraction
    do n = 2,nnmax
        nn = nnmax - n ! reverse                                ! Inverse iterations of computing infeld functions
        arrIq(nn) = 1.0q0/Ri/arrKq(nn+1) - arrIq(nn+1)*arrKq(nn)/arrKq(nn+1)
    end do
    !---------------------------------------------------------
    ! Calculation of derivatives using Vatson GN "A treatise on the theory of Bessel functions" page 93 formulas
    do n = 1,nnmax-1
        arrdIq(n) = 0.5q0 * ( arrIq(n-1) + arrIq(n+1) )
    end do
    arrdIq(0) = arrIq(0)/2.0q0/Ri + arrIq(1)
    arrdIq(nnmax) = arrIq(nnmax-1) - arrIq(nnmax)/2.0q0/Ri
    !---------------------------------------------------------
        write(*,*) Ri
    do n = 0,nnmax
        write(*,*) n+0.5,arrIq(n)
    end do
    open(343,file='Bessel_test.txt')
    write(343,*) 'n', 'k(z)', 'dK(z)', 'I(z)','dI(z)'
    write(343,*) 'z = ', Ri
    do n=0,nnmax
        !write(343,*) n,arrKq(n),arrdKq(n),arrIq(n),arrdIq(n)
        write(343,*) arrKq(n),arrdKq(n),arrIq(n),arrdIq(n)
    end do

    !---------------------------------------------------------
    ! write(*,*) resKnu,i0q!,resKnu,ri
    !write(*,*)'Gauss fraction = ', GaussFrac0(num,ri,400)

end subroutine
