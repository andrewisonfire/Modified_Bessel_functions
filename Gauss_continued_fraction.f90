real(16) function a0(nu,z)
    implicit none
    real(16) nu,z

    a0 = 2.0q0 * nu / z

    return
end function

real(16) function ak(nu,z,k)
    implicit none
    integer(4) k
    real(16) nu,z,kk

if (k.le.0)then

    write(*,*) 'Bad index in Gauss continued fraction'

    stop

    else

    kk = real(k,16)
    ak =z**2 / 4.0q0 / (nu + kk - 1.0q0) / (nu + kk)

end if


    return
end function

real(16) function GaussFrac0(nu,z,kmax)
    implicit none
    integer(4) k,kmax
    real(16) nu,z
    real(16) a0,ak,tmp

    tmp = ak(nu,z,kmax) + 1.0q0

    do k = kmax-1,1
        tmp = 1.0q0 + ak(k) / tmp
    end do

    GaussFrac0 = 1.0q0 / tmp / a0(nu,z)

    return
end function
