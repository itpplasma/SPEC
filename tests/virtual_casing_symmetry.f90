! An independently prescribed Fourier field on even and odd periodic grids.
program test_virtual_casing_symmetry
    use constants, only: pi2
    implicit none
    integer :: Nt, Nz, Nzwithsym, jj, kk, index
    real :: theta, zeta
    real, allocatable :: expected(:), field(:)

    do Nt = 7, 8
        do Nz = 1, 8
            allocate(expected(Nt*Nz), field(Nt*Nz))
            do kk = 0, Nz-1
                zeta = pi2*kk/Nz
                do jj = 0, Nt-1
                    theta = pi2*jj/Nt
                    index = 1 + jj + kk*Nt
                    expected(index) = sin(theta)*cos(zeta) + 0.3*cos(2*theta)*sin(zeta)
                enddo
            enddo
            Nzwithsym = min(Nz, (Nz+3)/2)
            ! Poison every sample that bnorml does not evaluate.
            field = huge(1.0)
            if (Nz.eq.1) then
                field(1:(Nt+1)/2+1) = expected(1:(Nt+1)/2+1)
            else
                field(1:Nt*Nzwithsym) = expected(1:Nt*Nzwithsym)
            endif
            call mirror_normal_field(Nt, Nz, Nzwithsym, field)
            if (maxval(abs(field-expected)).gt.2.e-14) error stop 'periodic symmetry mismatch'
            deallocate(expected, field)
        enddo
    enddo
end program test_virtual_casing_symmetry
