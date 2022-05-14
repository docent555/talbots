module fun
    use, intrinsic :: iso_c_binding
    use types

contains

    subroutine calc_idx()

        implicit none

        kk = 4.0d0*pi/lambda
        if (s2k .eq. .true.) then
            !lx = lx / sqrt(kk)

            !what depends on the new lx
            xcp = 0.5d0*lx! *norm** (1.0/6.0)
            alfa = 0.5d0*lx/100.0d0! *norm** (1.0/6.0)
            g_x0 = 0.35*lx! *norm** (1.0/6.0)
            g_x1 = 0.65*lx! *norm** (1.0/6.0)
            !xe = 0.166666666666667 * lx! *norm** (1.0/6.0)
            xe = 0.22*lx! *norm** (1.0/6.0)
            x_out = 0.5d0*lx! *norm** (1.0/6.0)

            s2k = .false.
        end if

        if (norm .ne. 0.0d0) then
            lx = lx*norm**(1.0/6.0)

            !what depends on the new lx
            xcp = 0.5d0*lx! *norm** (1.0/6.0)
            alfa = 0.5d0*lx/100.0d0! *norm** (1.0/6.0)
            g_x0 = 0.35*lx! *norm** (1.0/6.0)
            g_x1 = 0.65*lx! *norm** (1.0/6.0)
            !xe = 0.166666666666667 * lx! *norm** (1.0/6.0)
            xe = 0.22*lx! *norm** (1.0/6.0)
            x_out = 0.5d0*lx! *norm** (1.0/6.0)

            delta = delta/norm**(1.0/3.0)
            a0_peak = a0_peak/norm**(2.0/3.0)
            sigma = sigma/norm**(1.0/3.0)
            c = c/norm
            if (period == .true.) then
                lz = lz*(lx*lx)/lambda
            end if
            period = .false.
            norm = 0
        else
            if (period == .true.) then
                lz = lz*(lx*lx)/lambda
            end if
            period = .false.
        end if

        open (unit=1, file='input_fortran_real.dat')
        write (unit=1, nml=param)
        close (unit=1)

        !c3 = c * c * c
        c3 = c
        h = lz/nz
        nz = nz + 1

        hth = 2.0d0*pi/nth
        hx = lx/nx

        ixe1 = xe/hx + 1
        ixe2 = nx - ixe1 + 2
        xe = (ixe1 - 1)*hx !xe clarification

        !print *, 'ixe1 = ', ixe1
        !print *, 'ixe2 = ', ixe2
        !print *, 'xe = ', xe
        !pause

        ix_out = int(x_out/hx)
        if (ix_out <= 0) ix_out = 1
        if (ix_out > nx) ix_out = nx

        iimp_x0 = max(1, int(imp_x0/hx) + 1)
        imp_x0 = (iimp_x0 - 1)*hx !to (iimp_x0 - 1) * hx exactly equal to imp_x0
        iimp_xend = min(nx + 1, int((imp_x0 + imp_xsp)/hx) + 1) ! calculate for nx'= nx + 1
        imp_xsp = (iimp_xend - iimp_x0)*hx !for point precision
        if (iimp_xend == nx + 1) iimp_xend = iimp_xend - 1 ! last interval point not involved

        !what is the iteration
        if (it_flag == 0) it_made = 0
        it_doiter = it_made + it_todo
        !it_todo = 0
    end subroutine

    subroutine calc_theta(th, dthdz)

        implicit none

        real(c_double), intent(inout) :: th(:, :), dthdz(:, :)
        integer, dimension(size(th)) :: i
        integer ix

        i = (/1:size(th, 1)/)

        do ix = 1, 2
            th(:, ix) = hth*(i - 1)
            dthdz(:, ix) = delta
        end do

        !open(777, file='test.dat')
        !doix=1,nth
        ! write(777,*) ix-1, th(ix,1), th(ix,2)
        !enddo
        !close(777)
        !stop
    end subroutine

    function jf(th)
        implicit none
        real(c_double), intent(in) :: th(:)
        complex(c_double_complex) jf

        jf = 2.0d0/dble(nth)*sum(cdexp(-im1*th))
    end function jf

    function rhs(ak, th)
        implicit none
        complex(c_double_complex), intent(in) :: ak(:)
        real(c_double), intent(in) :: th(:, :)
        real(c_double), dimension(size(th, 1), size(th, 2)) :: rhs

        !rhs(:,1) = dreal(sum(fk1 * ak) * cdexp(im1 * th(:,1)))
        !rhs(:,2) = dreal(sum(fk2 * ak) * cdexp(im1 * th(:,2)))

        rhs(:, 1) = dreal(mysum(fk1*ak)*cdexp(im1*th(:, 1)))
        rhs(:, 2) = dreal(mysum(fk2*ak)*cdexp(im1*th(:, 2)))

        !atmp=ifs(ak)
        !rhs(:,1) = dreal(atmp(ixe1) * cdexp(im1 * th(:,1)))
        !rhs(:,2) = dreal(atmp(ixe2) * cdexp(im1 * th(:,2)))
    end function rhs

    function mysum(a)
        implicit none
        complex(c_double_complex), intent(in) :: a(:)
        complex(c_double_complex) :: mysum
        integer(c_int) i, n

        n = size(a)
        mysum = dcmplx(0)

        do i = n, 1, -1
            mysum = mysum + a(i)
        end do
    end function mysum

    function dmysum(a)
        implicit none
        real(c_double), intent(in) :: a(:)
        real(c_double) :: dmysum
        integer(c_int) i, n

        n = size(a)
        !dmysum = dcmplx(0)
        dmysum = 0.0d0

        do i = n, 1, -1
            dmysum = dmysum + a(i)
        end do
    end function dmysum

    subroutine init() bind(c, name='init')
        use, intrinsic :: iso_c_binding

        implicit none

        integer i

        !interface
        !    subroutine read_param() bind(c, name='read_param')
        !    end subroutine read_param
        !    function a0_fn_stat() result(a0_res)
        !        !import
        !        import
        !        complex(c_double_complex), dimension(nx) :: a0_res
        !    end function a0_fn_stat
        !    function fk_fn(xe) result(fk_res)
        !        use, intrinsic :: iso_c_binding, only: c_double
        !        import, only: nk
        !        real(c_double), dimension(nk) :: fk_res
        !        real(c_double) xe
        !    end function fk_fn
        !    function g_fn() result(g_res)
        !        import
        !        real(c_double), dimension(nx) :: g_res
        !    end function g_fn
        !    function k_fn() result(k)
        !    use, intrinsic :: iso_c_binding
        !        import, only: nx
        !        complex(c_double_complex), dimension(2*nx) :: k
        !    end function k_fn
        !    function k2_fn() result(k2_res)
        !    import
        !        complex(c_double_complex), dimension(nk) :: k2_res
        !    end function k2_fn
        !    function dn_fn() result(dn_res)
        !        import
        !        complex(c_double_complex), dimension(nk) :: dn_res
        !    end function dn_fn
        !end interface

        call read_param()
        call calc_idx()
        call allocate_arrays()
        call calc_zxit()
        call sincost_init(nx)
        call fft_init(2*nx)
        call dst_init(nx, lx)

        ! initial conditions for a (z = 0)
        if (cont .eq. .true.) then
            open (1, file='a0.bin', form='binary', err=101)
            read (1) a0
            close (1)
        else
            a0 = a0_fn_stat()
            a0z0 = a0
        end if

        ! open (17, file = 'test.dat')
        ! do i = 1, nx
        ! write (17, '(3f17.8)') (i-1)*hx, a0 (i)
        ! enddo
        ! close (17)
        ! stop

        ! smooth f
        fk1(:) = fk_fn(xe)
        fk2(:) = fk_fn(lx - xe)

        !smoothing function for a
        g = g_fn()

        ! k ** 2
        if (ops .eq. .false.) then
            k2 = k2_fn()
        else
            k2 = dn_fn()
        end if

        open (1, file='initag.dat')
        do i = 1, nx
            write (1, '(3e17.8,i10)') (i - 1)*hx, dreal(a0(i)), g(i)
        end do
        close (1)

        open (1, file='initfk.dat')
        do i = 1, nk
            write (1, '(2e17.8,i10)') fk1(i), fk2(i), int(k2(i))
        end do
        close (1)

        write (*, *) 'nz = ', nz
        write (*, *) 'h = ', h

        print *, 'lz = ', lz
        print *, 'lx = ', lx
        print *, 'c3 = ', c3

        !print *, size(k2), size(fk1), size(fk2), size(dlt), size(ex)
        !pause
        !stop

        return
101     stop 'error of file open.'
    end subroutine init

    subroutine finish() bind(c, name='finish')
        use fourier, only: sincost_destroy, fft_destroy

        call sincost_destroy()
        call fft_destroy()
        call deallocate_arrays()
    end subroutine finish

    subroutine write_result()
        use, intrinsic :: iso_c_binding
        import

        implicit none

        integer i, j

        call cpu_time(start_time)

        it_made = it_made - 1

        open (1, file='aend.dat', err=101)
        do i = 1, nx
            write (1, '(1p3e17.8)', err=103) (i - 1)*hx, a_amp_z0(i), a_amp_zl(i)
        end do
        close (1)

        !open(877, file = 'theta1.dat', err = 101)
        ! do i=1,nz
        ! write(877, '(e17.8,\)') (i - 1) * h
        ! do j=1,nth
        ! write(877, '(e17.8,\)') theta(j, 1, i)
        ! enddo
        ! write(877, '(/,\)')
        ! enddo
        !close(877)
        !
        !open(877, file = 'theta2.dat', err = 101)
        ! do i=1,nz
        ! write(877, '(e17.8,\)') (i - 1) * h
        ! do j=1,nth
        ! write(877, '(e17.8,\)') theta(j, 2, i)
        ! enddo
        ! write(877, '(/,\)')
        ! enddo
        !close(877)

        call cpu_time(finish_time)
        print *, 'writing time = ', finish_time - start_time, ' seconds'

        return
101     stop 'error of file open.'
102     stop 'error of file reading.'
103     stop 'error of file writing.'
    end subroutine write_result

    subroutine calc_zxit()
        implicit none

        integer i

        do i = 1, nz
            z(i) = (i - 1)*h
        end do
        do i = 1, nx
            x(i) = (i - 1)*hx
        end do
        do i = (it_made + 1), it_doiter
            it(i) = i
        end do

        open (1, file='z.dat', err=101)
        do i = 1, nz
            write (1, *, err=103) z(i)
        end do
        close (1)

        open (1, file='x.dat', err=101)
        do i = 1, nx
            write (1, *, err=103) x(i)
        end do
        close (1)

        open (1, file='k.dat', err=101)
        do i = 1, nk
            write (1, *, err=103) i
        end do
        close (1)

        !open(1, file = 'it.dat', err = 101)
        !do i=(it_made + 1),it_doiter
        ! write(1,*,err = 103) it(i)
        !enddo
        !close(1)

        return
101     stop 'error of file open.'
103     stop 'error of file writing.'
    end subroutine calc_zxit

    function a0_fn_stat() result(a0_res)
        use fourier

        import, only:nx, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, in_type, lx, central_mirror, xcp, alfa!, coeff

        implicit none

        complex(c_double_complex), dimension(nx) :: a0_res, c
        real(c_double), dimension(nx) :: a0env
        integer i, icp, ix(nx)

        if (in_type == 1) then
            !initial conditions for a (one pulse in the middle)
            if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
            a0_res(nx/2 + 2:) = 0.0d0
            do i = iimp_x0, nx/2 + 1
                a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
            end do
            a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
        elseif (in_type == 2) then
            !initial conditions for a (symmetric pulses at the edges)
            if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
            if (iimp_xend < nx) a0_res(iimp_xend + 1:nx) = 0.0d0
            do i = iimp_x0, iimp_xend
                a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
            end do
            a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
        elseif (in_type == 3) then
            !initial conditions from harmonics
            c = 0
            !seed = (/2147483562, 2147483398/)
            !seed = (/3, 2/)
            !call random_seed(size = n)
            !if (n /= 2) stop 'error of random at a0_fn_stat'
            !call random_seed (put = seed)
            !call random_number(rc)
            !c(2:10:2) = cmplx(0.1 * rc(1:5), 0.0d0)

            c(2) = 0.1
            c(3) = 0.1
            c(4) = 0.05
            c(5) = 0.05

            !print *, size(c(2:10:2))
            !do i=1,size(rc)
            ! write(*,'(a,i2,a,f6.4,a,i2,a,f6.4,a,f6.4)') 'rc(', i, ') = ', rc(i) , ' c(', 2*i, ') = ', dreal(c(2*i)), ' ', dimag(c(2*i))
            !enddo

            a0_res = ifs(c)
            a0_res = cmplx(dreal(a0_res), 0.0d0)

            !open(1, file = 'test.dat')
            !do i=1,nx
            ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i)), abs(a0_res(i))
            !enddo
            !close(1)
            !stop
        elseif (in_type == 4) then
            !test initial conditions for a (one pulse in the middle)
            do i = 1, nx
                a0_res(i) = a0_peak*dsin(1*pi/lx*(i - 1)*hx)
            end do
        elseif (in_type == 5) then
            !specialist. initial conditions for a

            c = dcmplx(0)

            !c(2) = dcmplx(0.1)
            c(3) = dcmplx(0.1)
            !c(4) = dcmplx(0.1)
            !c(6) = dcmplx(0.1)
            !c(8) = dcmplx(0.1)
            !c(10) = dcmplx(0.1)
            !c(12) = dcmplx(0.1)
            !c(14) = dcmplx(0.1)

            a0_res = ifs(c)

            !open(1, file = 'test.dat')
            !do i=1,nx
            ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
            !enddo
            !close(1)
            !stop
        elseif (in_type == 6) then
            !specialist. initial conditions for a
            !initial conditions for a (symmetric pulses at the edges)

            if (central_mirror == .false.) then
                icp = xcp/hx + 1
                iimp_xend = 2*icp - 1
                ix = 0; ix = (/0:nx/2 - 1/)

                a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
                a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)
            else
                icp = xcp/hx + 1
                iimp_xend = 2*icp - 1
                ix = 0; ix = (/0:nx - 1/)

                a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
            end if

            !open(1, file = 'test.dat')
            !do i=1,nx
            ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
            ! !write(1, *) ix(i)
            !enddo
            !close(1)
            !stop
        elseif (in_type == 7) then
            !specialist. initial conditions for a
            !initial conditions for a (symmetric pulses at the edges)

            c = dcmplx(0)

            !c(2) = dcmplx(0.1)
            c(4) = dcmplx(0.1)
            c(6) = dcmplx(0.1)
            c(8) = dcmplx(0.1)
            c(10) = dcmplx(0.1)
!c(12) = dcmplx(0.1)
            !c(14) = dcmplx(0.1)

            if (central_mirror == .false.) then
                icp = xcp/hx + 1
                iimp_xend = 2*icp - 1
                ix = 0; ix = (/0:nx/2 - 1/)

                a0env(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
                a0env(nx/2 + 2:nx) = a0env(nx/2:2:-1)
            else
                ix = (/1:nx/) - 1

                a0env = dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
            end if

            a0_res = ifs(c)*a0env

            !open(1, file = 'test.dat')
            !do i=1,nx
            ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
            ! !write(1, *) ix(i)
            !enddo
            !close(1)
            !stop
        else
            print *, 'error: wrong in_type'
            pause
            stop
        end if

    end function a0_fn_stat

    function fk_fn(xe) result(fk_res)
        use, intrinsic :: iso_c_binding, only: c_double, c_int
        import, only:nk, pi, lx

        implicit none

        real(c_double) :: fk_res(nk), xe
        integer(c_int) n(nk)

        n = (/0:nk - 1/)

        fk_res = dsin(pi*n*xe/lx)

        !open(1, file = 'test.dat')
        !do i=1,nx
        ! write(1,'(i,e17.8)') i-1, fk_res(i)
        !enddo
        ! close (1)
        ! stop
    end function fk_fn

! function g_fn () result (g_res)
! import
!
! implicit none
!
! real (c_double), dimension (nx) :: g_res
! integer i
!
! i = g_x0 / hx + 1
! ! i = g_x0 / hx
!
! if (central_mirror == .false.) then
! g_res (1: i) = 1.0d0
! g_res (nx-i+2: nx) = 1.0d0
! g_res (i+1: nx-i+1) = 0.0d0
! g_res = g_res * g_amp
! else
! g_res (1: i) = 0.0d0
! g_res (nx-i+2: nx) = 0.0d0
! g_res (i+1: nx-i+1) = 1.0d0
! g_res = g_res * g_amp
! endif
! end function g_fn

    function g_fn() result(g_res)
        import

        implicit none

        real(c_double), dimension(nx) :: g_res
        integer(c_int) i

        ig0 = g_x0/hx + 1
        ig1 = g_x1/hx + 2

        if (central_mirror == .false.) then
            g_res = 1.0d0
            g_res(ig0:ig1) = 0.0d0
            soutm = -1.0
            smirr = 0.0
        else
            g_res = 0.0d0
            g_res(ig0:ig1) = 1.0d0
            smirr = -1.0
            soutm = 0.0
        end if

        do i = 1, nx
            if (g_res(i) > 0.0) then
                smirr = smirr + 1.0
            else
                soutm = soutm + 1.0
            end if
        end do

        g_res = g_res*g_amp

        ! print *, 'smirr =', smirr, 'soutm =', soutm, 'ig0 =', ig0, 'ig1 =', ig1
        ! stop
        ! print *, i0, nx-i0+2, i1
        ! stop
    end function g_fn

    function k_fn() result(k)
        use, intrinsic :: iso_c_binding
        import, only:nx

        implicit none

        complex(c_double_complex), dimension(2*nx) :: k
        complex(c_double_complex) :: im1 = (0.0d0, 1.0d0)
        integer nn

        nn = 2*nx

        k = im1*(/0:nn/2 - 1, -nn/2:-1/)
    end function k_fn

    function k2_fn() result(k2_res)
        import

        implicit none

        complex(c_double_complex), dimension(nk) :: k2_res
        integer i
        real(c_double) w

        !k**2
        do i = 1, nk
            w = pi*(i - 1)/lx
            !k2_res(i) = w * w - was
            !k2_res(i) = - w * w ! become
            k2_res(i) = -w*w/kk! become
        end do

        open (1, file='k2_n.dat')
        do i = 1, nk
            write (1, '(i,2e17.8)') i, k2_res(i)
        end do
        close (1)
    end function k2_fn

    function dn_fn() result(dn_res)
        import, only:nk, c_double_complex, c_double, lambda, lx, im1, pi

        implicit none

        complex(c_double_complex), dimension(nk) :: dn_res
        complex(c_double_complex) k
        real(c_double) tmp
        integer i

        k = 2.0d0*pi/lambda

        dn_res(1) = dcmplx(1)
        do i = 1, nk
            tmp = 1.0d0 - (i - 1)*(i - 1)/4.0d0*(lambda/lx)*(lambda/lx)
            dn_res(i) = dsqrt(tmp) - 1.0d0
        end do

        dn_res = k*dn_res

        !tmp = 1.0d0 - 1.0d0 / 4.0d0 * (lambda / lx) * (lambda / lx)
        !if (tmp .ge. 0.0d0) then
        !    dn1 = dsqrt(tmp)
        !else
        !    dn1 = im1 * dsqrt(dabs(tmp))
        !endif
        !
        !do i=1,nk
        !    tmp = 1.0d0 - (i * i) / 4.0d0 * (lambda / lx) * (lambda / lx)
        !    if (tmp .ge. 0.0d0) then
        !        dn_res(i) = dsqrt(tmp) - dn1
        !    else
        !        dn_res(i) = im1 * dsqrt(dabs(tmp)) - dn1
        !    endif
        !end do

        open (1, file='delta_n.dat')
        do i = 1, nk
            write (1, '(i,2e17.8)') i, dn_res(i)
        end do
        close (1)
    end function dn_fn

    subroutine allocate_arrays()
        import

        implicit none

        integer(c_int) err_alloc

        allocate (g(nx), a1(nx), a0(nx), ak1(nk), ak0(nk), jk1(nk), jk0(nk), atmp(nx), &
                  th0(nth, 2), th1(nth, 2), dthdz(nth, 2), fk1(nk), fk2(nk), rhs0(nth, 2), z(nz), x(nx), k2(nk), &
                  a_amp_z0(nx), a_amp_zl(nx), aktmp(nx), akzl(nk), akz0(nk), a0z0(nx), a0z0cut(nx), &
                  a_spec_amp_z0(nx), a_spec_amp_zl(nx), it(it_todo), &
                  ex(nk), k(2*nk), dlt(nk), sum_abs2_a_plus_by_z(nx), sum_abs2_a_plus_by_z_k(nk), tmp(nx), &
                  !theta(nth, 2, nz), &
                  stat=err_alloc)

        if (err_alloc /= 0) then
            pause "allocation error"
            stop
        end if
    end subroutine allocate_arrays

    subroutine deallocate_arrays()
        import

        implicit none

        integer(c_int) err_dealloc

        deallocate (g, a1, a0, ak1, ak0, jk1, jk0, atmp, &
                    th0, th1, dthdz, fk1, fk2, rhs0, z, x, k2, &
                    a_amp_z0, a_amp_zl, aktmp, akzl, akz0, a0z0cut, &
                    a_spec_amp_z0, a_spec_amp_zl, it, &
                    ex, k, dlt, sum_abs2_a_plus_by_z, sum_abs2_a_plus_by_z_k, tmp, &
                    !theta, &
                    stat=err_dealloc)

        if (err_dealloc /= 0) stop "deallocation error"
    end subroutine deallocate_arrays

    subroutine read_param() bind(c, name='read_param')
        import

        implicit none

        open (unit=1, file='input_fortran.dat', status='old', err=101)
        read (unit=1, nml=param, err=102)
        close (unit=1)

        write (*, nml=param)

        return
101     print *, 'error of file open'; pause; stop
102     print *, 'error of reading file "input_fortran.dat"'; pause; stop
    end subroutine read_param

!subroutine d(dadx)
!    import only : k, nx, lx, pi, hx
!    use fourier
!    use, intrinsic :: iso_c_binding
!
!    implicit none
!
!    complex(c_double_complex), dimension(nx) :: dadx
!    complex(c_double_complex), dimension(2*nx) :: v
!
!    v(1) = (0.0d0,0.0d0)
!    v(2:nx) = dadx(2:nx)
!    v(nx+1) = (0.0d0,0.0d0)
!    v(2*nx:nx+2:-1) = -dadx(2:nx)
!
!    call fft(v)
!    v = 2.0 * pi / (2.0d0 * lx) * k * v
!    call ifft(v)
!    dadx = v(1:nx)
!
!end subroutine d

!function check_alloc()
!use, intrinsic :: iso_c_binding, only : c_bool
!    import only : a1, a0, ak1, ak0, atmp, jk1, jk0, dadx, k, ex, &
!        th0, th1, dthdz, fk1, fk2, rhs0, z, x, k2, &
!        a_amp_z0, a_amp_zl, a_spec_amp_z0, a_spec_amp_zl, g, sum_abs2_dadx_by_x
!
!    implicit none
!
!    logical(c_bool) check_alloc
!
!    !print *, allocated(a0), size(a0)
!
!    if (.not. allocated(a0)) then
!        check_alloc = .false.
!        print *, 'arrays were not allocated!'
!        pause
!        stop
!    endif
!
!    check_alloc = .true.
!end function check_alloc

    subroutine makea(atmp, ak)
        use fourier
        use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

        import, only:nk, nx

        implicit none

        complex(c_double_complex), dimension(:), intent(inout) :: atmp
        complex(c_double_complex), dimension(:), intent(in) :: ak
        integer(c_int) n1, n2

        n1 = size(ak)
        n2 = size(atmp)

        if (n1 .ne. nk .or. n2 .ne. nx) then
            print *, 'error in "makea"'
            pause
            stop
        end if

        atmp = dcmplx(0.0d0)
        atmp(1:nk) = ak

        call isint(atmp)
    end subroutine makea

    subroutine makeak(ak, atmp)
        use fourier
        use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
        import, only:nk, nx, aktmp

        implicit none

        complex(c_double_complex), dimension(:), intent(inout) :: ak
        complex(c_double_complex), dimension(:), intent(in) :: atmp
        integer(c_int) n1, n2

        n1 = size(ak)
        n2 = size(atmp)

        if (n1 .ne. nk .or. n2 .ne. nx) then
            print *, 'error in "makeak"'
            pause
            stop
        end if

        aktmp = atmp

        call sint(aktmp)
        ak = aktmp(1:nk)
    end subroutine makeak

    subroutine make_a0z0_ak1_atmp(a0z0, az0cut, ak0, ak)
        use fourier
        use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
        import, only:nk, nx, r0, g, makea

        implicit none

        interface
            !subroutine makea(atmp, ak)
            !    use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
            !    complex(c_double_complex), dimension(:), intent(inout) :: atmp
            !    complex(c_double_complex), dimension(:), intent(in) :: ak
            !end subroutine makea
            subroutine makeak(ak, atmp)
                use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
                complex(c_double_complex), dimension(:), intent(inout) :: ak
                complex(c_double_complex), dimension(:), intent(in) :: atmp
            end subroutine makeak
        end interface

        complex(c_double_complex), dimension(:), intent(inout) :: a0z0, ak0, az0cut
        complex(c_double_complex), dimension(:), intent(in) :: ak
        integer(c_int) n1, n2, n3, n4

        n1 = size(ak)
        n2 = size(az0cut)
        n3 = size(a0z0)
        n4 = size(ak0)

        if (n1 .ne. nk .or. n2 .ne. nx .or. n3 .ne. nx .or. n4 .ne. nk) then
            print *, 'error in "makea"'
            pause
            stop
        end if

        call makea(a0z0, ak) !before cutting

        az0cut = a0z0*r0*g !mirror reflection in z=0 and cut

        call sint(az0cut)

        ak0 = az0cut(1:nk)

        call makea(az0cut, ak0)
    end subroutine make_a0z0_ak1_atmp
end module fun
