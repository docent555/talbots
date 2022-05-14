module types
    use, intrinsic :: iso_c_binding
    use fourier

    implicit none

    namelist /param/ nz, period, lz, lx, nx, nk, nth, delta, a0_peak, xcp, alfa, g_amp, g_x0, g_x1, r0, r1, sigma, gamma, xe, c, &
        x_out, intrvl, in_type, recount, thout, central_mirror, cont, amp_only, it_todo, norm, ops, lambda, s2k
        real(c_double)   h, lz, lx, hx, hth, delta, a0_peak, imp_x0, imp_xsp, g_amp, g_x0, g_x1, r0, r1, x_out, sigma, xe, gamma, c, c3, xcp, alfa, norm, &
        kappa, lambda, kk
        integer(c_int) ::  nx, nk, nth, nz, iimp_x0, iimp_xend, ix_out, intrvl, it_todo, it_doiter, in_type, it_made = 0, it_flag =0, ie, ig0, ig1, ixe1, ixe2
    real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
    logical(c_bool) recount, thout, central_mirror, period, amp_only, cont, ops, s2k

!namelist /perenorm_param/ nz, period, lz, lx, nx, nk, nth, delta, a0_peak, xcp, alfa, g_amp, g_x0, g_x1, r0, r1, sigma, gamma, xe, c, &
!    x_out, intrvl, in_type, recount, thout, central_mirror, cont, amp_only, it_todo, norm

        complex(c_double_complex), allocatable :: a1(:), a0(:), ak1(:), ak0(:), atmp(:), jk1(:), jk0(:), k(:), ex(:), dlt(:), tmp(:), aktmp(:), akzl(:), &
                                              k2(:), akz0(:), a0z0(:), a0z0cut(:)
    real(c_double), allocatable :: th0(:, :), th1(:, :), dthdz(:, :), fk1(:), fk2(:), rhs0(:, :), z(:), x(:), &
                                   a_amp_z0(:), a_amp_zl(:), a_spec_amp_z0(:), a_spec_amp_zl(:), g(:), &
                                   sum_abs2_a_plus_by_z(:), sum_abs2_a_plus_by_z_k(:)!, theta(:,:,:)
    integer(c_int), allocatable :: it(:)

    complex(c_double_complex), parameter :: im1 = (0.0d0, 1.0d0)
    real(c_double) :: start_time, finish_time, smirr, soutm

!interface
!    subroutine fn_for_fortran_to_call(ptr) &
!        bind(c, name='fn_for_fortran_to_call')
!        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
!        implicit none
!        type(c_ptr), intent(in), value :: ptr
!    end subroutine fn_for_fortran_to_call
!end interface
end module types
