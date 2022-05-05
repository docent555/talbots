MODULE IC
USE, INTRINSIC :: ISO_C_BINDING
USE FOURIER

IMPLICIT NONE

NAMELIST /PARAM/ NZ, PERIOD, LZ, LX, NX, NK, NTH, DELTA, A0_PEAK, XCP, ALFA, G_AMP, G_X0, G_X1, R0, R1, SIGMA, GAMMA, XE, C, &
    X_OUT, INTRVL, IN_TYPE, RECOUNT, THOUT, CENTRAL_MIRROR, CONT, AMP_ONLY, IT_TODO, NORM, OPS, LAMBDA, S2K
        REAL(C_DOUBLE)   H, LZ, LX, HX, HTH, DELTA, A0_PEAK, IMP_X0, IMP_XSP, G_AMP, G_X0, G_X1, R0, R1, X_OUT, SIGMA, XE, GAMMA, C, C3, XCP, ALFA, NORM, &
    KAPPA, LAMBDA, KK
        INTEGER(C_INT) ::  NX, NK, NTH, NZ, IIMP_X0, IIMP_XEND, IX_OUT, INTRVL, IT_TODO, IT_DOITER, IN_TYPE, IT_MADE = 0, IT_FLAG =0, IE, IG0, IG1, IXE1, IXE2
REAL(C_DOUBLE), PARAMETER :: PI = 2.0D0*DACOS(0.0D0)
LOGICAL(C_BOOL) RECOUNT, THOUT, CENTRAL_MIRROR, PERIOD, AMP_ONLY, CONT, OPS, S2K

!NAMELIST /PERENORM_PARAM/ NZ, PERIOD, LZ, LX, NX, NK, NTH, DELTA, A0_PEAK, XCP, ALFA, G_AMP, G_X0, G_X1, R0, R1, SIGMA, GAMMA, XE, C, &
!    X_OUT, INTRVL, IN_TYPE, RECOUNT, THOUT, CENTRAL_MIRROR, CONT, AMP_ONLY, IT_TODO, NORM

        COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: A1(:), A0(:), AK1(:), AK0(:), ATMP(:), JK1(:), JK0(:), K(:), EX(:), DLT(:), TMP(:), AKTMP(:), AKZL(:), &
                                          K2(:), AKZ0(:), A0Z0(:), A0Z0CUT(:)
REAL(C_DOUBLE), ALLOCATABLE :: TH0(:, :), TH1(:, :), DTHDZ(:, :), FK1(:), FK2(:), RHS0(:, :), Z(:), X(:), &
                               A_AMP_Z0(:), A_AMP_ZL(:), A_SPEC_AMP_Z0(:), A_SPEC_AMP_ZL(:), G(:), &
                               SUM_ABS2_A_PLUS_BY_Z(:), SUM_ABS2_A_PLUS_BY_Z_K(:)!, THETA(:,:,:)
INTEGER(C_INT), ALLOCATABLE :: IT(:)

COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: IM1 = (0.0D0, 1.0D0)
REAL(C_DOUBLE) :: START_TIME, FINISH_TIME, SMIRR, SOUTM

!INTERFACE
!    SUBROUTINE FN_FOR_FORTRAN_TO_CALL(PTR) &
!        BIND(C, NAME='FN_FOR_FORTRAN_TO_CALL')
!        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT
!        IMPLICIT NONE
!        TYPE(C_PTR), INTENT(IN), VALUE :: PTR
!    END SUBROUTINE FN_FOR_FORTRAN_TO_CALL
!END INTERFACE

CONTAINS
SUBROUTINE CALC_IDX()

    IMPLICIT NONE

    KK = 4.0D0*PI/LAMBDA
    IF (S2K .EQ. .TRUE.) THEN
        !LX = LX / SQRT(KK)

        !ТО, ЧТО ЗАВИСИТ ОТ НОВОГО LX
        XCP = 0.5D0*LX! * NORM ** (1.0/6.0)
        ALFA = 0.5D0*LX/100.0D0! * NORM ** (1.0/6.0)
        G_X0 = 0.35*LX! * NORM ** (1.0/6.0)
        G_X1 = 0.65*LX! * NORM ** (1.0/6.0)
        !XE = 0.166666666666667 * LX! * NORM ** (1.0/6.0)
        XE = 0.22*LX! * NORM ** (1.0/6.0)
        X_OUT = 0.5D0*LX! * NORM ** (1.0/6.0)

        S2K = .FALSE.
    END IF

    IF (NORM .NE. 0.0D0) THEN
        LX = LX*NORM**(1.0/6.0)

        !ТО, ЧТО ЗАВИСИТ ОТ НОВОГО LX
        XCP = 0.5D0*LX! * NORM ** (1.0/6.0)
        ALFA = 0.5D0*LX/100.0D0! * NORM ** (1.0/6.0)
        G_X0 = 0.35*LX! * NORM ** (1.0/6.0)
        G_X1 = 0.65*LX! * NORM ** (1.0/6.0)
        !XE = 0.166666666666667 * LX! * NORM ** (1.0/6.0)
        XE = 0.22*LX! * NORM ** (1.0/6.0)
        X_OUT = 0.5D0*LX! * NORM ** (1.0/6.0)

        DELTA = DELTA/NORM**(1.0/3.0)
        A0_PEAK = A0_PEAK/NORM**(2.0/3.0)
        SIGMA = SIGMA/NORM**(1.0/3.0)
        C = C/NORM
        IF (PERIOD == .TRUE.) THEN
            LZ = LZ*(LX*LX)/LAMBDA
        END IF
        PERIOD = .FALSE.
        NORM = 0
    ELSE
        IF (PERIOD == .TRUE.) THEN
            LZ = LZ*(LX*LX)/LAMBDA
        END IF
        PERIOD = .FALSE.
    END IF

    OPEN (UNIT=1, FILE='input_fortran_REAL.DAT')
    WRITE (UNIT=1, NML=PARAM)
    CLOSE (UNIT=1)

    !C3 = C * C * C
    C3 = C
    H = LZ/NZ
    NZ = NZ + 1

    HTH = 2.0D0*PI/NTH
    HX = LX/NX

    IXE1 = XE/HX + 1
    IXE2 = NX - IXE1 + 2
    XE = (IXE1 - 1)*HX !УТОЧНЕНИЕ XE

    !PRINT *, 'IXE1 = ', IXE1
    !PRINT *, 'IXE2 = ', IXE2
    !PRINT *, 'XE = ', XE
    !PAUSE

    IX_OUT = INT(X_OUT/HX)
    IF (IX_OUT <= 0) IX_OUT = 1
    IF (IX_OUT > NX) IX_OUT = NX

    IIMP_X0 = MAX(1, INT(IMP_X0/HX) + 1)
    IMP_X0 = (IIMP_X0 - 1)*HX !ЧТОБЫ (IIMP_X0 - 1) * HX ТОЧНО РАВНЯЛОСЬ IMP_X0
    IIMP_XEND = MIN(NX + 1, INT((IMP_X0 + IMP_XSP)/HX) + 1) ! РАССЧИТЫВАЕМ ДЛЯ NX'= NX + 1
    IMP_XSP = (IIMP_XEND - IIMP_X0)*HX !ДЛЯ ТОЧНОСТИ ПОПАДАНИЯ В ТОЧКИ
    IF (IIMP_XEND == NX + 1) IIMP_XEND = IIMP_XEND - 1 ! ПОСЛЕДНЯЯ ТОЧКА ИНТЕРВАЛА НЕ УЧАСТВУЕТ

    !КАКАЯ ПО СЧЕТУ ИТЕРАЦИЯ
    IF (IT_FLAG == 0) IT_MADE = 0
    IT_DOITER = IT_MADE + IT_TODO
    !IT_TODO = 0
END SUBROUTINE

SUBROUTINE CALC_THETA(TH, DTHDZ)

    IMPLICIT NONE

    REAL(C_DOUBLE), INTENT(INOUT) :: TH(:, :), DTHDZ(:, :)
    INTEGER, DIMENSION(SIZE(TH)) :: I
    INTEGER IX

    I = (/1:SIZE(TH, 1)/)

    DO IX = 1, 2
        TH(:, IX) = HTH*(I - 1)
        DTHDZ(:, IX) = DELTA
    END DO

    !OPEN(777, FILE='TEST.DAT')
    !DO IX=1,NTH
    !    WRITE(777,*) IX-1, TH(IX,1), TH(IX,2)
    !ENDDO
    !CLOSE(777)
    !STOP
END SUBROUTINE
END MODULE IC

PROGRAM ELEKTRON2DSI
    USE IC; USE FOURIER

    IMPLICIT NONE

    INTERFACE
        SUBROUTINE INIT() BIND(C, NAME='INIT')
        END SUBROUTINE INIT
        SUBROUTINE FINISH() BIND(C, NAME='FINISH')
        END SUBROUTINE FINISH
        SUBROUTINE MAKEA(ATMP, AK)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: ATMP
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: AK
        END SUBROUTINE MAKEA
        SUBROUTINE MAKEAK(AK, ATMP)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: AK
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: ATMP
        END SUBROUTINE MAKEAK
        SUBROUTINE MAKE_A0Z0_AK1_ATMP(A0Z0, AK1, ATMP, AK)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: A0Z0, AK1, ATMP
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: AK
        END SUBROUTINE MAKE_A0Z0_AK1_ATMP
    END INTERFACE

    !TYPE(C_PTR), INTENT(IN), VALUE :: PTR
    REAL(C_DOUBLE) EFF_TMP, EFF_TMP_K, EFF_TMP_B, EFF_TMP_K_B, EFF(2), &
        SUM_EFF, LOSS_ON_THE_WAY_PLUS, LOSS_ON_THE_WAY_PLUS_K, LOSS_ON_THE_WAY_MINUS, &
     INT_ABS2_A_PLUS_AT_Z0, INT_ABS2_A_PLUS_AT_Z0_K, INT_ABS2_A_PLUS_AT_ZL, INT_ABS2_A_PLUS_AT_ZL_K, INT_ABS2_A_PLUS_AT_ZL_ON_MIR, &
        INT_ABS2_A_PLUS_AT_ZL_OUT_MIR, INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K, INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K, &
            INT_ABS2_A_MINUS_AT_Z0_K, INT_ABS2_A_MINUS_AT_Z0_ON_MIR, INT_ABS2_A_MINUS_AT_Z0_OUT_MIR, INT_ABS2_A_MINUS_AT_ZL_ON_MIR, LOSS_ON_THE_WAY_MINUS_K, &
        INT_ABS2_A_MINUS_AT_ZL_ON_MIR_K, INT_OBREZKI_Z0, INT_OBREZKI_ZL, INT_ABS2_A_MINUS_AT_Z0, INT_ABS2_A_MINUS_AT_ZL, OBREZKI
    INTEGER(C_INT) IZ, PERCENT, REC_LENGTH, NR, FIRST_IT, ITH
    CHARACTER*6 STR
    INTEGER I

    !CHK_ARR = CHECK_ALLOC()
    !IF (CHK_ARR /= .TRUE.) STOP 'ALLOCATION ERROR!'

    !СЧИТЫВАЕМ ПАРАМЕТРЫ, ВЫЧИСЛЯЕМ A0(X), F(X), K2, Z, X И Т.Д.
    CALL INIT()

    !ТО, ЧТО НЕ ЗАВИСИТ ОТ J
    DLT = IM1*GAMMA*K2 + SIGMA
    EX = CDEXP(-DLT*H)

    !OPEN(17, FILE='TEST.DAT')
    !DO I=1,NK
    !    WRITE(17, '(I, 2F17.8)') I, DLT(I)
    !ENDDO
    !CLOSE(17)
    !PRINT *, SIZE(DLT), SIZE(EX), NK, SIZE(K2)
    !PAUSE
    !STOP

    REC_LENGTH = C_DOUBLE_COMPLEX*2*NX/4
    !ОТКРЫВАЕМ ФАЙЛЫ
    OPEN (221, FILE='$$$z0.bin', ACCESS='DIRECT', RECL=REC_LENGTH, ERR=101)
    OPEN (222, FILE='$$$zl.bin', ACCESS='DIRECT', RECL=REC_LENGTH, ERR=101)
    OPEN (3, FILE='eff.dat', ERR=101)
    OPEN (35, FILE='eff_b.dat', ERR=101)
    OPEN (53, FILE='eff_new.dat', ERR=101)
    OPEN (33, FILE='eff_k.dat', ERR=101)
    OPEN (335, FILE='eff_k_b.dat', ERR=101)
    OPEN (533, FILE='eff_new_k.dat', ERR=101)
    OPEN (788, FILE='it.dat', ERR=101)
    OPEN (777, FILE='for_graphics.dat', ERR=101)

    CALL CPU_TIME(START_TIME)

    FIRST_IT = IT_MADE + 1
    DO_ITER: DO IT_MADE = FIRST_IT, IT_DOITER
        PERCENT = 0
        WRITE (*, '(A11,\,I5,\,A2,I5,A,\)') 'ITERATION #', IT_MADE, ': ', PERCENT, '%'

        !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ ФАЗЫ И ЕЁ РАССТРОЙКИ ПРИ Z=0
        CALL CALC_THETA(TH0, DTHDZ)

        !THETA(:, :, 1) = TH0 !НАЧАЛЬНОЕ ПОЛОЖЕНИЕ ЭЛЕКТРОНОВ

        !РАСЧЕТ КПД В ТОЧКЕ Z=0
        !DO IX=1,NX
        !    EFF(IX) = 1.0D0 / NTH * SUM(DTHDZ(:,IX) - DELTA) !СЧИТАЕМ КПД
        !END DO
        EFF(1) = 1.0D0/NTH*SUM(DTHDZ(:, 1) - DELTA) !(XCP1)
        EFF(2) = 1.0D0/NTH*SUM(DTHDZ(:, 2) - DELTA) !(XCP2)

        !НАЧАЛЬНЫЙ ТОК (Z=0)
        JK0 = FK1*JF(TH0(:, 1)) + FK2*JF(TH0(:, 2))

        !АК0
        !AK0 = FS(A0)
        CALL MAKEAK(AK0, A0)

        !УСТАНАВЛИВАЕМ ФАЙЛОВЫЙ УКАЗАТЕЛЬ
        INQUIRE (222, NEXTREC=NR)
        !ЗАПИСЫВАЕМ A(Z=0,X)
        WRITE (221, REC=NR) A0
        !PRINT *, 'NR = ', NR

        IF (((INTRVL > 0) .AND. (MOD(IT_MADE, INTRVL) == 0)) &
            .OR. (IT_MADE == FIRST_IT) &
            .OR. (IT_MADE == IT_DOITER)) THEN
            A_AMP_Z0 = CDABS(A0)
            A_SPEC_AMP_Z0 = CDABS(FS(A0))
            WRITE (STR, 105) IT_MADE
105         FORMAT(I0.6)
            STR = TRIM(STR)
            OPEN (1, FILE='a_'//STR//'.bin', FORM='BINARY', ERR=101)
            OPEN (2, FILE='ak_'//STR//'.bin', FORM='BINARY', ERR=101)
            OPEN (88, FILE='jk_'//STR//'.bin', FORM='BINARY', ERR=101)
            WRITE (1, ERR=103) A0Z0
            WRITE (2, ERR=103) AK0
            WRITE (88, ERR=103) JK0
            IF (AMP_ONLY == .FALSE.) THEN
                OPEN (8, FILE='eff_'//STR//'.bin', FORM='BINARY', ERR=101)
                WRITE (8, ERR=103) EFF
            END IF
            IF (THOUT == .TRUE.) THEN
                OPEN (10, FILE='th_'//STR//'.dat', ERR=101)
                WRITE (10, '(E17.8,\)', ERR=103) 0.0
                DO ITH = 1, NTH
                    WRITE (10, '(E17.8,\)', ERR=103) TH0(ITH, 1)
                END DO
                DO ITH = 1, NTH
                    WRITE (10, '(E17.8,\)', ERR=103) TH0(ITH, 2)
                END DO
                WRITE (10, '(/,\)')
            END IF
        END IF

        !IF (((INTRVL > 0) .AND. (MOD(IT_MADE, INTRVL) == 0)) &
        !    .OR. (IT_MADE == FIRST_IT) &
        !    .OR. (IT_MADE == IT_DOITER)) THEN
        !    ATMP = IFS(ATMP)
        !    OPEN(777, FILE='TEST_' // STR // '.DAT')
        !    DO IX=1,NX
        !        WRITE(777,'(9E27.18)') (IX-1)*HX, CDABS(A0(IX)), G(IX), ATMP(IX), ATMP(IX) * G(IX), A0(IX)
        !    ENDDO
        !    CLOSE(777)
        !ENDIF

        !ДЛЯ РАСЧЕТА КПД ВЫЧИСЛЯЕМ СУММУ ABS(A(Z=0))**2 ПО X
        INT_ABS2_A_PLUS_AT_Z0 = SUM(CDABS(A0)*CDABS(A0))*HX
        INT_ABS2_A_PLUS_AT_Z0_K = 0.5D0*SUM(CDABS(AK0)*CDABS(AK0))

        !СУММА АММПЛИТУД ПОЛЯ ПРИ Z=0 (ТОЖЕ ДЛЯ КПД)
        SUM_ABS2_A_PLUS_BY_Z = CDABS(A0)*CDABS(A0)
        SUM_ABS2_A_PLUS_BY_Z_K = 0.5D0*CDABS(AK0)*CDABS(AK0)

        DO_Z: DO IZ = 1, NZ - 1
            RHS0 = RHS(AK0, TH0)
            TH1 = TH0 + DTHDZ*H + H/2.0D0*RHS0*H !ПРЕДИКТОР THETA
            JK1 = FK1*JF(TH1(:, 1)) + FK2*JF(TH1(:, 2)) !ПРЕДИКТОР ТОКА

            !ПРЕДИКТОР A (ИНТЕРПОЛЯЦИЯ)
            AK1(1) = AK0(1) + H/2.0D0*(JK0(1) + JK1(1))
            AK1(2:NK) = AK0(2:NK)*EX(2:NK) + &
                        C3*(JK0(2:NK) + JK1(2:NK)*(-1.0D0 + DLT(2:NK)*H) + &
                            EX(2:NK)*(JK1(2:NK) - JK0(2:NK)*(1.0D0 + DLT(2:NK)*H)))/DLT(2:NK)/DLT(2:NK)/H

            !ПРЕДИКТОР A (ТРАПЕЦИИ)
            !ATMP = (AK0 + C3 * H / 2.0D0 * JK0) * CDEXP(-DLT * H) !ЧАСТЬ A
            !AK1 = ATMP + C3 * H / 2.0D0 * JK1 !ПРЕДИКТОР A

            !A1 = IFS(AK1) !ВОЗВРАЩАЕМСЯ В РЕАЛЬНОСТЬ
            CALL MAKEA(A1, AK1)

            !КОРРЕКТОР THETA - 1
            TH1 = TH0 + DTHDZ*H + H/6.0D0*RHS0*H &
                  + H/3.0D0*RHS(AK1, TH1)*H

            !THETA(:, :, IZ + 1) = TH1 !ЗАПИШЕМ ТРАЕКТОРИЮ ЭЛЕКТРОНА

            JK1 = FK1*JF(TH1(:, 1)) + FK2*JF(TH1(:, 2)) !КОРРЕКТОР ТОКА

            !ТО, ЧТО ЗАВИСИТ ОТ ТОКА
            !JKD = JK1 - JK0

            !КОРРЕКТОР A (ИНТЕРПОЛЯЦИЯ)
            AK1(1) = AK0(1) + H/2.0D0*(JK0(1) + JK1(1))
            AK1(2:NK) = AK0(2:NK)*EX(2:NK) + &
                        C3*(JK0(2:NK) + JK1(2:NK)*(-1.0D0 + DLT(2:NK)*H) + &
                            EX(2:NK)*(JK1(2:NK) - JK0(2:NK)*(1.0D0 + DLT(2:NK)*H)))/DLT(2:NK)/DLT(2:NK)/H

            !КОРРЕКТОР A (ТРАПЕЦИИ)
            !ATMP = (AK0 + C3 * H / 2.0D0 * JK0) * CDEXP(-DLT * H) !ЧАСТЬ A
            !AK1 = ATMP + C3 * H / 2.0D0 * JK1 !КОРРЕКТОР A

            DTHDZ = DTHDZ + H/2.0D0*(RHS0 + RHS(AK1, TH1))
            !КОРРЕКТОР THETA - 1 END

            !РАСЧЕТ КПД В ТОЧКЕ Z = IZ*H
            !DO IX=1,NX
            !    EFF(IX) = 1.0D0 / NTH * SUM(DTHDZ(:,IX) - DELTA) !СЧИТАЕМ КПД
            !END DO
            EFF(1) = 1.0D0/NTH*SUM(DTHDZ(:, 1) - DELTA) !СЧИТАЕМ КПД (XCP1)
            EFF(2) = 1.0D0/NTH*SUM(DTHDZ(:, 2) - DELTA) !СЧИТАЕМ КПД (XCP2)

            !ВОЗВРАЩАЕМСЯ В РЕАЛЬНОСТЬ
            !A1 = IFS(AK1)
            CALL MAKEA(A1, AK1)

            !СУММА АММПЛИТУД ПОЛЯ ПРИ Z=IZ*H
            SUM_ABS2_A_PLUS_BY_Z = SUM_ABS2_A_PLUS_BY_Z + CDABS(A1)*CDABS(A1)
            SUM_ABS2_A_PLUS_BY_Z_K = SUM_ABS2_A_PLUS_BY_Z_K + 0.5D0*CDABS(AK1)*CDABS(AK1)

            IF (((INTRVL > 0) .AND. (MOD(IT_MADE, INTRVL) == 0)) &
                .OR. (IT_MADE == FIRST_IT) &
                .OR. (IT_MADE == IT_DOITER)) THEN
                WRITE (1, ERR=103) A1
                WRITE (2, ERR=103) AK1
                WRITE (88, ERR=103) JK1
                IF (AMP_ONLY == .FALSE.) THEN
                    WRITE (8, ERR=103) EFF
                END IF
                IF (THOUT == .TRUE.) THEN
                    WRITE (10, '(E17.8,\)', ERR=103) IZ*H
                    DO ITH = 1, NTH
                        WRITE (10, '(E17.8,\)', ERR=103) TH1(ITH, 1)
                    END DO
                    DO ITH = 1, NTH
                        WRITE (10, '(E17.8,\)', ERR=103) TH1(ITH, 2)
                    END DO
                    WRITE (10, '(/,\)')
                END IF
            END IF

            TH0 = TH1 !ДЛЯ СЛЕДУЮЩЕГО ШАГА ПО Z
            JK0 = JK1 !ДЛЯ СЛЕДУЮЩЕГО ШАГА ПО Z
            AK0 = AK1 !ДЛЯ СЛЕДУЮЩЕГО ШАГА ПО Z
            !DTHDZ0 = DTHDZ1 !ДЛЯ СЛЕДУЮЩЕГО ШАГА ПО Z

            PERCENT = INT(REAL(IZ)/REAL(NZ - 2)*100 + 0.5)
            WRITE (*, '(\,A6,I5,A)') '\B\B\B\B\B\B'C, PERCENT, '%'
        END DO DO_Z

        !РАСЧЕТ СУММЫ КПД ПО X В ТОЧКЕ Z=LZ
        SUM_EFF = (EFF(1) + EFF(2))/2.0D0

        !ДЛЯ РАСЧЕТА КПД ВЫЧИСЛЯЕМ СУММУ ABS(A(Z=LZ))**2 ПО X
        INT_ABS2_A_PLUS_AT_ZL = SUM(CDABS(A1)*CDABS(A1))*HX
        INT_ABS2_A_PLUS_AT_ZL_K = 0.5D0*SUM(CDABS(AK1)*CDABS(AK1))

        INT_ABS2_A_PLUS_AT_ZL_ON_MIR = SUM(CDABS(A1(IG0:IG1 - 1))*CDABS(A1(IG0:IG1 - 1)))*HX
    INT_ABS2_A_PLUS_AT_ZL_OUT_MIR = (SUM(CDABS(A1(1:IG0-1))*CDABS(A1(1:IG0-1))) + SUM(CDABS(A1(IG1+1:NX))*CDABS(A1(IG1+1:NX)))) * HX

        AKTMP = A1*G
        CALL SINT(AKTMP)
        INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K = INT_ABS2_A_PLUS_AT_ZL_K - 0.5D0*SUM(CDABS(AKTMP)*CDABS(AKTMP)) !ВЫЧИТАЕМ НЕОБРЕЗАННЫЕ МОДЫ ОБРЕЗАННОГО ОТРАЖЕНИЯ

        AKTMP(1:NK) = DCMPLX(0)
        INT_OBREZKI_ZL = 0.5D0*DMYSUM(CDABS(AKTMP)*CDABS(AKTMP)) ! ЭНЕРГИЯ В ОБРЕЗКАЭ НА ZL

        !ДВОЙНАЯ СУММА АМПЛИТУД ПОЛЯ ПО Z И ПО X
        LOSS_ON_THE_WAY_PLUS = 2.0D0*SIGMA*SUM(SUM_ABS2_A_PLUS_BY_Z)*HX*H
        LOSS_ON_THE_WAY_PLUS_K = 2.0D0*SIGMA*SUM(SUM_ABS2_A_PLUS_BY_Z_K)*H

        !ЗАПИСЫВАЕМ A(Z=L,X)
        WRITE (222, REC=NR) A1
        WRITE (788, *) IT_MADE

        IF (IT_MADE == IT_DOITER) THEN
            A_AMP_ZL = CDABS(A1)
            A_SPEC_AMP_ZL = CDABS(FS(A1))
        END IF

        !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ СЛЕДУЮЩЕЙ ИТЕРАЦИИ
        !IF (IT_MADE < IT_DOITER) THEN
        IF (RECOUNT == .FALSE.) THEN
            A0 = R1*A1*G*R0
        ELSE
            ATMP = R1*A1*G !!А МИНУС НА АППЕРТУРЕ Z=LZ (РЕАЛ)

            !INT_ABS2_A_MINUS_AT_ZL_ON_MIR = SUM(CDABS(ATMP(IG0:IG1-1))*CDABS(ATMP(IG0:IG1-1))) * HX
            !INT_ABS2_A_MINUS_AT_ZL_ON_MIR = SUM(CDABS(A1(IG0:IG1-1))*CDABS(A1(IG0:IG1-1))) * HX

            !ФУРЬЕ
            CALL MAKEAK(AKZL, ATMP)

            !AKTMP = DCMPLX(0)
            !AKTMP(1:NK) = AKZL
            !CALL ISINT(AKTMP)
            CALL MAKEA(ATMP, AKZL)

            INT_ABS2_A_MINUS_AT_ZL_ON_MIR_K = 0.5D0*DMYSUM(CDABS(AKZL)*CDABS(AKZL)) !ДЛЯ КПД (ТО, ЧТО УЖЕ ОТРАЗИЛОСЬ ОТ ЗЕРКАЛА)
            INT_ABS2_A_MINUS_AT_ZL_ON_MIR = SUM(CDABS(ATMP)*CDABS(ATMP))*HX

            !OPEN(777, FILE='TEST.DAT')
            !DO I=1,NX
            !    WRITE(777, '(3E17.8)') (I-1)*HX, ATMP(I)
            !ENDDO
            !CLOSE(777)
            !
            !PRINT *, INT_ABS2_A_MINUS_AT_ZL_ON_MIR_K
            !PRINT *, INT_ABS2_A_MINUS_AT_ZL
            !PAUSE
            !STOP

            !-------------------------------------------------------------------------------------------------------------------------------------------------------

            AKZ0 = AKZL*CDEXP(-DLT*LZ) !ТО, ЧТО ВЕРНУЛОСЬ В Z0

            INT_ABS2_A_MINUS_AT_Z0_K = 0.5D0*DMYSUM(CDABS(AKZ0)*CDABS(AKZ0)) !А МИНУС НА АППЕРТУРЕ Z=0 (РЕАЛ)

            CALL MAKE_A0Z0_AK1_ATMP(A0Z0, A0Z0CUT, AK0, AKZ0)

            INT_ABS2_A_MINUS_AT_Z0 = SUM(CDABS(A0Z0)*CDABS(A0Z0))*HX

            !AKZ0 МОДЫ ТОГО ЧТО ВВЕРНУЛОСЬ В НАЧАЛО
            !AK0 ОБРЕЗАННЫЕ МОДЫ ОБРЕЗАННОГО ОТРАЖЕНИЯ

            LOSS_ON_THE_WAY_MINUS_K = INT_ABS2_A_MINUS_AT_ZL_ON_MIR_K - INT_ABS2_A_MINUS_AT_Z0_K !ВЫЧИТАЕМ ТО ЧТО ВЕРНУЛОСЬ В НАЧАЛО, ОТСТАЛЬНОЕ - ПОТЕРИ

            AKTMP = A0Z0*G
            CALL SINT(AKTMP) ! НЕОБРЕЗАННЫЕ МОДЫ ОБРЕЗАННОГО ОТРАЖЕНИЯ
            INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K = INT_ABS2_A_MINUS_AT_Z0_K - 0.5D0*SUM(CDABS(AKTMP)*CDABS(AKTMP)) !ВЫЧИТАЕМ НЕОБРЕЗАННЫЕ МОДЫ ОБРЕЗАННОГО ОТРАЖЕНИЯ

            AKTMP(1:NK) = DCMPLX(0)
            INT_OBREZKI_Z0 = 0.5D0*DMYSUM(CDABS(AKTMP)*CDABS(AKTMP)) ! ЭНЕРГИЯ В ОБРЕЗКАЭ НА Z0

            INT_ABS2_A_MINUS_AT_Z0_ON_MIR = SUM(CDABS(A0Z0(IG0:IG1 - 1))*CDABS(A0Z0(IG0:IG1 - 1)))*HX
                    INT_ABS2_A_MINUS_AT_Z0_OUT_MIR = (SUM(CDABS(A0Z0(1:IG0-1))*CDABS(A0Z0(1:IG0-1))) + SUM(CDABS(A0Z0(IG1+1:NX))*CDABS(A0Z0(IG1+1:NX)))) * HX
            LOSS_ON_THE_WAY_MINUS = INT_ABS2_A_MINUS_AT_ZL_ON_MIR - (INT_ABS2_A_MINUS_AT_Z0_ON_MIR + INT_ABS2_A_MINUS_AT_Z0_OUT_MIR)

            !ОТРАЖЕНИЕ И ОБРЕЗАНИЕ У ПЕРВОГО ЗЕРКАЛА
            A0 = A0Z0CUT
            OPEN (388, FILE='a0.bin', FORM='BINARY', ERR=101)
            WRITE (388) A0
            CLOSE (388)
        END IF
        !ENDIF

        !ЗАПИСЬ КПД ПРИ ДАННОЙ ИТЕРАЦИИ
        EFF_TMP = (LOSS_ON_THE_WAY_PLUS + INT_ABS2_A_PLUS_AT_ZL - INT_ABS2_A_PLUS_AT_Z0)/LX - 4.0D0*C3*SUM_EFF
        EFF_TMP_K = LOSS_ON_THE_WAY_PLUS_K + INT_ABS2_A_PLUS_AT_ZL_K - INT_ABS2_A_PLUS_AT_Z0_K - 4.0D0*C3*SUM_EFF

        OBREZKI = INT_OBREZKI_Z0 + INT_OBREZKI_ZL

        WRITE (3, 104, ERR=103) & ! eff.dat
            IT_MADE, &
            INT_ABS2_A_PLUS_AT_ZL/LX, &
            LOSS_ON_THE_WAY_PLUS/LX, &
            INT_ABS2_A_PLUS_AT_Z0/LX, &
            SUM_EFF, &
            INT_OBREZKI_Z0, &
            INT_OBREZKI_ZL, &
            EFF_TMP, &
            EFF_TMP/(4.0D0*C3*SUM_EFF), &
            CDABS(A1(IX_OUT))

        EFF_TMP_B = (LOSS_ON_THE_WAY_MINUS + INT_ABS2_A_MINUS_AT_Z0 - INT_ABS2_A_MINUS_AT_ZL_ON_MIR)/LX

        WRITE (35, 108, ERR=103) & ! eff_b.dat
            IT_MADE, &
            LOSS_ON_THE_WAY_MINUS/LX, &
            INT_ABS2_A_MINUS_AT_ZL_ON_MIR/LX, &
            INT_ABS2_A_MINUS_AT_Z0/LX, &
            EFF_TMP_B

108     FORMAT(I, 4E17.8)

        WRITE (33, 104, ERR=103) & ! eff_k.dat
            IT_MADE, &
            INT_ABS2_A_PLUS_AT_ZL_K, &
            LOSS_ON_THE_WAY_PLUS_K, &
            INT_ABS2_A_PLUS_AT_Z0_K, &
            SUM_EFF, &
            INT_OBREZKI_Z0, &
            INT_OBREZKI_ZL, &
            EFF_TMP_K, &
            EFF_TMP_K/(4.0D0*C3*SUM_EFF), &
            CDABS(A1(IX_OUT))

104     FORMAT(I, 9E17.8)

        EFF_TMP_K_B = LOSS_ON_THE_WAY_MINUS_K + INT_ABS2_A_MINUS_AT_Z0_K - INT_ABS2_A_MINUS_AT_ZL_ON_MIR_K

        WRITE (335, 108, ERR=103) & ! eff_k_b.dat
            IT_MADE, &
            LOSS_ON_THE_WAY_MINUS_K, &
            INT_ABS2_A_MINUS_AT_ZL_ON_MIR_K, &
            INT_ABS2_A_MINUS_AT_Z0_K, &
            EFF_TMP_K_B

        WRITE (53, 107, ERR=103) IT_MADE, &   ! eff_new.dat
            INT_ABS2_A_MINUS_AT_Z0_OUT_MIR/LX, &
            INT_ABS2_A_PLUS_AT_ZL_OUT_MIR/LX, &
            (LOSS_ON_THE_WAY_MINUS_K + LOSS_ON_THE_WAY_PLUS/LX), &
            OBREZKI, &
            4.0D0*C3*SUM_EFF, &
            LOSS_ON_THE_WAY_MINUS_K + (INT_ABS2_A_PLUS_AT_ZL_OUT_MIR + &
                                       LOSS_ON_THE_WAY_PLUS + &
                                       INT_ABS2_A_MINUS_AT_Z0_OUT_MIR)/LX + OBREZKI - &
            4.0D0*C3*SUM_EFF, &
            (LOSS_ON_THE_WAY_MINUS_K + (INT_ABS2_A_PLUS_AT_ZL_OUT_MIR + &
                                        LOSS_ON_THE_WAY_PLUS + &
                                        INT_ABS2_A_MINUS_AT_Z0_OUT_MIR)/LX + OBREZKI - &
             4.0D0*C3*SUM_EFF)/4.0D0/C3/SUM_EFF

        WRITE (533, 107, ERR=103) IT_MADE, &   ! eff_new_k.dat
            INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K, &
            INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K, &
            (LOSS_ON_THE_WAY_MINUS_K + LOSS_ON_THE_WAY_PLUS_K), &
            OBREZKI, &
            4.0D0*C3*SUM_EFF, &
            LOSS_ON_THE_WAY_MINUS_K + (INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K + &
                                       LOSS_ON_THE_WAY_PLUS_K + &
                                       INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K + OBREZKI) - &
            4.0D0*C3*SUM_EFF, &
            (LOSS_ON_THE_WAY_MINUS_K + (INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K + &
                                        LOSS_ON_THE_WAY_PLUS_K + &
                                        INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K + OBREZKI) - &
             4.0D0*C3*SUM_EFF)/4.0D0/C3/SUM_EFF

107     FORMAT(I, 7E17.8)

        WRITE (777, '(8E17.8)', ERR=103) & ! for_graphics.dat
            LZ, &
            INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K, &
            INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K, &
            (LOSS_ON_THE_WAY_MINUS_K + LOSS_ON_THE_WAY_PLUS_K), &
            OBREZKI, &
            4.0D0*C3*SUM_EFF, &
            (INT_ABS2_A_MINUS_AT_Z0_OUT_MIR_K + INT_ABS2_A_PLUS_AT_ZL_OUT_MIR_K)/(4.0D0*C3*SUM_EFF), &
            DELTA

        !IX_OUT = 65
        IF (RECOUNT == .FALSE.) THEN
            WRITE (*, '(A,\)') CHAR(13)
                WRITE(*,'(A,I6,A,F7.3,A,E17.8,A,E17.8,A)') 'ITERATION #', IT_MADE, ':     A(', X_OUT, ') = ', CDABS(A1(IX_OUT)), '   EFF = ', SUM_EFF, '     CARRYOVER'
        ELSE
            WRITE (*, '(A,\)') CHAR(13)
                WRITE(*,'(A,I6,A,F7.3,A,E17.8,A,E17.8,A)') 'ITERATION #', IT_MADE, ':     A(', X_OUT, ') = ', CDABS(A1(IX_OUT)), '   EFF = ', SUM_EFF, '     RECOUNT'                                                    
        END IF

        !CALL FN_FOR_FORTRAN_TO_CALL(PTR)
    END DO DO_ITER

    WRITE (*, '(/,/)')

    !ЗАКРЫВАЕМ ФАЙЛЫ ДЛЯ ЗАПИСИ A(Z=0,X) И A(Z=L,X)
    CLOSE (1)
    CLOSE (3)
    CLOSE (53)
    CLOSE (221)
    CLOSE (222)
    CLOSE (788)
    IF (AMP_ONLY == .FALSE.) THEN
        CLOSE (2)
        CLOSE (8)
    END IF
    IF (THOUT == .TRUE.) CLOSE (10)

    CALL CPU_TIME(FINISH_TIME)
    PRINT *, 'COMPUTATION TIME = ', FINISH_TIME - START_TIME, ' SECONDS'

    CALL WRITE_RESULT()
    CALL FINISH()

    PRINT *, 'CALCULATION FINISHED.'
    PAUSE
    STOP
101 STOP 'ERROR OF FILE OPEN.'
102 STOP 'ERROR OF FILE READING.'
103 STOP 'ERROR OF FILE WRITING.'

CONTAINS
    FUNCTION JF(TH)
        IMPLICIT NONE
        REAL(C_DOUBLE), INTENT(IN) :: TH(:)
        COMPLEX(C_DOUBLE_COMPLEX) JF

        JF = 2.0D0/DBLE(NTH)*SUM(CDEXP(-IM1*TH))
    END FUNCTION JF

    FUNCTION RHS(AK, TH)
        IMPLICIT NONE
        COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: AK(:)
        REAL(C_DOUBLE), INTENT(IN) :: TH(:, :)
        REAL(C_DOUBLE), DIMENSION(SIZE(TH, 1), SIZE(TH, 2)) :: RHS

        !RHS(:,1) = DREAL(SUM(FK1 * AK) * CDEXP(IM1 * TH(:,1)))
        !RHS(:,2) = DREAL(SUM(FK2 * AK) * CDEXP(IM1 * TH(:,2)))

        RHS(:, 1) = DREAL(MYSUM(FK1*AK)*CDEXP(IM1*TH(:, 1)))
        RHS(:, 2) = DREAL(MYSUM(FK2*AK)*CDEXP(IM1*TH(:, 2)))

        !ATMP = IFS(AK)
        !RHS(:,1) = DREAL(ATMP(IXE1) * CDEXP(IM1 * TH(:,1)))
        !RHS(:,2) = DREAL(ATMP(IXE2) * CDEXP(IM1 * TH(:,2)))
    END FUNCTION RHS

    FUNCTION MYSUM(A)
        IMPLICIT NONE
        COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: A(:)
        COMPLEX(C_DOUBLE_COMPLEX) :: MYSUM
        INTEGER(C_INT) I, N

        N = SIZE(A)
        MYSUM = DCMPLX(0)

        DO I = N, 1, -1
            MYSUM = MYSUM + A(I)
        END DO
    END FUNCTION MYSUM

    FUNCTION DMYSUM(A)
        IMPLICIT NONE
        REAL(C_DOUBLE), INTENT(IN) :: A(:)
        REAL(C_DOUBLE) :: DMYSUM
        INTEGER(C_INT) I, N

        N = SIZE(A)
        !DMYSUM = DCMPLX(0)
        DMYSUM = 0.0D0

        DO I = N, 1, -1
            DMYSUM = DMYSUM + A(I)
        END DO
    END FUNCTION DMYSUM
    !END SUBROUTINE CALCULATE_FORTRAN
END PROGRAM ELEKTRON2DSI

SUBROUTINE INIT() BIND(C, NAME='INIT')
    USE, INTRINSIC :: ISO_C_BINDING
    USE IC

    IMPLICIT NONE

    INTEGER I

    INTERFACE
        SUBROUTINE READ_PARAM() BIND(C, NAME='READ_PARAM')
        END SUBROUTINE READ_PARAM
        FUNCTION A0_FN_STAT() RESULT(A0_RES)
            USE IC
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX) :: A0_RES
        END FUNCTION A0_FN_STAT
        FUNCTION FK_FN(XE) RESULT(FK_RES)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
            USE IC, ONLY: NK
            REAL(C_DOUBLE), DIMENSION(NK) :: FK_RES
            REAL(C_DOUBLE) XE
        END FUNCTION FK_FN
        FUNCTION G_FN() RESULT(G_RES)
            USE IC
            REAL(C_DOUBLE), DIMENSION(NX) :: G_RES
        END FUNCTION G_FN
        FUNCTION K_FN() RESULT(K)
            USE IC, ONLY: NX
            USE, INTRINSIC :: ISO_C_BINDING
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(2*NX) :: K
        END FUNCTION K_FN
        FUNCTION K2_FN() RESULT(K2_RES)
            USE IC
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NK) :: K2_RES
        END FUNCTION K2_FN
        FUNCTION DN_FN() RESULT(DN_RES)
            USE IC
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NK) :: DN_RES
        END FUNCTION DN_FN
    END INTERFACE

    CALL READ_PARAM()
    CALL CALC_IDX()
    CALL ALLOCATE_ARRAYS()
    CALL CALC_ZXIT()
    CALL SINCOST_INIT(NX)
    CALL FFT_INIT(2*NX)
    CALL dst_init(NX, LX)

    !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A(Z=0)
    IF (CONT .EQ. .TRUE.) THEN
        OPEN (1, FILE='a0.bin', FORM='BINARY', ERR=101)
        READ (1) A0
        CLOSE (1)
    ELSE
        A0 = A0_FN_STAT()
        A0Z0 = A0
    END IF

    !OPEN(17, FILE='TEST.DAT')
    !DO I=1,NX
    !    WRITE(17, '(3F17.8)') (I-1)*HX, A0(I)
    !ENDDO
    !CLOSE(17)
    !STOP

    !ГЛАДКАЯ F
    FK1(:) = FK_FN(XE)
    FK2(:) = FK_FN(LX - XE)

    !СГЛАЖИВАЮЩАЯ Ф-ЦИЯ ДЛЯ A
    G = G_FN()

    !K**2
    IF (OPS .EQ. .FALSE.) THEN
        K2 = K2_FN()
    ELSE
        K2 = DN_FN()
    END IF

    OPEN (1, FILE='initAG.dat')
    DO I = 1, NX
        WRITE (1, '(3E17.8,I10)') (I - 1)*HX, DREAL(A0(I)), G(I)
    END DO
    CLOSE (1)

    OPEN (1, FILE='initFK.dat')
    DO I = 1, NK
        WRITE (1, '(2E17.8,I10)') FK1(I), FK2(I), INT(K2(I))
    END DO
    CLOSE (1)

    WRITE (*, *) 'NZ = ', NZ
    WRITE (*, *) 'H = ', H

    PRINT *, 'LZ = ', LZ
    PRINT *, 'LX = ', LX
    PRINT *, 'C3 = ', C3

    !PRINT *, SIZE(K2), SIZE(FK1), SIZE(FK2), SIZE(DLT), SIZE(EX)
    !PAUSE
    !STOP

    RETURN
101 STOP 'ERROR OF FILE OPEN.'
END SUBROUTINE INIT

SUBROUTINE FINISH() BIND(C, NAME='FINISH')
    USE FOURIER, ONLY: SINCOST_DESTROY, FFT_DESTROY

    CALL SINCOST_DESTROY()
    CALL FFT_DESTROY()
    CALL DEALLOCATE_ARRAYS()
END SUBROUTINE FINISH

SUBROUTINE WRITE_RESULT()
    USE IC
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTEGER I, J

    CALL CPU_TIME(START_TIME)

    IT_MADE = IT_MADE - 1

    OPEN (1, FILE='aend.dat', ERR=101)
    DO I = 1, NX
        WRITE (1, '(1P3E17.8)', ERR=103) (I - 1)*HX, A_AMP_Z0(I), A_AMP_ZL(I)
    END DO
    CLOSE (1)

    !OPEN(877, FILE = 'theta1.dat', ERR = 101)
    !    DO I=1,NZ
    !        WRITE(877, '(E17.8,\)') (I - 1) * H
    !        DO J=1,NTH
    !            WRITE(877, '(E17.8,\)') THETA(J, 1, I)
    !        ENDDO
    !        WRITE(877, '(/,\)')
    !    ENDDO
    !CLOSE(877)
    !
    !OPEN(877, FILE = 'theta2.dat', ERR = 101)
    !    DO I=1,NZ
    !        WRITE(877, '(E17.8,\)') (I - 1) * H
    !        DO J=1,NTH
    !            WRITE(877, '(E17.8,\)') THETA(J, 2, I)
    !        ENDDO
    !        WRITE(877, '(/,\)')
    !    ENDDO
    !CLOSE(877)

    CALL CPU_TIME(FINISH_TIME)
    PRINT *, 'WRITING TIME = ', FINISH_TIME - START_TIME, ' SECONDS'

    RETURN
101 STOP 'ERROR OF FILE OPEN.'
102 STOP 'ERROR OF FILE READING.'
103 STOP 'ERROR OF FILE WRITING.'
END SUBROUTINE WRITE_RESULT

SUBROUTINE CALC_ZXIT()
    USE IC

    IMPLICIT NONE

    INTEGER I

    DO I = 1, NZ
        Z(I) = (I - 1)*H
    END DO
    DO I = 1, NX
        X(I) = (I - 1)*HX
    END DO
    DO I = (IT_MADE + 1), IT_DOITER
        IT(I) = I
    END DO

    OPEN (1, FILE='z.dat', ERR=101)
    DO I = 1, NZ
        WRITE (1, *, ERR=103) Z(I)
    END DO
    CLOSE (1)

    OPEN (1, FILE='x.dat', ERR=101)
    DO I = 1, NX
        WRITE (1, *, ERR=103) X(I)
    END DO
    CLOSE (1)

    OPEN (1, FILE='k.dat', ERR=101)
    DO I = 1, NK
        WRITE (1, *, ERR=103) I
    END DO
    CLOSE (1)

    !OPEN(1, FILE = 'it.dat', ERR = 101)
    !DO I=(IT_MADE + 1),IT_DOITER
    !    WRITE(1,*,ERR = 103) IT(I)
    !ENDDO
    !CLOSE(1)

    RETURN
101 STOP 'ERROR OF FILE OPEN.'
103 STOP 'ERROR OF FILE WRITING.'
END SUBROUTINE CALC_ZXIT

FUNCTION A0_FN_STAT() RESULT(A0_RES)
    USE IC, ONLY: NX, HX, IIMP_X0, A0_PEAK, PI, IMP_XSP, IMP_X0, IIMP_XEND, IN_TYPE, LX, CENTRAL_MIRROR, XCP, ALFA!, COEFF
    USE FOURIER

    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX) :: A0_RES, C
    REAL(C_DOUBLE), DIMENSION(NX) :: A0ENV
    INTEGER I, ICP, IX(NX)

    IF (IN_TYPE == 1) THEN
        !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A (ОДИН ИМПУЛЬС В СЕРЕДИНЕ)
        IF (IIMP_X0 > 1) A0_RES(1:IIMP_X0 - 1) = 0.0D0
        A0_RES(NX/2 + 2:) = 0.0D0
        DO I = IIMP_X0, NX/2 + 1
            A0_RES(I) = A0_PEAK*DSIN(PI/IMP_XSP*((I - 1)*HX - IMP_X0))*DSIN(PI/IMP_XSP*((I - 1)*HX - IMP_X0))
        END DO
        A0_RES(NX:NX/2 + 2:-1) = A0_RES(2:NX/2) !ОТРАЖАЕМ ИМПУЛЬС
    ELSEIF (IN_TYPE == 2) THEN
        !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A (СИММЕТРИЧНЫЕ ИМПУЛЬСЫ НА КРАЯХ)
        IF (IIMP_X0 > 1) A0_RES(1:IIMP_X0 - 1) = 0.0D0
        IF (IIMP_XEND < NX) A0_RES(IIMP_XEND + 1:NX) = 0.0D0
        DO I = IIMP_X0, IIMP_XEND
            A0_RES(I) = A0_PEAK*DSIN(PI/IMP_XSP*((I - 1)*HX - IMP_X0))*DSIN(PI/IMP_XSP*((I - 1)*HX - IMP_X0))
        END DO
        A0_RES(NX:NX/2 + 2:-1) = A0_RES(2:NX/2) !ОТРАЖАЕМ ИМПУЛЬС
    ELSEIF (IN_TYPE == 3) THEN
        !НАЧАЛЬНЫЕ УСЛОВИЯ ИЗ ГАРМОНИК
        C = 0
        !SEED = (/2147483562, 2147483398/)
        !SEED = (/3, 2/)
        !CALL RANDOM_SEED(SIZE = N)
        !IF (N /= 2) STOP 'ERROR OF RANDOM AT A0_FN_STAT'
        !CALL RANDOM_SEED (PUT = SEED)
        !CALL RANDOM_NUMBER(RC)
        !C(2:10:2) = CMPLX(0.1 * RC(1:5), 0.0D0)

        C(2) = 0.1
        C(3) = 0.1
        C(4) = 0.05
        C(5) = 0.05

        !PRINT *, SIZE(C(2:10:2))
        !DO I=1,SIZE(RC)
        !    WRITE(*,'(A,I2,A,F6.4,A,I2,A,F6.4,A,F6.4)') 'RC(', I, ') = ', RC(I), '     C(', 2*I, ') = ', DREAL(C(2*I)), '   ', DIMAG(C(2*I))
        !ENDDO

        A0_RES = IFS(C)
        A0_RES = CMPLX(DREAL(A0_RES), 0.0D0)

        !OPEN(1, FILE = 'TEST.DAT')
        !DO I=1,NX
        !    WRITE(1, '(4E17.8)') (I-1)*HX, DREAL(A0_RES(I)), DIMAG(A0_RES(I)), ABS(A0_RES(I))
        !ENDDO
        !CLOSE(1)
        !STOP
    ELSEIF (IN_TYPE == 4) THEN
        !ТЕСТОВЫЕ НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A (ОДИН ИМПУЛЬС В СЕРЕДИНЕ)
        DO I = 1, NX
            A0_RES(I) = A0_PEAK*DSIN(1*PI/LX*(I - 1)*HX)
        END DO
    ELSEIF (IN_TYPE == 5) THEN
        !СПЕЦ. НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A

        C = DCMPLX(0)

        !C(2) = DCMPLX(0.1)
        C(3) = DCMPLX(0.1)
        !C(4) = DCMPLX(0.1)
        !C(6) = DCMPLX(0.1)
        !C(8) = DCMPLX(0.1)
        !C(10) = DCMPLX(0.1)
        !C(12) = DCMPLX(0.1)
        !C(14) = DCMPLX(0.1)

        A0_RES = IFS(C)

        !OPEN(1, FILE = 'TEST.DAT')
        !DO I=1,NX
        !    WRITE(1, '(4E17.8)') (I-1)*HX, DREAL(A0_RES(I)), DIMAG(A0_RES(I))
        !ENDDO
        !CLOSE(1)
        !STOP
    ELSEIF (IN_TYPE == 6) THEN
        !СПЕЦ. НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A
        !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A (СИММЕТРИЧНЫЕ ИМПУЛЬСЫ НА КРАЯХ)

        IF (CENTRAL_MIRROR == .FALSE.) THEN
            ICP = XCP/HX + 1
            IIMP_XEND = 2*ICP - 1
            IX = 0; IX = (/0:NX/2 - 1/)

            A0_RES(1:NX/2) = A0_PEAK*DEXP(-(IX*HX - XCP)**2/ALFA) !+ DEXP(-(IX * HX - XCP**2)**2/ALFA)
            A0_RES(NX/2 + 2:NX) = A0_RES(NX/2:2:-1)
        ELSE
            ICP = XCP/HX + 1
            IIMP_XEND = 2*ICP - 1
            IX = 0; IX = (/0:NX - 1/)

            A0_RES = A0_PEAK*DEXP(-(IX*HX - XCP)**2/ALFA) !+ DEXP(-(IX * HX - XCP**2)**2/ALFA)
        END IF

        !OPEN(1, FILE = 'TEST.DAT')
        !DO I=1,NX
        !    WRITE(1, '(4E17.8)') (I-1)*HX, DREAL(A0_RES(I)), DIMAG(A0_RES(I))
        !    !WRITE(1, *) IX(I)
        !ENDDO
        !CLOSE(1)
        !STOP
    ELSEIF (IN_TYPE == 7) THEN
        !СПЕЦ. НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A
        !НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ A (СИММЕТРИЧНЫЕ ИМПУЛЬСЫ НА КРАЯХ)

        C = DCMPLX(0)

        !C(2) = DCMPLX(0.1)
        C(4) = DCMPLX(0.1)
        C(6) = DCMPLX(0.1)
        C(8) = DCMPLX(0.1)
        C(10) = DCMPLX(0.1)
        !C(12) = DCMPLX(0.1)
        !C(14) = DCMPLX(0.1)

        IF (CENTRAL_MIRROR == .FALSE.) THEN
            ICP = XCP/HX + 1
            IIMP_XEND = 2*ICP - 1
            IX = 0; IX = (/0:NX/2 - 1/)

            A0ENV(1:NX/2) = A0_PEAK*DEXP(-(IX*HX - XCP)**2/ALFA) !+ DEXP(-(IX * HX - XCP**2)**2/ALFA)
            A0ENV(NX/2 + 2:NX) = A0ENV(NX/2:2:-1)
        ELSE
            IX = (/1:NX/) - 1

            A0ENV = DEXP(-(IX*HX - XCP)**2/ALFA) !+ DEXP(-(IX * HX - XCP**2)**2/ALFA)
        END IF

        A0_RES = IFS(C)*A0ENV

        !OPEN(1, FILE = 'TEST.DAT')
        !DO I=1,NX
        !    WRITE(1, '(4E17.8)') (I-1)*HX, DREAL(A0_RES(I)), DIMAG(A0_RES(I))
        !    !WRITE(1, *) IX(I)
        !ENDDO
        !CLOSE(1)
        !STOP
    ELSE
        PRINT *, 'ERROR: WRONG IN_TYPE'
        PAUSE
        STOP
    END IF

END FUNCTION A0_FN_STAT

FUNCTION FK_FN(XE) RESULT(FK_RES)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE, C_INT
    USE IC, ONLY: NK, PI, LX

    IMPLICIT NONE

    REAL(C_DOUBLE) :: FK_RES(NK), XE
    INTEGER(C_INT) N(NK)

    N = (/0:NK - 1/)

    FK_RES = DSIN(PI*N*XE/LX)

    !OPEN(1, FILE = 'TEST.DAT')
    !DO I=1,NX
    !    WRITE(1,'(I,E17.8)') I-1, FK_RES(I)
    !ENDDO
    !CLOSE(1)
    !STOP
END FUNCTION FK_FN

!FUNCTION G_FN() RESULT(G_RES)
!    USE IC
!
!    IMPLICIT NONE
!
!    REAL(C_DOUBLE), DIMENSION(NX) :: G_RES
!    INTEGER I
!
!    I = G_X0 / HX + 1
!    !I = G_X0 / HX
!
!    IF (CENTRAL_MIRROR == .FALSE.) THEN
!        G_RES(1:I) = 1.0D0
!        G_RES(NX-I+2:NX) = 1.0D0
!        G_RES(I+1:NX-I+1) = 0.0D0
!        G_RES = G_RES * G_AMP
!    ELSE
!        G_RES(1:I) = 0.0D0
!        G_RES(NX-I+2:NX) = 0.0D0
!        G_RES(I+1:NX-I+1) = 1.0D0
!        G_RES = G_RES * G_AMP
!    ENDIF
!END FUNCTION G_FN

FUNCTION G_FN() RESULT(G_RES)
    USE IC

    IMPLICIT NONE

    REAL(C_DOUBLE), DIMENSION(NX) :: G_RES
    INTEGER(C_INT) I

    IG0 = G_X0/HX + 1
    IG1 = G_X1/HX + 2

    IF (CENTRAL_MIRROR == .FALSE.) THEN
        G_RES = 1.0D0
        G_RES(IG0:IG1) = 0.0D0
        SOUTM = -1.0
        SMIRR = 0.0
    ELSE
        G_RES = 0.0D0
        G_RES(IG0:IG1) = 1.0D0
        SMIRR = -1.0
        SOUTM = 0.0
    END IF

    DO I = 1, NX
        IF (G_RES(I) > 0.0) THEN
            SMIRR = SMIRR + 1.0
        ELSE
            SOUTM = SOUTM + 1.0
        END IF
    END DO

    G_RES = G_RES*G_AMP

    !PRINT *, 'SMIRR = ', SMIRR, 'SOUTM = ', SOUTM, 'IG0 = ', IG0, 'IG1 = ', IG1
    !STOP
    !PRINT *, I0, NX-I0+2,  I1
    !STOP
END FUNCTION G_FN

FUNCTION K_FN() RESULT(K)
    USE IC, ONLY: NX
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(2*NX) :: K
    COMPLEX(C_DOUBLE_COMPLEX) :: IM1 = (0.0D0, 1.0D0)
    INTEGER NN

    NN = 2*NX

    K = IM1*(/0:NN/2 - 1, -NN/2:-1/)
END FUNCTION K_FN

FUNCTION K2_FN() RESULT(K2_RES)
    USE IC

    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NK) :: K2_RES
    INTEGER I
    REAL(C_DOUBLE) W

    !K**2
    DO I = 1, NK
        W = PI*(I - 1)/LX
        !K2_RES(I) = W * W - БЫЛО
        !K2_RES(I) = - W * W ! СТАЛО
        K2_RES(I) = -W*W/KK! СТАЛО
    END DO

    OPEN (1, FILE='K2_N.DAT')
    DO I = 1, NK
        WRITE (1, '(I,2E17.8)') I, K2_RES(I)
    END DO
    CLOSE (1)
END FUNCTION K2_FN

FUNCTION DN_FN() RESULT(DN_RES)
    USE IC, ONLY: NK, C_DOUBLE_COMPLEX, C_DOUBLE, LAMBDA, LX, IM1, PI

    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NK) :: DN_RES
    COMPLEX(C_DOUBLE_COMPLEX) K
    REAL(C_DOUBLE) TMP
    INTEGER I

    K = 2.0D0*PI/LAMBDA

    DN_RES(1) = DCMPLX(1)
    DO I = 1, NK
        TMP = 1.0D0 - (I - 1)*(I - 1)/4.0D0*(LAMBDA/LX)*(LAMBDA/LX)
        DN_RES(I) = DSQRT(TMP) - 1.0D0
    END DO

    DN_RES = K*DN_RES

    !TMP = 1.0D0 - 1.0D0 / 4.0D0 * (LAMBDA / LX) * (LAMBDA / LX)
    !IF (TMP .GE. 0.0D0) THEN
    !    DN1 = DSQRT(TMP)
    !ELSE
    !    DN1 = IM1 * DSQRT(DABS(TMP))
    !ENDIF
    !
    !DO I=1,NK
    !    TMP = 1.0D0 - (I * I) / 4.0D0 * (LAMBDA / LX) * (LAMBDA / LX)
    !    IF (TMP .GE. 0.0D0) THEN
    !        DN_RES(I) = DSQRT(TMP) - DN1
    !    ELSE
    !        DN_RES(I) = IM1 * DSQRT(DABS(TMP)) - DN1
    !    ENDIF
    !END DO

    OPEN (1, FILE='DELTA_N.DAT')
    DO I = 1, NK
        WRITE (1, '(I,2E17.8)') I, DN_RES(I)
    END DO
    CLOSE (1)
END FUNCTION DN_FN

SUBROUTINE ALLOCATE_ARRAYS()
    USE IC

    IMPLICIT NONE

    INTEGER(C_INT) ERR_ALLOC

    ALLOCATE (G(NX), A1(NX), A0(NX), AK1(NK), AK0(NK), JK1(NK), JK0(NK), ATMP(NX), &
              TH0(NTH, 2), TH1(NTH, 2), DTHDZ(NTH, 2), FK1(NK), FK2(NK), RHS0(NTH, 2), Z(NZ), X(NX), K2(NK), &
              A_AMP_Z0(NX), A_AMP_ZL(NX), AKTMP(NX), AKZL(NK), AKZ0(NK), A0Z0(NX), A0Z0CUT(NX), &
              A_SPEC_AMP_Z0(NX), A_SPEC_AMP_ZL(NX), IT(IT_TODO), &
              EX(NK), K(2*NK), DLT(NK), SUM_ABS2_A_PLUS_BY_Z(NX), SUM_ABS2_A_PLUS_BY_Z_K(NK), TMP(NX), &
              !THETA(NTH, 2, NZ), &
              STAT=ERR_ALLOC)

    IF (ERR_ALLOC /= 0) THEN
        PAUSE "ALLOCATION ERROR"
        STOP
    END IF
END SUBROUTINE ALLOCATE_ARRAYS

SUBROUTINE DEALLOCATE_ARRAYS()
    USE IC

    IMPLICIT NONE

    INTEGER(C_INT) ERR_DEALLOC

    DEALLOCATE (G, A1, A0, AK1, AK0, JK1, JK0, ATMP, &
                TH0, TH1, DTHDZ, FK1, FK2, RHS0, Z, X, K2, &
                A_AMP_Z0, A_AMP_ZL, AKTMP, AKZL, AKZ0, A0Z0CUT, &
                A_SPEC_AMP_Z0, A_SPEC_AMP_ZL, IT, &
                EX, K, DLT, SUM_ABS2_A_PLUS_BY_Z, SUM_ABS2_A_PLUS_BY_Z_K, TMP, &
                !THETA, &
                STAT=ERR_DEALLOC)

    IF (ERR_DEALLOC /= 0) STOP "DEALLOCATION ERROR"
END SUBROUTINE DEALLOCATE_ARRAYS

SUBROUTINE READ_PARAM() BIND(C, NAME='READ_PARAM')
    USE IC

    IMPLICIT NONE

    OPEN (UNIT=1, FILE='input_fortran.dat', STATUS='OLD', ERR=101)
    READ (UNIT=1, NML=PARAM, ERR=102)
    CLOSE (UNIT=1)

    WRITE (*, NML=PARAM)

    RETURN
101 PRINT *, 'ERROR OF FILE OPEN'; PAUSE; STOP
102 PRINT *, 'ERROR OF READING FILE "INPUT_FORTRAN.DAT"'; PAUSE; STOP
END SUBROUTINE READ_PARAM

!SUBROUTINE D(DADX)
!    USE IC, ONLY : K, NX, LX, PI, HX
!    USE FOURIER
!    USE, INTRINSIC :: ISO_C_BINDING
!
!    IMPLICIT NONE
!
!    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX) :: DADX
!    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(2*NX) :: V
!
!    V(1) = (0.0D0,0.0D0)
!    V(2:NX) = DADX(2:NX)
!    V(NX+1) = (0.0D0,0.0D0)
!    V(2*NX:NX+2:-1) = -DADX(2:NX)
!
!    CALL FFT(V)
!    V = 2.0 * PI / (2.0D0 * LX) * K * V
!    CALL IFFT(V)
!    DADX = V(1:NX)
!
!END SUBROUTINE D

!FUNCTION CHECK_ALLOC()
!USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_BOOL
!    USE IC, ONLY : A1, A0, AK1, AK0, ATMP, JK1, JK0, DADX, K, EX, &
!        TH0, TH1, DTHDZ, FK1, FK2, RHS0, Z, X, K2, &
!        A_AMP_Z0, A_AMP_ZL, A_SPEC_AMP_Z0, A_SPEC_AMP_ZL, G, SUM_ABS2_DADX_BY_X
!
!    IMPLICIT NONE
!
!    LOGICAL(C_BOOL) CHECK_ALLOC
!
!    !PRINT *, ALLOCATED(A0), SIZE(A0)
!
!    IF (.NOT. ALLOCATED(A0)) THEN
!        CHECK_ALLOC = .FALSE.
!        PRINT *, 'ARRAYS WERE NOT ALLOCATED!'
!        PAUSE
!        STOP
!    ENDIF
!
!    CHECK_ALLOC = .TRUE.
!END FUNCTION CHECK_ALLOC

SUBROUTINE MAKEA(ATMP, AK)
    USE IC, ONLY: NK, NX
    USE FOURIER
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT

    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: ATMP
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: AK
    INTEGER(C_INT) N1, N2

    N1 = SIZE(AK)
    N2 = SIZE(ATMP)

    IF (N1 .NE. NK .OR. N2 .NE. NX) THEN
        PRINT *, 'ERROR IN "MAKEA"'
        PAUSE
        STOP
    END IF

    ATMP = DCMPLX(0.0D0)
    ATMP(1:NK) = AK

    CALL ISINT(ATMP)
END SUBROUTINE MAKEA

SUBROUTINE MAKEAK(AK, ATMP)
    USE IC, ONLY: NK, NX, AKTMP
    USE FOURIER
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT

    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: AK
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: ATMP
    INTEGER(C_INT) N1, N2

    N1 = SIZE(AK)
    N2 = SIZE(ATMP)

    IF (N1 .NE. NK .OR. N2 .NE. NX) THEN
        PRINT *, 'ERROR IN "MAKEAK"'
        PAUSE
        STOP
    END IF

    AKTMP = ATMP

    CALL SINT(AKTMP)
    AK = AKTMP(1:NK)
END SUBROUTINE MAKEAK

SUBROUTINE MAKE_A0Z0_AK1_ATMP(A0Z0, AZ0CUT, AK0, AK)
    USE IC, ONLY: NK, NX, R0, G
    USE FOURIER
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT

    IMPLICIT NONE

    INTERFACE
        SUBROUTINE MAKEA(ATMP, AK)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: ATMP
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: AK
        END SUBROUTINE MAKEA
        SUBROUTINE MAKEAK(AK, ATMP)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE_COMPLEX, C_INT
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: AK
            COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: ATMP
        END SUBROUTINE MAKEAK
    END INTERFACE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(INOUT) :: A0Z0, AK0, AZ0CUT
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), INTENT(IN) :: AK
    INTEGER(C_INT) N1, N2, N3, N4

    N1 = SIZE(AK)
    N2 = SIZE(AZ0CUT)
    N3 = SIZE(A0Z0)
    N4 = SIZE(AK0)

    IF (N1 .NE. NK .OR. N2 .NE. NX .OR. N3 .NE. NX .OR. N4 .NE. NK) THEN
        PRINT *, 'ERROR IN "MAKEA"'
        PAUSE
        STOP
    END IF

    CALL MAKEA(A0Z0, AK) !ПЕРЕД ОБРЕЗАНИЕМ

    AZ0CUT = A0Z0*R0*G !ОТРАЖЕНИЕ ОТ ЗЕРКАЛА В Z=0 И ОБРЕЗАНИЕ

    CALL SINT(AZ0CUT)

    AK0 = AZ0CUT(1:NK)

    CALL MAKEA(AZ0CUT, AK0)
END SUBROUTINE MAKE_A0Z0_AK1_ATMP
