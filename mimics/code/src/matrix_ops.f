        SUBROUTINE EXP_MAT_test(LAMBDA,Q,QINV,z,ANSWER)
        save
c
        INTEGER I,J
        REAL Z
        COMPLEX LAMBDA(4),Q(4,4),QINV(4,4),D(4,4),ANSWER(4,4)
c
        DO I=1,4
         DO J=1,4
            D(J,I) = CMPLX(0.0,0.0)
	 END DO
	END DO
C
        DO I=1,4
            D(I,I) = CEXP(LAMBDA(I)*Z)
	END DO
C
        CALL CMATMULT3(Q,D,QINV,ANSWER)
C
        RETURN
        END
c
        SUBROUTINE EXP_MAT(LAMBDA,Q,QINV,z,ANSWER)
        save
c
        INTEGER I,J,K
        REAL Z ,ANSWER(4,4)
        REAL Qr(4,4),QINVr(4,4),Lr(4,4),Qi(4,4),QINVi(4,4),Li(4,4)
        REAL FACTRR(4,4), FACTRI(4,4), FACTII(4,4), FACTIR(4,4)
        REAL DIFF(4,4), SUM(4,4)
        COMPLEX LAMBDA(4),Q(4,4),QINV(4,4),D(4,4)
c
        DO I=1,4
         DO J=1,4
            D(J,I) = CMPLX(0.0,0.0)
            Lr(J,I) = 0.0
            Li(J,I) = 0.0
            Qr(J,I) = REAL(Q(J,I))
            Qi(J,I) = AIMAG(Q(J,I))
            QINVr(J,I) = REAL(QINV(J,I))
            QINVi(J,I) = AIMAG(QINV(J,I))
            ANSWER(J,I) = 0.0
	 END DO
	END DO
C
        DO I=1,4
            D(I,I) = CEXP(LAMBDA(I)*Z)
            Lr(I,I) = REAL(D(I,I))
            Li(I,I) = AIMAG(D(I,I))
	END DO
C
        DO I=1,4
         DO J=1,4
            FACTRR(J,I) = QINVR(J,I)*LR(J,J) 
            FACTRI(J,I) = QINVR(J,I)*LI(J,J) 
            FACTII(J,I) = QINVI(J,I)*LI(J,J) 
            FACTIR(J,I) = QINVI(J,I)*LR(J,J) 
	 END DO
	END DO
C
        DO I=1,4
         DO J=1,4
            DIFF(J,I) = FACTRR(J,I) - FACTII(J,I)
            SUM(J,I)  = FACTRI(J,I) + FACTIR(J,I)
	 END DO
	END DO
C
        DO K=1,4
        DO J=1,4
        DO I=1,4
         ANSWER(K,J) = ANSWER(K,J) + QR(K,I)*DIFF(I,J)-QI(K,I)*SUM(I,J)
	END DO
	END DO
	END DO
C
        RETURN
        END
c
        SUBROUTINE INVERSE_MAT(INPUT,OUTPUT)
        save
c
        INTEGER IV(4),JOB,I,J
        REAL RC
        COMPLEX SV(4),DET(2)
        COMPLEX INPUT(4,4),OUTPUT(4,4)
c

        DO I=1,4
         DO J=1,4
            OUTPUT(J,I)=INPUT(J,I)
	 END DO
	END DO 
c
        CALL CGECO(OUTPUT,4,4,IV,RC,SV)   
        JOB = 1
        CALL CGEDI(OUTPUT,4,4,IV,DET,SV,JOB) 
C
        RETURN
        END
c
        SUBROUTINE CMATMULT2(MAT1,MAT2,OUTPUT)
        save
c
        INTEGER I,J,K
        COMPLEX MAT1(4,4),MAT2(4,4),OUTPUT(4,4)
c
        DO I=1,4
         DO J=1,4
            OUTPUT(J,I) = CMPLX(0.0,0.0)
	 END DO
	END DO
C
        DO K=1,4
          DO J=1,4
            DO I=1,4
              OUTPUT(J,K) = OUTPUT(J,K) + MAT1(J,I)*MAT2(I,K)
	  END DO
	 END DO
	END DO
C
        RETURN
        END
c
        SUBROUTINE CMATMULT3(MAT1,MAT2,MAT3,OUTPUT)
        save
c
        COMPLEX MAT1(4,4),MAT2(4,4),MAT3(4,4),MAT4(4,4),OUTPUT(4,4)
c
        CALL CMATMULT2(MAT1,MAT2,MAT4)
        CALL CMATMULT2(MAT4,MAT3,OUTPUT)
C
        RETURN
        END
c
        SUBROUTINE CrMATMULT3(MAT1,MAT2,MAT3,OUTPUT)
        save
c
        REAL MAT2(4,4)
        COMPLEX MAT1(4,4),MAT3(4,4),MAT4(4,4),OUTPUT(4,4)
c
        CALL CrMATMULT2(MAT1,MAT2,MAT4)
        CALL CMATMULT2(MAT4,MAT3,OUTPUT)
C
        RETURN
        END
c
        SUBROUTINE CrMATMULT2(MAT1,MAT2,OUTPUT)
        save
c
        INTEGER I,J,K
        REAL MAT2(4,4)
        COMPLEX MAT1(4,4),OUTPUT(4,4)
c
        DO I=1,4
         DO J=1,4
            OUTPUT(J,I) = CMPLX(0.0,0.0)
	 END DO
	END DO
C
        DO K=1,4
          DO J=1,4
            DO I=1,4
              OUTPUT(J,K) = OUTPUT(J,K) + MAT1(J,I)*MAT2(I,K)
	  END DO
	 END DO
	END DO
C
        RETURN
        END
c
        SUBROUTINE MATMULT3(MAT1,MAT2,MAT3,OUTPUT)
        save
c
        REAL MAT1(4,4),MAT2(4,4),MAT3(4,4),MAT4(4,4),OUTPUT(4,4)
c
        CALL MATMULT2(MAT1,MAT2,MAT4)
        CALL MATMULT2(MAT4,MAT3,OUTPUT)
C
        RETURN
        END
c
        SUBROUTINE MATMULT2(MAT1,MAT2,OUTPUT)
        save
c
        INTEGER I,J,K
        REAL MAT1(4,4),MAT2(4,4),OUTPUT(4,4)
c
        DO I=1,4
         DO J=1,4
            OUTPUT(J,I) = 0.0
	 END DO
	END DO
C
        DO K=1,4
          DO J=1,4
            DO I=1,4
              OUTPUT(J,K) = OUTPUT(J,K) + MAT1(J,I)*MAT2(I,K)
	  END DO
	 END DO
	END DO
C
        RETURN
        END
c
        SUBROUTINE MATMULT3Cr(Q,A,QINV,ANSWER)
        save
c
        INTEGER I,J,K
        REAL ANSWER(4,4)
        REAL Qr(4,4),QINVr(4,4),Lr(4,4),Qi(4,4),QINVi(4,4),Li(4,4)
        REAL FACTRR(4,4), FACTRI(4,4), FACTII(4,4), FACTIR(4,4)
        REAL DIFF(4,4), SUM(4,4)
        COMPLEX Q(4,4),QINV(4,4),A(4,4)
c
        DO I=1,4
         DO J=1,4
            Lr(J,I) = REAL(A(J,I))
            Li(J,I) = AIMAG(A(J,I))
            Qr(J,I) = REAL(Q(J,I))
            Qi(J,I) = AIMAG(Q(J,I))
            QINVr(J,I) = REAL(QINV(J,I))
            QINVi(J,I) = AIMAG(QINV(J,I))
            ANSWER(J,I) = 0.0
            FACTRR(I,J) = 0.0
            FACTII(I,J) = 0.0
            FACTRI(I,J) = 0.0
            FACTIR(I,J) = 0.0
	 END DO
	END DO
c
        DO J = 1, 4
         DO I = 1, 4
          DO K = 1, 4
            FACTRR(I,J) = FACTRR(I,J) + QINVR(K,J)*LR(I,K)
            FACTII(I,J) = FACTII(I,J) + QINVI(K,J)*LI(I,K)
            FACTRI(I,J) = FACTRI(I,J) + QINVR(K,J)*LI(I,K)
            FACTIR(I,J) = FACTIR(I,J) + QINVI(K,J)*LR(I,K)
          END DO
         END DO
        END DO
C
        DO J = 1, 4
         DO I = 1, 4
          DIFF(I,J) = FACTRR(I,J) - FACTII(I,J)
          SUM(I,J)  = FACTRI(I,J) + FACTIR(I,J)
         END DO
        END DO
C
        DO J = 1,4
         DO I = 1, 4
          DO K = 1, 4
           ANSWER(I,J) = ANSWER(I,J) +
     &                   QR(I,K)*DIFF(K,J) - QI(I,K)*SUM(K,J)
          END DO
         END DO
        END DO
C
        RETURN
        END
c
        SUBROUTINE EXTINCT_SUB(M,KAPPA)
        save
c
        REAL KAPPA(4,4)
c
        COMPLEX M(2,2)
c
        KAPPA(1,1) = -2.0*REAL(M(1,1))
        KAPPA(2,1) =  0.0
        KAPPA(3,1) = -2.0*REAL(M(2,1))
        KAPPA(4,1) =  2.0*AIMAG(M(2,1))
        KAPPA(1,2) =  0.0
        KAPPA(2,2) = -2.0*REAL(M(2,2))
        KAPPA(3,2) = -2.0*REAL(M(1,2))
        KAPPA(4,2) = -2.0*AIMAG(M(1,2))
        KAPPA(1,3) = -REAL(M(1,2))
        KAPPA(2,3) = -REAL(M(2,1))
        KAPPA(3,3) = -(REAL(M(1,1)) + REAL(M(2,2))) 
        KAPPA(4,3) = -(AIMAG(M(1,1)) - AIMAG(M(2,2))) 
        KAPPA(1,4) = -2.0*REAL(M(1,2))
        KAPPA(2,4) = -2.0*REAL(M(2,1))
        KAPPA(3,4) =  (AIMAG(M(1,1)) - AIMAG(M(2,2))) 
        KAPPA(4,4) = -(REAL(M(1,1)) + REAL(M(2,2))) 
C
        RETURN
        END
c
        SUBROUTINE MAT_EIGEN_SUB(M,LAMBDA,Q,QINV)
        save
c
        REAL TOL, WORK1, WORK2, CHECK1, CHECK2
        REAL FREQ_HERTZ, WAVELENGTH, k0
c
        COMPLEX M(2,2),LAMBDA(4),Q(4,4),QINV(4,4),KV1(2),KV2(2),K1,K2
        COMPLEX r,b1,b2,CDIFF,CSUM,CWORK1,CWORK2
        COMPLEX J, MYCSQRT
C
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
c
        J = CMPLX(0.0,1.0)
        TOL = 1.0E-3
c
        IF(CABS(M(1,2)).EQ.0.0)THEN
            CHECK1 = 0.0
        ELSE IF(CABS(M(1,1)).NE.0.0)THEN
            CHECK1 = CABS(M(1,2)/M(1,1))
        ELSE
            CHECK1 = TOL + 1.0
        ENDIF
C
        IF(CABS(M(2,1)).EQ.0.0)THEN
            CHECK2 = 0.0
        ELSE IF(CABS(M(2,2)).NE.0.0)THEN
            CHECK2 = CABS(M(2,1)/M(2,2))
        ELSE
            CHECK2 = TOL + 1.0
        ENDIF
C
        IF((CHECK1.LT.TOL).AND.(CHECK2.LT.TOL)) THEN
c
            K1 = k0 - J*M(1,1)
            K2 = k0 - J*M(2,2)
C
            KV1(1) = CMPLX(1.0,0.0)
            KV1(2) = CMPLX(0.0,0.0)
            KV2(1) = CMPLX(0.0,0.0)
            KV2(2) = CMPLX(1.0,0.0)
C
            WORK1 = -REAL(M(1,1) + M(2,2))
            WORK2 = AIMAG(M(1,1) - M(2,2))
C
            LAMBDA(1) = -2.0*REAL(M(1,1))
            LAMBDA(2) =  CMPLX(WORK1,-WORK2)
            LAMBDA(3) =  CMPLX(WORK1,WORK2)
            LAMBDA(4) = -2.0*REAL(M(2,2))
C
            Q(1,1) = CMPLX(1.0,0.0)
            Q(2,1) = CMPLX(0.0,0.0)
            Q(3,1) = CMPLX(0.0,0.0)
            Q(4,1) = CMPLX(0.0,0.0)
            Q(1,2) = CMPLX(0.0,0.0)
            Q(2,2) = CMPLX(0.0,0.0)
            Q(3,2) = CMPLX(1.0,0.0)
            Q(4,2) = CMPLX(0.0,-1.0)
            Q(1,3) = CMPLX(0.0,0.0)
            Q(2,3) = CMPLX(0.0,0.0)
            Q(3,3) = CMPLX(1.0,0.0)
            Q(4,3) = CMPLX(0.0,1.0)
            Q(1,4) = CMPLX(0.0,0.0)
            Q(2,4) = CMPLX(1.0,0.0)
            Q(3,4) = CMPLX(0.0,0.0)
            Q(4,4) = CMPLX(0.0,0.0)
c
            QINV(1,1) = CMPLX(1.0,0.0)
            QINV(2,1) = CMPLX(0.0,0.0)
            QINV(3,1) = CMPLX(0.0,0.0)
            QINV(4,1) = CMPLX(0.0,0.0)
            QINV(1,2) = CMPLX(0.0,0.0)
            QINV(2,2) = CMPLX(0.0,0.0)
            QINV(3,2) = CMPLX(0.0,0.0)
            QINV(4,2) = CMPLX(1.0,0.0)
            QINV(1,3) = CMPLX(0.0,0.0)
            QINV(2,3) = CMPLX(0.5,0.0)
            QINV(3,3) = CMPLX(0.5,0.0)
            QINV(4,3) = CMPLX(0.0,0.0)
            QINV(1,4) = CMPLX(0.0,0.0)
            QINV(2,4) = CMPLX(0.0,0.5)
            QINV(3,4) = CMPLX(0.0,-0.5)
            QINV(4,4) = CMPLX(0.0,0.0)
c
        else
c
            CDIFF = M(1,1)-M(2,2)
            CSUM =  M(1,1)+M(2,2)
c
            r = MYCSQRT(CDIFF*CDIFF +4.0*M(2,1)*M(1,2))
            K1 = k0 - (J/2.)*(CSUM + r)
            K2 = k0 - (J/2.)*(CSUM - r)
C
            b1 = 2.0*M(2,1)/(CDIFF+r)
            b2 = 2.0*M(1,2)/(-CDIFF-r)
C
            K1 = k0 - J*(M(1,1)+M(2,2)+r)/2.0
            K2 = k0 - J*(M(1,1)+M(2,2)-r)/2.0
C
            KV1(1) = CMPLX(1.0,0.0)
c            KV1(2) = CMPLX(0.0,b1)
            KV1(2) = b1
c            KV2(1) = CMPLX(b2,0.0)
            KV2(1) = b2
            KV2(2) = CMPLX(1.0,0.0)
C
            LAMBDA(1) = 2.0*AIMAG(K1)
            LAMBDA(2) = J*(CONJG(K2) - K1)
            LAMBDA(3) = J*(CONJG(K1) - K2)
            LAMBDA(4) = 2.0*AIMAG(K2)
C
            WORK1 = CABS(b1)
            WORK2 = CABS(b2)
            CWORK1 = CONJG(b1)
            CWORK2 = CONJG(b2)
C
            Q(1,1) =  CMPLX(1.0,0.0)
            Q(2,1) =  WORK1*WORK2
            Q(3,1) =  2.0*REAL(b1)
            Q(4,1) = -2.0*AIMAG(b1)
            Q(1,2) =  CWORK2
            Q(2,2) =  b1
            Q(3,2) =  1.0 + b1*CWORK2
            Q(4,2) = -J*(1.0 - b1*CWORK2)
            Q(1,3) =  b2
            Q(2,3) =  CWORK1
            Q(3,3) =  1.0 + b2*CWORK1
            Q(4,3) =  J*(1.0 - b2*CWORK1)
            Q(1,4) =  WORK2*WORK2
            Q(2,4) =  CMPLX(1.0,0.0)
            Q(3,4) =  2.0*REAL(b2)
            Q(4,4) =  2.0*AIMAG(b2)
C
            CALL INVERSE_MAT(Q,QINV)
        ENDIF
C
        RETURN
        END
c
        subroutine modify_mat(mat) 
        save
c
        integer i, j
        real mat(4,4), u(4,4), utranspose(4,4), mat_mod(4,4)
c
        u(1,1) =  1.0
        u(2,1) =  1.0
        u(3,1) =  0.0
        u(4,1) =  0.0
        u(1,2) =  1.0
        u(2,2) = -1.0
        u(3,2) =  0.0
        u(4,2) =  0.0
        u(1,3) =  0.0
        u(2,3) =  0.0
        u(3,3) =  1.0
        u(4,3) =  0.0
        u(1,4) =  0.0
        u(2,4) =  0.0
        u(3,4) =  0.0
        u(4,4) =  1.0
c
        do i=1,4
          do j=1,4
            utranspose(j,i) = u(i,j)
          enddo
        enddo
c
        call MATMULT3(utranspose,mat,u,mat_mod)
c
        do i=1,4
          do j=1,4
            mat(j,i) = mat_mod(j,i)
          enddo
        enddo
c
        return
        end
