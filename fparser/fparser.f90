! Copyright (c) 2000-2008, Roland Schmehl.
!
! All rights reserved.
!
! * Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! * Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! * Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! * Neither the name of the copyright holder nor the names of its
! contributors may be used to endorse or promote products derived from
! this software without specific prior written permission.
!
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

MODULE fparser_mod
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! Fortran 90 function parser v1.1
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  !
  ! This function parser module is intended for applications where a set of mathematical
  ! fortran-style expressions is specified at runtime and is then evaluated for a large 
  ! number of variable values. This is done by compiling the set of function strings 
  ! into byte code, which is interpreted efficiently for the various variable values. 
  !
  ! The source code is available from http://fparser.sourceforge.net
  !
  ! Please send comments, corrections or questions to the author:
  ! Roland Schmehl <roland.schmehl@alumni.uni-karlsruhe.de>
  !
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! The function parser concept is based on a C++ class library written by  Juha 
  ! Nieminen <warp@iki.fi> available from http://warp.povusers.org/FunctionParser/
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  USE prec_mod
  USE constants_mod
  USE fparser_dat_mod
  IMPLICIT NONE
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PUBLIC                     :: initf,    &     ! Initialize function parser for n functions
                                parsef,   &     ! Parse single function string
                                evalf,    &     ! Evaluate single function
                                evald,    &     ! Evaluate first derivative of single function
                                evaldd,   &     ! Evaluate second derivative of single function
                                EvalfErrMsg, &  ! Error message (Use only when EvalfErrType\='')
                                EvaldErrMsg     ! Error message (Use only when EvaldErrType\='')
  PUBLIC                     :: RemoveAllSpaces,      &
                                CheckVariables,       &
                                CheckVarExpressions,  &
                                InsertVarExpressions
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PRIVATE

  INTEGER :: EvalfErrType ! = 0: no error occured, else: evaluation error
  INTEGER :: EvaldErrType ! = 0: no error occured, else: evaluation error

  INTERFACE evalf
     MODULE PROCEDURE evalf1, evalfn
  END INTERFACE

  INTERFACE evald
     MODULE PROCEDURE evald1, evaldm, evaldn, evaldmn
  END INTERFACE

  INTERFACE evaldd
     MODULE PROCEDURE evaldd1
  END INTERFACE

  SAVE
  !
CONTAINS
  !
  SUBROUTINE help_operators
    IMPLICIT NONE

    WRITE(*,*)
    WRITE(*,*) 'fparser - available operators: '
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '+ (unary plus): -arg1          '
    WRITE(*,*) '- (unary minus): -arg1         '
    WRITE(*,*) '+ (addition): arg1 + arg2      '
    WRITE(*,*) '- (subtraction): arg1 - arg2   '
    WRITE(*,*) '* (multiplication): arg1 * arg2'
    WRITE(*,*) '/ (division): arg1 / arg2      '
    WRITE(*,*) '^ (power): arg1 ^ arg2         '
    WRITE(*,*) '** (power): arg1 ** arg2       '
    WRITE(*,*)

  END SUBROUTINE help_operators
  !
  SUBROUTINE help_functions
    IMPLICIT NONE

    WRITE(*,*)
    WRITE(*,*) 'fparser - available functions:                                         '
    WRITE(*,*) '-----------------------------------------------------------------------'
    WRITE(*,*) 'abs(arg1):       absolute value                                        '
    WRITE(*,*) 'exp(arg1):       exponential                                           '
    WRITE(*,*) 'log10(arg1):     common logarithm                                      '        
    WRITE(*,*) 'log(arg1):       natural logarithm                                     '
    WRITE(*,*) 'sqrt(arg1):      square root                                           '
    WRITE(*,*) 'sinh(arg1):      hyperbolic sine                                       '
    WRITE(*,*) 'cosh(arg1):      hyperbolic cosine                                     '
    WRITE(*,*) 'tanh(arg1):      hyperbolic tangent                                    '
    WRITE(*,*) 'sin(arg1):       sine                                                  '
    WRITE(*,*) 'cos(arg1):       cosine                                                '
    WRITE(*,*) 'tan(arg1):       tangent                                               '
    WRITE(*,*) 'asinh(arg1):     inverse of hyperbolic sine                            '
    WRITE(*,*) 'acosh(arg1):     inverse of hyperbolic cosine                          '
    WRITE(*,*) 'atanh(arg1):     inverse of hyperbolic tangent                         '
    WRITE(*,*) 'asin(arg1):      inverse of sine                                       '
    WRITE(*,*) 'acos(arg1):      inverse of cosine                                     '
    WRITE(*,*) 'atan(arg1):      inverse of tangent                                    '
    WRITE(*,*) 'erfc(arg1):      inverse of error function                             '
    WRITE(*,*) 'erf(arg1):       error function                                        '
    WRITE(*,*) 'sgn(arg1):       signum function                                       '
    WRITE(*,*) 'sign(arg1,arg2): signum function with two arguments                    '
    WRITE(*,*) 'logb(arg1,arg2): general logarithm of arg1 with base arg2              '
    WRITE(*,*) 'sum(arg1,...,argN): sum of N>=1 arguments                              '
    WRITE(*,*) 'ave(arg1,...,argN): average of N>=1 arguments (=mean, =amean)          '
    WRITE(*,*) 'mean(arg1,...,argN): arithemtic mean of N>=1 arguments (=ave, =amean)  '
    WRITE(*,*) 'amean(arg1,...,argN): arithemtic mean of N>=1 arguments (=ave, =mean)  '
    WRITE(*,*) 'prod(arg1,...,argN): product of N>=1 arguments                         '
    WRITE(*,*) 'gmean(arg1,...,argN): geometric mean of N>=1 arguments                 '
    WRITE(*,*) 'hmean(arg1,...,argN): harmonic mean of N>=1 arguments                  '
    WRITE(*,*) 'min(arg1,...,argN): minimum of N>=1 arguments                          '
    WRITE(*,*) 'max(arg1,...,argN): maximum of N>=1 arguments                          '
    WRITE(*,*) 'dot(arg1,...,arg6): dot product of two vectors                         '
    WRITE(*,*) 'crossx(arg1,...,arg6): first component of cross product of two vectors '
    WRITE(*,*) 'crossy(arg1,...,arg6): second component of cross product of two vectors'
    WRITE(*,*) 'crossz(arg1,...,arg6): third component of cross product of two vectors '
    WRITE(*,*) 'norm(arg1,...,arg3): norm of a vector                                  '
    WRITE(*,*) 'dis(arg1,...,arg6): distance between two points                        '
    WRITE(*,*) 'ang(arg1,...,arg9): angle betweeb three points                         '
    WRITE(*,*) 'dih(arg1,...,arg12): dihedral angle between four points                '
    WRITE(*,*)

  END SUBROUTINE help_functions
  !
  SUBROUTINE initf (Comp, n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Initialize function parser for n functions
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(INOUT) :: Comp
    INTEGER, INTENT(in) :: n                                 ! Number of functions
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    Comp => NULL()
    ALLOCATE (Comp(n))
  END SUBROUTINE initf
  !
  SUBROUTINE parsef (Comp, i, FuncStr, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Parse ith function string FuncStr and compile it into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(INOUT) :: Comp
    INTEGER,                         INTENT(in) :: i         ! Function identifier
    CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    CHARACTER (LEN=LEN(FuncStr))                :: Func      ! Function string, local use
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (i < 1 .OR. i > SIZE(Comp)) THEN
       WRITE(*,*) '*** Parser error: Function number ',i,' out of range'
       STOP
    END IF
    ALLOCATE (ipos(LEN_TRIM(FuncStr)))                       ! Char. positions in orig. string
    Func = FuncStr                                           ! Local copy of function string
    CALL Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
    CALL RemoveSpaces (Func)                                 ! Condense function string
    CALL CheckSyntax (Func,FuncStr,Var)
    DEALLOCATE (ipos)
    CALL Compile (Comp,i,Func,Var)                           ! Compile into bytecode
  END SUBROUTINE parsef
  !
  FUNCTION evalf1(Comp, i, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function for the values passed in array Val(:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                INTENT(in) :: i                  ! Function identifier
    REAL(dp), DIMENSION(:), INTENT(in) :: Val                ! Variable values
    REAL(dp)                           :: res                ! Result
    INTEGER                            :: InP,             & ! Instruction pointer
                                          DaP,             & ! Data pointer
                                          StP                ! Stack pointer
    REAL(dp)                           :: Stack(Comp(i)%StackSize)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    EvalfErrType = 0
    DaP = 1
    StP = 0
    DO InP=1,Comp(i)%ByteCodeSize
       SELECT CASE (Comp(i)%ByteCode(InP))

       CASE(cImmed)
            StP=StP+1
            Stack(StP)=Comp(i)%Immed(DaP)
            DaP=DaP+1
       CASE(cNeg)
            Stack(StP)=-Stack(StP)
       CASE(cAdd)
            Stack(StP-1)=Stack(StP-1)+Stack(StP)
            StP=StP-1
       CASE(cSub)
            Stack(StP-1)=Stack(StP-1)-Stack(StP)
            StP=StP-1
       CASE(cMul)
            Stack(StP-1)=Stack(StP-1)*Stack(StP)
            StP=StP-1
       CASE(cDiv)
            IF(Stack(StP)==zero) THEN
               EvalfErrType=1
               res=zero
               RETURN
            ENDIF
            Stack(StP-1)=Stack(StP-1)/Stack(StP)
            StP=StP-1
       CASE(cPow)
            IF(Stack(StP-1)==zero) THEN
               IF(Stack(StP)<=zero) THEN
                  EvalfErrType=1
                  res=zero
                  RETURN
               ELSE
                  Stack(StP-1)=zero
               END IF
            ELSE IF(Stack(StP-1)<zero) THEN
               Stack(StP-1)=Stack(StP-1)**NINT(Stack(StP)) ! lam81
            ELSE
               Stack(StP-1)=Stack(StP-1)**Stack(StP)
            ENDIF
            StP=StP-1
       CASE(cAbs)
            Stack(StP)=ABS(Stack(StP))
       CASE(cExp)
            Stack(StP)=EXP(Stack(StP))
       CASE(cLog10)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=2
               res=zero
               RETURN
            ENDIF
            Stack(StP)=LOG10(Stack(StP))
       CASE(cLog)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=3
               res=zero
               RETURN
            ENDIF
            Stack(StP)=LOG(Stack(StP))
       CASE(cSqrt)
            IF(Stack(StP)<zero) THEN
               EvalfErrType=4
               res=zero
               RETURN
            ENDIF
            Stack(StP)=SQRT(Stack(StP))
       CASE(cSinh)
            Stack(StP)=SINH(Stack(StP))
       CASE(cCosh)
            Stack(StP)=COSH(Stack(StP))
       CASE(cTanh)
            Stack(StP)=TANH(Stack(StP))
       CASE(cSin)
            Stack(StP)=SIN(Stack(StP))
       CASE(cCos)
            Stack(StP)=COS(Stack(StP))
       CASE(cTan)
            Stack(StP)=TAN(Stack(StP))
       CASE(cAsinh)
            Stack(StP)=ASINH(Stack(StP))
       CASE(cAcosh)
            IF(Stack(StP)<one) THEN
               EvalfErrType=5
               res=zero
               RETURN
            ENDIF
            Stack(StP)=ACOSH(Stack(StP))
       CASE(cAtanh)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=6
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=7
               res=zero
               RETURN
            ENDIF
            Stack(StP)=ATANH(Stack(StP))
       CASE(cAsin)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=8
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=9
               res=zero
               RETURN
            ENDIF
            Stack(StP)=ASIN(Stack(StP))
       CASE(cAcos)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=10
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=11
               res=zero
               RETURN
            ENDIF
            Stack(StP)=ACOS(Stack(StP))
       CASE(cAtan)
            Stack(StP)=ATAN(Stack(StP))
       CASE(cErfc)
            Stack(StP)=ERFC(Stack(StP))
       CASE(cErf)
            Stack(StP)=ERF(Stack(StP))
       CASE(cSgn)
            Stack(StP)=SIGN(one,Stack(StP))
       CASE(cSign)
            Stack(StP-1)=SIGN(Stack(StP-1),Stack(StP))
       CASE(cLogb)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=12
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvalfErrType=13
               res=zero
               RETURN
            END IF
            IF(Stack(StP-1)<=zero) THEN
               EvalfErrType=14
               res=zero
               RETURN
            END IF
            Stack(StP-1)=LOG(Stack(StP-1))/LOG(Stack(StP))
            StP=StP-1
       CASE(cSum)
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cAve, cMean, cAmean)
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(Stack(StP-Comp(i)%NArg(InP)+1:StP))/REAL(Comp(i)%NArg(InP),dp)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cProd)
            Stack(StP-Comp(i)%NArg(InP)+1)=PRODUCT(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cGmean)
            Stack(StP-Comp(i)%NArg(InP)+1)=PRODUCT(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            IF(Stack(StP-Comp(i)%NArg(InP)+1)<zero) THEN
                EvalfErrType=15
                res=zero
                RETURN
            ENDIF
            Stack(StP-Comp(i)%NArg(InP)+1)=Stack(StP-Comp(i)%NArg(InP)+1)**(one/REAL(Comp(i)%NArg(InP),dp))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cHmean)
            IF( ANY(Stack(StP-Comp(i)%NArg(InP)+1:StP)==zero) ) THEN
                EvalfErrType=16
                res=zero
                RETURN
            END IF
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(one/Stack(StP-Comp(i)%NArg(InP)+1:StP))
            IF( Stack(StP-Comp(i)%NArg(InP)+1)==zero ) THEN
                EvalfErrType=17
                res=zero
                RETURN
            END IF
            Stack(StP-Comp(i)%NArg(InP)+1)=REAL(Comp(i)%NArg(InP),dp)/Stack(StP-Comp(i)%NArg(InP)+1)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cMin)
            Stack(StP-Comp(i)%NArg(InP)+1)=MINVAL(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cMax)
            Stack(StP-Comp(i)%NArg(InP)+1)=MAXVAL(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cWrap)
            Stack(StP-1)=Stack(StP-1)-NINT(Stack(StP-1)/Stack(StP))*Stack(StP)
            StP=StP-1
       CASE(cDot)
            Stack(StP-5)=dotf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossx)
            Stack(StP-5)=crossxf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossy)
            Stack(StP-5)=crossyf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossz)
            Stack(StP-5)=crosszf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cNorm)
            Stack(StP-2)=normf(Stack(StP-2:StP))
            StP=StP-2
       CASE(cDis)
            Stack(StP-5)=disf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cAng)
            Stack(StP-8)=angf(Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-8
       CASE(cDih)
            Stack(StP-11)=dihf(Stack(StP-11:StP-9),Stack(StP-8:StP-6), &
                                       Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-11
       CASE DEFAULT
             StP=StP+1
             Stack(StP)=Val(Comp(i)%ByteCode(InP)-VarBegin+1)
       END SELECT
    END DO
    res = Stack(1)
  END FUNCTION evalf1
  !
  FUNCTION evalfn(Comp, i, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function for the values passed in array Val(:,:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                  INTENT(in) :: i                ! Function identifier
    REAL(dp), DIMENSION(:,:), INTENT(in) :: Val              ! Variable values
    REAL(dp), DIMENSION(SIZE(Val,2))     :: res              ! Result

    INTEGER                              :: n
    !----- -------- --------- --------- --------- --------- --------- --------- -------

    DO n=1, SIZE(Val,2)
       res(n) = evalf1(Comp, i, Val(:,n))
    END DO
  END FUNCTION evalfn
  !
  FUNCTION evald1(Comp, i, j, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function's jth derivative for the values passed in array Val(:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                INTENT(in) :: i                  ! Function identifier
    INTEGER,                INTENT(in) :: j                  ! Variable identifier
    REAL(dp), DIMENSION(:), INTENT(in) :: Val                ! Variable values
    REAL(dp)                           :: res                ! Result
    INTEGER                            :: InP,             & ! Instruction pointer
                                          DaP,             & ! Data pointer
                                          StP                ! Stack pointer
    INTEGER                            :: n
    INTEGER                            :: t
    REAL(dp)                           :: d
    REAL(dp)                           :: tmp1, tmp2
    REAL(dp)                           :: Stack(Comp(i)%StackSize)
    REAL(dp)                           :: DStack(1,Comp(i)%StackSize)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF(Comp(i)%Var(j) == 0) THEN
       res = zero
       RETURN
    END IF
    EvalfErrType = 0
    EvaldErrType = 0
    DaP = 1
    StP = 0
    DO InP=1,Comp(i)%ByteCodeSize
       SELECT CASE (Comp(i)%ByteCode(InP))

       CASE(cImmed)
            StP=StP+1
            DStack(1,StP)=zero
            Stack(StP)=Comp(i)%Immed(DaP)
            DaP=DaP+1
       CASE(cNeg)
            DStack(1,StP)=-DStack(1,StP)
            Stack(StP)=-Stack(StP)
       CASE(cAdd)
            DStack(1,StP-1)=DStack(1,StP-1)+DStack(1,StP)
            Stack(StP-1)=Stack(StP-1)+Stack(StP)
            StP=StP-1
       CASE(cSub)
            DStack(1,StP-1)=DStack(1,StP-1)-DStack(1,StP)
            Stack(StP-1)=Stack(StP-1)-Stack(StP)
            StP=StP-1
       CASE(cMul)
            DStack(1,StP-1)=DStack(1,StP-1)*Stack(StP)+Stack(StP-1)*DStack(1,StP)
            Stack(StP-1)=Stack(StP-1)*Stack(StP)
            StP=StP-1
       CASE(cDiv)
            IF(Stack(StP)==zero) THEN
               EvalfErrType=1
               res=zero
               RETURN
            ENDIF
            DStack(1,StP-1)=( DStack(1,StP-1)*Stack(StP)-Stack(StP-1)*DStack(1,StP) ) / &
                                 ( Stack(StP)**2 )
            Stack(StP-1)=Stack(StP-1)/Stack(StP)
            StP=StP-1
       CASE(cPow)
            IF(Stack(StP-1)==zero) THEN
               IF(Stack(StP)<=zero) THEN
                  EvalfErrType=1
                  res=zero
                  RETURN
               ELSE IF(Stack(StP)==one) THEN ! special case when exponent=1
                  DStack(1,StP-1)=DStack(1,StP-1)
                  Stack(StP-1)=zero
               ELSE
                  DStack(1,StP-1)=zero
                  Stack(StP-1)=zero
               ENDIF
            ELSE IF(Stack(StP-1)<zero) THEN
               tmp1=Stack(StP-1)**NINT(Stack(StP))
               DStack(1,StP-1)=tmp1 * ( DStack(1,StP)*LOG(-Stack(StP-1)) + &
                                       Stack(StP)/Stack(StP-1)*DStack(1,StP-1) )
               Stack(StP-1)=tmp1
            ELSE
               tmp1=Stack(StP-1)**Stack(StP)
               DStack(1,StP-1)=tmp1 * ( DStack(1,StP)*LOG(Stack(StP-1)) + &
                                       Stack(StP)/Stack(StP-1)*DStack(1,StP-1) )
               Stack(StP-1)=tmp1
            ENDIF
            StP=StP-1
       CASE(cAbs)
            IF(Stack(StP)==zero) THEN
               EvaldErrType=1
               res=zero
               RETURN
            ENDIF
            d=Stack(StP)/ABS(Stack(StP))
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ABS(Stack(StP))
       CASE(cExp)
            tmp1=EXP(Stack(StP))
            d=tmp1
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=tmp1
       CASE(cLog10)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=2
               res=zero
               RETURN
            ENDIF
            d=one/(LOG(ten)*Stack(StP))
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=LOG10(Stack(StP))
       CASE(cLog)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=3
               res=zero
               RETURN
            ENDIF
            d=one/Stack(StP)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=LOG(Stack(StP))
       CASE(cSqrt)
            IF(Stack(StP)<zero) THEN
               EvalfErrType=4
               res=zero
               RETURN
            ELSE IF(Stack(StP)==zero) THEN
               EvaldErrType=2
               res=zero
               RETURN
            ENDIF
            tmp1=SQRT(Stack(StP))
            d=one/(two*tmp1)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=tmp1
       CASE(cSinh)
            d=COSH(Stack(StP))
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=SINH(Stack(StP))
       CASE(cCosh)
            d=SINH(Stack(StP))
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=COSH(Stack(StP))
       CASE(cTanh)
            tmp1=TANH(Stack(StP))
            d=one-tmp1**2
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=tmp1
       CASE(cSin)
            d=COS(Stack(StP))
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=SIN(Stack(StP))
       CASE(cCos)
            d=-SIN(Stack(StP))
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=COS(Stack(StP))
       CASE(cTan)
            d=one/(COS(Stack(StP))**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=TAN(Stack(StP))
       CASE(cAsinh)
            d=one/SQRT(Stack(StP)**2+one)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ASINH(Stack(StP))
       CASE(cAcosh)
            IF(Stack(StP)<one) THEN
               EvalfErrType=5
               res=zero
               RETURN
            ENDIF
            IF(Stack(StP)==-one) THEN
               EvaldErrType=3
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=4
               res=zero
               RETURN
            END IF
            d=one/SQRT(Stack(StP)**2-one)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ACOSH(Stack(StP))
       CASE(cAtanh)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=6
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=7
               res=zero
               RETURN
            ELSE IF(Stack(StP)==-one) THEN
               EvaldErrType=5
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=6
               res=zero
               RETURN
            ENDIF
            d=one/(one-Stack(StP)**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ATANH(Stack(StP))
       CASE(cAsin)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=8
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=9
               res=zero
               RETURN
            ELSE IF(Stack(StP)==-one) THEN
               EvaldErrType=7
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=8
               res=zero
               RETURN
            ENDIF
            d=one/SQRT(one-Stack(StP)**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ASIN(Stack(StP))
       CASE(cAcos)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=10
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=11
               res=zero
               RETURN
            ELSE IF(Stack(StP)==-one) THEN
               EvaldErrType=9
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=10
               res=zero
               RETURN
            ENDIF
            d=-one/SQRT(one-Stack(StP)**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ACOS(Stack(StP))
       CASE(cAtan)
            d=one/(one+Stack(StP)**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ATAN(Stack(StP))
       CASE(cErfc)
            d=-two/SQRT(pi)*EXP(-Stack(StP)**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ERFC(Stack(StP))
       CASE(cErf)
            d=two/SQRT(pi)*EXP(-Stack(StP)**2)
            DStack(1,StP)=d*DStack(1,StP)
            Stack(StP)=ERF(Stack(StP))
       CASE(cSgn)
            DStack(1,StP)=zero
            Stack(StP)=SIGN(one,Stack(StP))
       CASE(cSign)
            DStack(1,StP-1)=zero
            Stack(StP-1)=SIGN(Stack(StP-1),Stack(StP))
       CASE(cLogb)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=12
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvalfErrType=13
               res=zero
               RETURN
            END IF
            IF(Stack(StP-1)<=zero) THEN
               EvalfErrType=14
               res=zero
               RETURN
            END IF
            DStack(1,StP-1)=( DStack(1,StP-1)/Stack(StP-1)*LOG(Stack(StP)) - &
                                    DStack(1,StP)/Stack(StP)*LOG(Stack(StP-1)) ) / &
                                    (LOG(Stack(StP))**2)
            Stack(StP-1)=LOG(Stack(StP-1))/LOG(Stack(StP))
            StP=StP-1
       CASE(cSum)
            DStack(1,StP-Comp(i)%NArg(InP)+1)=SUM(DStack(1,StP-Comp(i)%NArg(InP)+1:StP))
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cAve, cMean, cAmean)
            DStack(1,StP-Comp(i)%NArg(InP)+1)=SUM(DStack(1,StP-Comp(i)%NArg(InP)+1:StP))/REAL(Comp(i)%NArg(InP),dp)
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(Stack(StP-Comp(i)%NArg(InP)+1:StP))/REAL(Comp(i)%NArg(InP),dp)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cProd)
            tmp1=PRODUCT(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            tmp2=zero
            DO n=StP-Comp(i)%NArg(InP)+1,StP
               tmp2=tmp2+DStack(1,n)/Stack(n)*tmp1
            END DO
            DStack(1,StP-Comp(i)%NArg(InP)+1)=tmp2
            Stack(StP-Comp(i)%NArg(InP)+1)=tmp1
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cGmean)
            tmp1=PRODUCT(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            IF( tmp1<zero ) THEN
                EvalfErrType=15
                res=zero
                RETURN
            ENDIF
            tmp2=zero
            DO n=StP-Comp(i)%NArg(InP)+1,StP
               tmp2=tmp2+DStack(1,n)/Stack(n)*tmp1
            END DO
            DStack(1,StP-Comp(i)%NArg(InP)+1)=one/REAL(Comp(i)%NArg(InP),dp)*tmp1**(one/REAL(Comp(i)%NArg(InP),dp)-one)*tmp2
            Stack(StP-Comp(i)%NArg(InP)+1)=tmp1
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cHmean)
            IF( ANY(Stack(StP-Comp(i)%NArg(InP)+1:StP)==zero) ) THEN
                EvalfErrType=16
                res=zero
                RETURN
            END IF
            tmp1=SUM(one/Stack(StP-Comp(i)%NArg(InP)+1:StP))
            IF( tmp1==zero ) THEN
                EvalfErrType=17
                res=zero
                RETURN
            END IF
            tmp2=zero
            DO n=StP-Comp(i)%NArg(InP)+1,StP
               tmp2=tmp2+DStack(1,n)/(Stack(n)**2)
            END DO
            DStack(1,StP-Comp(i)%NArg(InP)+1)=REAL(Comp(i)%NArg(InP),dp)/(tmp1**2)*tmp2
            Stack(StP-Comp(i)%NArg(InP)+1)=REAL(Comp(i)%NArg(InP),dp)/tmp1
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cMin)
            t=MINLOC(Stack(StP-Comp(i)%NArg(InP)+1:StP),DIM=1)
            DStack(1,StP-Comp(i)%NArg(InP)+1)=DStack(1,StP-Comp(i)%NArg(InP)+t)
            Stack(StP-Comp(i)%NArg(InP)+1)=Stack(StP-Comp(i)%NArg(InP)+t)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cMax)
            t=MAXLOC(Stack(StP-Comp(i)%NArg(InP)+1:StP),DIM=1)
            DStack(1,StP-Comp(i)%NArg(InP)+1)=DStack(1,StP-Comp(i)%NArg(InP)+t)
            Stack(StP-Comp(i)%NArg(InP)+1)=Stack(StP-Comp(i)%NArg(InP)+t)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cWrap)
            DStack(1,StP-1)=DStack(1,StP-1)
            Stack(StP-1)=Stack(StP-1)-NINT(Stack(StP-1)/Stack(StP))*Stack(StP)
            StP=StP-1
       CASE(cDot)
            DStack(1,StP-5)=dotd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                       DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=dotf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossx)
            DStack(1,StP-5)=crossxd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                          DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=crossxf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossy)
            DStack(1,StP-5)=crossyd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                          DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=crossyf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossz)
            DStack(1,StP-5)=crosszd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                          DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=crosszf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cNorm)
            tmp1=normf(Stack(StP-2:StP))
            IF(tmp1==zero) THEN
               EvaldErrType=11
               res=zero
               RETURN
            END IF
            DStack(1,StP-2)=normd(Stack(StP-2:StP),DStack(1,StP-2:StP))
            Stack(StP-2)=tmp1
            StP=StP-2
       CASE(cDis)
            tmp1=disf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            IF(tmp1==zero) THEN
               EvaldErrType=12
               res=zero
               RETURN
            END IF
            DStack(1,StP-5)=disd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                       DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=tmp1
            StP=StP-5
       CASE(cAng)
            tmp1=angf(Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP))
            IF(tmp1==-pi) THEN
               EvaldErrType=13
               res=zero
               RETURN
            ELSE IF(tmp1==pi) THEN
               EvaldErrType=14
               res=zero
               RETURN
            END IF
            DStack(1,StP-8)=angd(Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                       DStack(1,StP-8:StP-6),DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-8)=tmp1
            StP=StP-8
       CASE(cDih)
            tmp1=dihf(Stack(StP-11:StP-9),Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP))

            DStack(1,StP-11)=dihd(Stack(StP-11:StP-9),Stack(StP-8:StP-6), &
                                        Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                        DStack(1,StP-11:StP-9),DStack(1,StP-8:StP-6), &
                                        DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-11)=tmp1
            StP=StP-11
       CASE  DEFAULT
            StP=StP+1
            IF(Comp(i)%ByteCode(InP)-VarBegin+1 == j) THEN
               DStack(1,StP)=one
            ELSE
               DStack(1,StP)=zero
            ENDIF
            Stack(StP)=Val(Comp(i)%ByteCode(InP)-VarBegin+1)
       END SELECT
    END DO
    res = DStack(1,1)
  END FUNCTION evald1
  !
  FUNCTION evaldm(Comp, i, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function's all derivatives for the values passed in array Val(:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                INTENT(in) :: i                  ! Function identifier
    REAL(dp), DIMENSION(:), INTENT(in) :: Val                ! Variable values
    REAL(dp), DIMENSION(SIZE(Val))     :: res                ! Result
    INTEGER                            :: m
    !----- -------- --------- --------- --------- --------- --------- --------- -------

    DO m=1, SIZE(Val)
       res(m) = evald1(Comp, i, m, Val(:))
    END DO
  END FUNCTION evaldm
  !
  FUNCTION evaldn(Comp, i, j, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function's jth derivative for the values passed in array Val(:,:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                  INTENT(in) :: i                  ! Function identifier
    INTEGER,                  INTENT(in) :: j                  ! Variable identifier
    REAL(dp), DIMENSION(:,:), INTENT(in) :: Val                ! Variable values
    REAL(dp), DIMENSION(SIZE(Val,2))     :: res                ! Result
    INTEGER                              :: n
    !----- -------- --------- --------- --------- --------- --------- --------- -------

    DO n=1, SIZE(Val,2)
       res(n) = evald1(Comp, i, j, Val(:,n))
    END DO
  END FUNCTION evaldn
  !
  FUNCTION evaldmn(Comp, i, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function's all derivatives for the values passed in array Val(:,:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                  INTENT(in)         :: i          ! Function identifier
    REAL(dp), DIMENSION(:,:), INTENT(in)         :: Val        ! Variable values
    REAL(dp), DIMENSION(SIZE(Val,1),SIZE(Val,2)) :: res        ! Result
    INTEGER                              :: n
    INTEGER                              :: m
    !----- -------- --------- --------- --------- --------- --------- --------- -------

    DO m=1, SIZE(Val,1)
       DO n=1, SIZE(Val,2)
          res(m,n) = evald1(Comp, i, m, Val(:,n))
       END DO
    END DO
  END FUNCTION evaldmn
  !
  FUNCTION evaldd1(Comp, i, j, k, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function's jth-kth second derivative for the values passed in array Val(:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(IN) :: Comp
    INTEGER,                INTENT(in) :: i                  ! Function identifier
    INTEGER,                INTENT(in) :: j                  ! Variable identifier
    INTEGER,                INTENT(in) :: k                  ! Variable identifier
    REAL(dp), DIMENSION(:), INTENT(in) :: Val                ! Variable values
    REAL(dp)                           :: res                ! Result
    INTEGER                            :: InP,             & ! Instruction pointer
                                          DaP,             & ! Data pointer
                                          StP                ! Stack pointer
    INTEGER                            :: n              
    INTEGER                            :: t
    REAL(dp)                           :: tmp1, tmp2, tmp3
    REAL(dp)                           :: d, dd
    REAL(dp)                           :: Stack(Comp(i)%StackSize)
    REAL(dp)                           :: DStack(2,Comp(i)%StackSize)
    REAL(dp)                           :: DDStack(Comp(i)%StackSize)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF((Comp(i)%Var(j) == 0) .OR. (Comp(i)%Var(k) == 0)) THEN
       res = zero
       RETURN
    END IF
    EvalfErrType = 0
    EvaldErrType = 0
    DaP = 1
    StP = 0
    DO InP=1,Comp(i)%ByteCodeSize
       SELECT CASE (Comp(i)%ByteCode(InP))

       CASE(cImmed)
            StP=StP+1
            DDStack(StP)=zero
            DStack(:,StP)=zero
            Stack(StP)=Comp(i)%Immed(DaP)
            DaP=DaP+1
       CASE(cNeg)
            DDStack(StP)=-DDStack(StP)
            DStack(:,StP)=-DStack(:,StP)
            Stack(StP)=-Stack(StP)
       CASE(cAdd)
            DDStack(StP-1)=DDStack(StP-1)+DDStack(StP)
            DStack(:,StP-1)=DStack(:,StP-1)+DStack(:,StP)
            Stack(StP-1)=Stack(StP-1)+Stack(StP)
            StP=StP-1
       CASE(cSub)
            DDStack(StP-1)=DDStack(StP-1)-DDStack(StP)
            DStack(:,StP-1)=DStack(:,StP-1)-DStack(:,StP)
            Stack(StP-1)=Stack(StP-1)-Stack(StP)
            StP=StP-1
       CASE(cMul)
            DDStack(StP-1)=DDStack(StP-1)*Stack(StP)+DStack(1,StP-1)*DStack(2,StP)+ &
                                   DStack(2,StP-1)*DStack(1,StP)+Stack(StP-1)*DDStack(StP)
            DStack(:,StP-1)=DStack(:,StP-1)*Stack(StP)+Stack(StP-1)*DStack(:,StP)
            Stack(StP-1)=Stack(StP-1)*Stack(StP)
            StP=StP-1
       CASE(cDiv)
            IF(Stack(StP)==zero) THEN
               EvalfErrType=1
               res=zero
               RETURN
            ENDIF
            DDStack(StP-1)=(Stack(StP)**2* &
                                    (Stack(StP)*DDStack(StP-1)- &
                                     DStack(1,StP-1)*DStack(2,StP)- &
                                     DStack(2,StP-1)*DStack(1,StP)- &
                                     Stack(StP-1)*DDStack(StP))+ &
                                    two*Stack(StP-1)*Stack(StP)*DStack(1,StP)*DStack(2,StP))/ &
                                    Stack(StP)**4
            DStack(:,StP-1)=( DStack(:,StP-1)*Stack(StP)-Stack(StP-1)*DStack(:,StP) ) / &
                                 ( Stack(StP)**2 )
            Stack(StP-1)=Stack(StP-1)/Stack(StP)
            StP=StP-1
       CASE(cPow)
            IF(Stack(StP-1)==zero) THEN
               IF(Stack(StP)<=zero) THEN
                  EvalfErrType=1
                  res=zero
                  RETURN
               ELSE IF(Stack(StP)==two) THEN ! special case when exponent=2
                  DDStack(StP-1)=two*DStack(1,StP-1)*DStack(2,StP-1)
                  DStack(:,StP-1)=zero
                  Stack(StP-1)=zero
               ELSE
                  DDStack(StP-1)=zero
                  DStack(:,StP-1)=zero
                  Stack(StP-1)=zero
               ENDIF
            ELSE IF(Stack(StP-1)<zero) THEN
               tmp1=Stack(StP-1)**NINT(Stack(StP))
               tmp2=( DStack(1,StP)*LOG(-Stack(StP-1)) + &
                    Stack(StP)/Stack(StP-1)*DStack(1,StP-1) )
               tmp3=( DStack(2,StP)*LOG(-Stack(StP-1)) + &
                    Stack(StP)/Stack(StP-1)*DStack(2,StP-1) )
               DDStack(StP-1)=tmp1*(tmp2*tmp3+DDStack(StP)*LOG(-(Stack(StP-1)))+ &
                                            (DStack(1,StP-1)*DStack(2,StP)+ &
                                             DStack(2,StP-1)*DStack(1,StP)- &
                                             Stack(StP)/Stack(StP-1)* &
                                             DStack(1,StP-1)*DStack(2,StP-1)+ &
                                             Stack(StP)*DDStack(StP-1))/Stack(StP-1))
               DStack(1,StP-1)=tmp1*tmp2
               DStack(2,StP-1)=tmp1*tmp3
               Stack(StP-1)=tmp1
            ELSE
               tmp1=Stack(StP-1)**Stack(StP)
               tmp2=( DStack(1,StP)*LOG(Stack(StP-1)) + &
                    Stack(StP)/Stack(StP-1)*DStack(1,StP-1) )
               tmp3=( DStack(2,StP)*LOG(Stack(StP-1)) + &
                    Stack(StP)/Stack(StP-1)*DStack(2,StP-1) )
               DDStack(StP-1)=tmp1*(tmp2*tmp3+DDStack(StP)*LOG(Stack(StP-1))+ &
                                            (DStack(1,StP-1)*DStack(2,StP)+ &
                                             DStack(2,StP-1)*DStack(1,StP)- &
                                             Stack(StP)/Stack(StP-1)* &
                                             DStack(1,StP-1)*DStack(2,StP-1)+ &  
                                             Stack(StP)*DDStack(StP-1))/Stack(StP-1))
               DStack(1,StP-1)=tmp1*tmp2
               DStack(2,StP-1)=tmp1*tmp3
               Stack(StP-1)=tmp1
            ENDIF
            StP=StP-1
       CASE(cAbs)
            IF(Stack(StP)==zero) THEN
               EvaldErrType=1
               res=zero
               RETURN
            ENDIF
            d=Stack(StP)/ABS(Stack(StP))
            DDStack(StP)=d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ABS(Stack(StP))
       CASE(cExp)
            tmp1=EXP(Stack(StP))
            d=tmp1
            dd=tmp1
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cLog10)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=2
               res=zero
               RETURN
            ENDIF
            d=one/(LOG(ten)*Stack(StP))
            dd=-d/Stack(StP)
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=LOG10(Stack(StP))
       CASE(cLog)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=3
               res=zero
               RETURN
            ENDIF
            d=one/Stack(StP)
            dd=-d/Stack(StP)
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=LOG(Stack(StP))
       CASE(cSqrt)
            IF(Stack(StP)<zero) THEN
               EvalfErrType=4
               res=zero
               RETURN
            ELSE IF(Stack(StP)==zero) THEN
               EvaldErrType=2
               res=zero
               RETURN
            ENDIF
            tmp1=SQRT(Stack(StP))
            d=one/tmp1/two
            dd=-d/(tmp1*tmp1*two)
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cSinh)
            tmp1=SINH(Stack(StP))
            d=COSH(Stack(StP))
            dd=tmp1
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cCosh)
            tmp1=COSH(Stack(StP))
            d=SINH(Stack(StP))
            dd=tmp1
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cTanh)
            tmp1=TANH(Stack(StP))
            d=one-tmp1**2
            dd=-two*tmp1*d
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cSin)
            tmp1=SIN(Stack(StP))
            d=COS(Stack(StP))
            dd=-tmp1
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cCos)
            tmp1=COS(Stack(StP))
            d=-SIN(Stack(StP))
            dd=-tmp1
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cTan)
            tmp1=TAN(Stack(StP))
            tmp2=COS(Stack(StP))
            d=one/(tmp2**2)
            dd=two*tmp1*d
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=tmp1
       CASE(cAsinh)
            d=one/SQRT(Stack(StP)**2+one)
            dd=-Stack(StP)*d**3
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ASINH(Stack(StP))
       CASE(cAcosh)
            IF(Stack(StP)<one) THEN
               EvalfErrType=5
               res=zero
               RETURN
            ENDIF
            IF(Stack(StP)==-one) THEN
               EvaldErrType=3
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=4
               res=zero
               RETURN
            END IF
            d=one/SQRT(Stack(StP)**2-one)
            dd=-Stack(StP)*d**3
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ACOSH(Stack(StP))
       CASE(cAtanh)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=6
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=7
               res=zero
               RETURN
            ELSE IF(Stack(StP)==-one) THEN
               EvaldErrType=5
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=6
               res=zero
               RETURN
            ENDIF
            d=one/(one-Stack(StP)**2)
            dd=two*Stack(StP)*d**2
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ATANH(Stack(StP))
       CASE(cAsin)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=8
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=9
               res=zero
               RETURN
            ELSE IF(Stack(StP)==-one) THEN
               EvaldErrType=7
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=8
               res=zero
               RETURN
            ENDIF
            d=one/SQRT(one-Stack(StP)**2)
            dd=Stack(StP)*d**3
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ASIN(Stack(StP))
       CASE(cAcos)
            IF(Stack(StP)<-one) THEN
               EvalfErrType=10
               res=zero
               RETURN
            ELSE IF(Stack(StP)>one) THEN
               EvalfErrType=11
               res=zero
               RETURN
            ELSE IF(Stack(StP)==-one) THEN
               EvaldErrType=9
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvaldErrType=10
               res=zero
               RETURN
            ENDIF
            d=-one/SQRT(one-Stack(StP)**2)
            dd=Stack(StP)*d**3
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ACOS(Stack(StP))
       CASE(cAtan)
            d=one/(one+Stack(StP)**2)
            dd=-two*Stack(StP)*d**2
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ATAN(Stack(StP))
       CASE(cErfc)
            d=-two/SQRT(pi)*EXP(-Stack(StP)**2)
            dd=-two*Stack(StP)*d
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ERFC(Stack(StP))
       CASE(cErf)
            d=two/SQRT(pi)*EXP(-Stack(StP)**2)
            dd=-two*Stack(StP)*d
            DDStack(StP)=dd*DStack(1,StP)*DStack(2,StP)+d*DDStack(StP)
            DStack(:,StP)=d*DStack(:,StP)
            Stack(StP)=ERF(Stack(StP))
       CASE(cSgn)
            DDStack(StP)=zero
            DStack(:,StP)=zero
            Stack(StP)=SIGN(one,Stack(StP))
       CASE(cSign)
            DDStack(StP-1)=zero
            DStack(:,StP-1)=zero
            Stack(StP-1)=SIGN(Stack(StP-1),Stack(StP))
       CASE(cLogb)
            IF(Stack(StP)<=zero) THEN
               EvalfErrType=12
               res=zero
               RETURN
            ELSE IF(Stack(StP)==one) THEN
               EvalfErrType=13
               res=zero
               RETURN
            END IF
            IF(Stack(StP-1)<=zero) THEN
               EvalfErrType=14
               res=zero
               RETURN
            END IF
            tmp1=LOG(Stack(StP))
            tmp2=LOG(Stack(StP-1))
            DDStack(StP-1)=( (-DStack(1,StP-1)*tmp1/Stack(StP-1)**2 &
                                      +DDStack(StP-1)*tmp1/Stack(StP-1) &
                                      +DStack(1,StP-1)*DStack(2,StP)/Stack(StP-1)/Stack(StP) &
                                      -DStack(2,StP-1)*DStack(1,StP)/Stack(StP-1)/Stack(StP) &
                                      +DStack(1,StP)*DStack(2,StP)*tmp2/Stack(StP)**2 &
                                      +DDStack(StP)*tmp2/Stack(StP))*tmp1**2 &
                                     -(DStack(1,StP-1)*tmp1/Stack(StP-1) &
                                       -DStack(1,StP)*tmp2/Stack(StP)) &
                                      *two*tmp1*DStack(2,StP)/Stack(StP) ) / tmp1**4
            DDStack(StP-1)=( DDStack(StP-1)*tmp1**3 - tmp1* &
                                     (DStack(1,StP-1)*DStack(2,StP) - &
                                      DStack(2,StP-1)*DStack(1,StP) - &
                                      DDStack(StP)*tmp2) + &
                                     two*DStack(1,StP-1)*DStack(2,StP-1)*tmp1*tmp2 ) / tmp1**4
            DStack(:,StP-1)=( DStack(:,StP-1)/Stack(StP-1)*tmp1 - &
                                    DStack(:,StP)/Stack(StP)*tmp2 ) / tmp1**2
            Stack(StP-1)=tmp2/tmp1
            StP=StP-1
       CASE(cSum)
            DDStack(StP-Comp(i)%NArg(InP)+1)=SUM(DDStack(StP-Comp(i)%NArg(InP)+1:StP))
            DStack(:,StP-Comp(i)%NArg(InP)+1)=SUM(DStack(:,StP-Comp(i)%NArg(InP)+1:StP))
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cAve, cMean, cAmean)
            DDStack(StP-Comp(i)%NArg(InP)+1)=SUM(DDStack(StP-Comp(i)%NArg(InP)+1:StP))/REAL(Comp(i)%NArg(InP),dp)
            DStack(:,StP-Comp(i)%NArg(InP)+1)=SUM(DStack(:,StP-Comp(i)%NArg(InP)+1:StP))/REAL(Comp(i)%NArg(InP),dp)
            Stack(StP-Comp(i)%NArg(InP)+1)=SUM(Stack(StP-Comp(i)%NArg(InP)+1:StP))/REAL(Comp(i)%NArg(InP),dp)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cProd)
!LAM
            tmp1=PRODUCT(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            tmp2=zero
            DO n=StP-Comp(i)%NArg(InP)+1,StP
               tmp2=tmp2+DStack(1,n)/Stack(n)*tmp1
            END DO
            DStack(1,StP-Comp(i)%NArg(InP)+1)=tmp2
            Stack(StP-Comp(i)%NArg(InP)+1)=tmp1
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cGmean)
!LAM
            tmp1=PRODUCT(Stack(StP-Comp(i)%NArg(InP)+1:StP))
            IF( tmp1<zero ) THEN
                EvalfErrType=15
                res=zero
                RETURN
            ENDIF
            tmp2=zero
            DO n=StP-Comp(i)%NArg(InP)+1,StP
               tmp2=tmp2+DStack(1,n)/Stack(n)*tmp1
            END DO
            DStack(1,StP-Comp(i)%NArg(InP)+1)=one/REAL(Comp(i)%NArg(InP),dp)*tmp1**(one/REAL(Comp(i)%NArg(InP),dp)-one)*tmp2
            Stack(StP-Comp(i)%NArg(InP)+1)=tmp1
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cHmean)
!LAM
            IF( ANY(Stack(StP-Comp(i)%NArg(InP)+1:StP)==zero) ) THEN
                EvalfErrType=16
                res=zero
                RETURN
            END IF
            tmp1=SUM(one/Stack(StP-Comp(i)%NArg(InP)+1:StP))
            IF( tmp1==zero ) THEN
                EvalfErrType=17
                res=zero
                RETURN
            END IF
            tmp2=zero
            DO n=StP-Comp(i)%NArg(InP)+1,StP
               tmp2=tmp2+DStack(1,n)/(Stack(n)**2)
            END DO
            DStack(1,StP-Comp(i)%NArg(InP)+1)=REAL(Comp(i)%NArg(InP),dp)/(tmp1**2)*tmp2
            Stack(StP-Comp(i)%NArg(InP)+1)=REAL(Comp(i)%NArg(InP),dp)/tmp1
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cMin)
            t=MINLOC(Stack(StP-Comp(i)%NArg(InP)+1:StP),DIM=1)
            DDStack(StP-Comp(i)%NArg(InP)+1)=DDStack(StP-Comp(i)%NArg(InP)+t)
            DStack(:,StP-Comp(i)%NArg(InP)+1)=DStack(:,StP-Comp(i)%NArg(InP)+t)
            Stack(StP-Comp(i)%NArg(InP)+1)=Stack(StP-Comp(i)%NArg(InP)+t)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cMax)
            t=MAXLOC(Stack(StP-Comp(i)%NArg(InP)+1:StP),DIM=1)
            DDStack(StP-Comp(i)%NArg(InP)+1)=DDStack(StP-Comp(i)%NArg(InP)+t)
            DStack(:,StP-Comp(i)%NArg(InP)+1)=DStack(:,StP-Comp(i)%NArg(InP)+t)
            Stack(StP-Comp(i)%NArg(InP)+1)=Stack(StP-Comp(i)%NArg(InP)+t)
            StP=StP-Comp(i)%NArg(InP)+1
       CASE(cWrap)
            DDStack(StP-1)=DDStack(StP-1)
            DStack(:,StP-1)=DStack(:,StP-1)
            Stack(StP-1)=Stack(StP-1)-NINT(Stack(StP-1)/Stack(StP))*Stack(StP)
            StP=StP-1
       CASE(cDot)
!LAM
            DStack(1,StP-5)=dotd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                       DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=dotf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossx)
!LAM
            DStack(1,StP-5)=crossxd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                          DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=crossxf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossy)
!LAM
            DStack(1,StP-5)=crossyd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                          DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=crossyf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cCrossz)
!LAM
            DStack(1,StP-5)=crosszd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                          DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=crosszf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            StP=StP-5
       CASE(cNorm)
!LAM
            tmp1=normf(Stack(StP-2:StP))
            IF(tmp1==zero) THEN
               EvaldErrType=11
               res=zero
               RETURN
            END IF
            DStack(1,StP-2)=normd(Stack(StP-2:StP),DStack(1,StP-2:StP))
            Stack(StP-2)=tmp1
            StP=StP-2
       CASE(cDis)
!LAM
            tmp1=disf(Stack(StP-5:StP-3),Stack(StP-2:StP))
            IF(tmp1==zero) THEN
               EvaldErrType=12
               res=zero
               RETURN
            END IF
            DStack(1,StP-5)=disd(Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                       DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-5)=tmp1
            StP=StP-5
       CASE(cAng)
!LAM
            tmp1=angf(Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP))
            IF(tmp1==-pi) THEN
               EvaldErrType=13
               res=zero
               RETURN
            ELSE IF(tmp1==pi) THEN
               EvaldErrType=14
               res=zero
               RETURN
            END IF
            DStack(1,StP-8)=angd(Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                       DStack(1,StP-8:StP-6),DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-8)=tmp1
            StP=StP-8
       CASE(cDih)
!LAM
            tmp1=dihf(Stack(StP-11:StP-9),Stack(StP-8:StP-6),Stack(StP-5:StP-3),Stack(StP-2:StP))

            DStack(1,StP-11)=dihd(Stack(StP-11:StP-9),Stack(StP-8:StP-6), &
                                        Stack(StP-5:StP-3),Stack(StP-2:StP), &
                                        DStack(1,StP-11:StP-9),DStack(1,StP-8:StP-6), &
                                        DStack(1,StP-5:StP-3),DStack(1,StP-2:StP))
            Stack(StP-11)=tmp1
            StP=StP-11
       CASE  DEFAULT
            StP=StP+1
            DDStack(StP)=zero
            IF(Comp(i)%ByteCode(InP)-VarBegin+1 == j) THEN
               DStack(1,StP)=one
            ELSE
               DStack(1,StP)=zero
            ENDIF
            IF(Comp(i)%ByteCode(InP)-VarBegin+1 == k) THEN
               DStack(2,StP)=one
            ELSE
               DStack(2,StP)=zero
            ENDIF
            Stack(StP)=Val(Comp(i)%ByteCode(InP)-VarBegin+1)
       END SELECT
    END DO
    res = DDStack(1)
  END FUNCTION evaldd1
  !
  SUBROUTINE CheckSyntax (Func,FuncStr,Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check syntax of function string,  returns 0 if syntax is ok
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),               INTENT(in) :: Func      ! Function string without spaces
    CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Original function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n
    CHARACTER (LEN=1)                           :: c
    REAL(dp)                                    :: r
    LOGICAL                                     :: err
    INTEGER                                     :: ParCnt    ! Parenthesis counter
    INTEGER                                     :: ArgCnt    ! Argument counter
    INTEGER                                     :: j,k,ib,in,lFunc,pc
    CHARACTER (LEN=length)                      :: ch1, ch2
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    j = 1
    ParCnt = 0
    lFunc = LEN_TRIM(Func)
    step: DO
       IF (j > lFunc) CALL ParseErrMsg (j, FuncStr)
       c = Func(j:j)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Check for valid operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (c == '-' .OR. c == '+') THEN                      ! Check for leading - or +
          j = j+1
          IF (j > lFunc) CALL ParseErrMsg (j-1, FuncStr, 'Missing operand')
          c = Func(j:j)
          IF (ANY(c == Ops)) CALL ParseErrMsg (j, FuncStr, 'Multiple operators')
          IF (c == ',') CALL ParseErrMsg (j, FuncStr, 'Missing argument') ! lam81
       END IF
       n = MathFunctionIndex (Func(j:))
       IF (n > 0) THEN                                       ! Check for math function
          j = j+LEN_TRIM(Funcs(n))
          IF (j > lFunc) CALL ParseErrMsg (j-1, FuncStr, 'Missing function argument')
          c = Func(j:j)
          IF (c /= '(') CALL ParseErrMsg (j, FuncStr, 'Missing opening parenthesis')
          ! check number of arguments
          ArgCnt = 1         ! lam81
          k=j                ! lam81
          pc = 1             ! lam81
          DO WHILE (pc /= 0) ! lam81
             k=k+1                             ! lam81
             IF (k > lFunc) CALL ParseErrMsg (k-1, FuncStr, 'Missing closing parenthesis')
             IF (Func(k:k) == '(') THEN        ! lam81
                pc = pc + 1                    ! lam81
             ELSE IF (Func(k:k) == ')') THEN   ! lam81
                pc = pc - 1                    ! lam81
             ELSE IF (Func(k:k) == ',' .AND. pc==1) THEN   ! lam81
                ArgCnt = ArgCnt + 1            ! lam81
             END IF                            ! lam81
          END DO                               ! lam81
          IF (NArgFuncs(n) > 0) THEN           ! lam81
             IF (ArgCnt /= NArgFuncs(n)) THEN  ! lam81
                 WRITE(ch1,*) ArgCnt           ! lam81
                 WRITE(ch2,*) NArgFuncs(n)     ! lam81
                 CALL ParseErrMsg (k, FuncStr, 'Number of specified argument('//TRIM(ADJUSTL(ch1))// &
                                   ') of '//TRIM(Funcs(n))//' is not equal to '//TRIM(ADJUSTL(ch2)))
             END IF                            ! lam81
          ELSE IF (NArgFuncs(n) < 0) THEN      ! lam81
             IF (ArgCnt < -NArgFuncs(n)) THEN  ! lam81
                 WRITE(ch1,*) ArgCnt           ! lam81
                 WRITE(ch2,*) NArgFuncs(n)     ! lam81
                 CALL ParseErrMsg (k, FuncStr, 'Number of specified argument('//TRIM(ADJUSTL(ch1))// &
                                   ') of '//TRIM(Funcs(n))//' is less than the minimum of '//TRIM(ADJUSTL(ch2)))
             END IF                            ! lam81
          END IF                               ! lam81
       END IF
       IF (c == '(') THEN                                    ! Check for opening parenthesis
          ParCnt = ParCnt+1
          j = j+1
          c = Func(j:j) ! lam81
          IF (c == ',') CALL ParseErrMsg (j, FuncStr, 'Missing argument') ! lam81
          CYCLE step
       END IF
       IF (SCAN(c,'0123456789.') > 0) THEN                   ! Check for number
          r = RealNum (Func(j:),ib,in,err)
          IF (err) CALL ParseErrMsg (j, FuncStr, 'Invalid number format:  '//Func(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       ELSE                                                  ! Check for variable
          n = VariableIndex (Func(j:),Var,ib,in)
          IF (n == 0) CALL ParseErrMsg (j, FuncStr, 'Invalid element: '//Func(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       END IF
       DO WHILE (c == ')')                                   ! Check for closing parenthesis
          ParCnt = ParCnt-1
          IF (ParCnt < 0) CALL ParseErrMsg (j, FuncStr, 'Mismatched parenthesis')
          IF (Func(j-1:j-1) == '(') CALL ParseErrMsg (j-1, FuncStr, 'Empty parentheses')
          j = j+1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       END DO
       IF (c == ',') THEN  ! lam81
          j = j+1          ! lam81
          c = Func(j:j)    ! lam81
          IF (j > lFunc) CALL ParseErrMsg (j-1, FuncStr, 'Missing argument') ! lam81
          IF (c == ',') CALL ParseErrMsg (j, FuncStr, 'Missing argument') ! lam81
          CYCLE step       ! lam81
       END IF              ! lam81
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have a legal operand: A legal operator or end of string must follow
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (j > lFunc) EXIT
       IF (ANY(c == Ops)) THEN                               ! Check for multiple operators
          IF (j+1 > lFunc) CALL ParseErrMsg (j, FuncStr)
          IF (ANY(Func(j+1:j+1) == Ops)) CALL ParseErrMsg (j+1, FuncStr, 'Multiple operators')
       ELSE                                                  ! Check for next operand
          CALL ParseErrMsg (j, FuncStr, 'Missing operator')
       END IF
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have an operand and an operator: the next loop will check for another 
       ! operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       j = j+1
    END DO step
    IF (ParCnt > 0) CALL ParseErrMsg (j, FuncStr, 'Missing )')
  END SUBROUTINE CheckSyntax
  !
  FUNCTION EvalfErrMsg () RESULT (msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return error message
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=40)                     :: msg
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (EvalfErrType == 0) THEN
       msg = ''
    ELSE
       msg = TRIM(EvalfErrText(EvalfErrType))
    ENDIF
  END FUNCTION EvalfErrMsg
    !
  FUNCTION EvaldErrMsg () RESULT (msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return error message
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=40)                     :: msg
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (EvaldErrType == 0) THEN
       msg = ''
    ELSE
       msg = TRIM(EvaldErrText(EvaldErrType))
    ENDIF
  END FUNCTION EvaldErrMsg
  !
  SUBROUTINE ParseErrMsg (j, FuncStr, Msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Print error message and terminate program
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                     INTENT(in) :: j
    CHARACTER (LEN=*),           INTENT(in) :: FuncStr       ! Original function string
    CHARACTER (LEN=*), OPTIONAL, INTENT(in) :: Msg
    INTEGER                                 :: k
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (PRESENT(Msg)) THEN
       WRITE(*,*) '*** Error in syntax of function string: '//TRIM(Msg)
    ELSE
       WRITE(*,*) '*** Error in syntax of function string:'
    ENDIF
    WRITE(*,*)
    WRITE(*,'(A)') ' '//TRIM(FuncStr)
    DO k=1,ipos(j)
       WRITE(*,'(A)',ADVANCE='NO') ' '                       ! Advance to the jth position
    END DO
    WRITE(*,'(A)') '?'
    STOP
  END SUBROUTINE ParseErrMsg
  !
  FUNCTION OperatorIndex (c) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return operator index
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=1), INTENT(in) :: c
    INTEGER(is)                   :: n,j
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    DO j=cAdd,cPow
       IF (c == Ops(j)) THEN
          n = j
          EXIT
       END IF
    END DO
  END FUNCTION OperatorIndex
  !
  FUNCTION MathFunctionIndex (str) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of math function beginnig at 1st position of string str
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(in) :: str
    INTEGER(is)                   :: n,j
    INTEGER                       :: k
    CHARACTER (LEN=LEN(Funcs))    :: fun
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    ! find the end of possible function string
    DO k=2,LEN(str)
       IF (SCAN(str(k:k),'+-*/^(), ') > 0) EXIT
    END DO
    CALL LowCase (str(1:k-1), fun)
    DO j=cPow+1,VarBegin-1                                   ! Check all math functions
       IF (fun == Funcs(j)) THEN                             ! Compare lower case letters
          n = j                                              ! Found a matching function         
          EXIT
       END IF 
    END DO
!   DO j=cPow+1,VarBegin-1                                   ! Check all math functions
!      k = MIN (LEN_TRIM(Funcs(j)), LEN(str))   
!      CALL LowCase (str(1:k), fun)
!      IF (fun  == Funcs(j)) THEN                             ! Compare lower case letters
!         n = j                                              ! Found a matching function
!         EXIT
!      END IF
!   END DO
  END FUNCTION MathFunctionIndex
  !
  FUNCTION VariableIndex (str, Var, ibegin, inext) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of variable at begin of string str (returns 0 if no variable found)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),               INTENT(in) :: str       ! String
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n         ! Index of variable
    INTEGER, OPTIONAL,              INTENT(out) :: ibegin, & ! Start position of variable name
                                                   inext     ! Position of character after name
    INTEGER                                     :: j,ib,in,lstr
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    lstr = LEN_TRIM(str)
    IF (lstr > 0) THEN
       DO ib=1,lstr                                          ! Search for first character in str
          IF (str(ib:ib) /= ' ') EXIT                        ! When lstr>0 at least 1 char in str
       END DO                        
       DO in=ib,lstr                                         ! Search for name terminators
!         IF (SCAN(str(in:in),'+-*/^) ') > 0) EXIT  ! lam81
          IF (SCAN(str(in:in),'+-*/^), ') > 0) EXIT ! lam81
       END DO
       DO j=1,SIZE(Var)
          IF (str(ib:in-1) == Var(j)) THEN                     
             n = j                                           ! Variable name found
             EXIT
          END IF
       END DO
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
  END FUNCTION VariableIndex
  !
  SUBROUTINE RemoveSpaces (str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Remove Spaces from string, remember positions of characters in old string
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(inout) :: str
    INTEGER                          :: k,lstr
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    lstr = LEN_TRIM(str)
    ipos = (/ (k,k=1,lstr) /)
    k = 1
    DO WHILE (str(k:lstr) /= ' ')                             
       IF (str(k:k) == ' ') THEN
          str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
          ipos(k:lstr) = (/ ipos(k+1:lstr), 0 /)             ! Move 1 element to left
          k = k-1
       END IF
       k = k+1
    END DO
  END SUBROUTINE RemoveSpaces
  !
  SUBROUTINE RemoveAllSpaces (str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Remove Spaces from string
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(inout) :: str
    INTEGER                          :: k,lstr
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    lstr = LEN_TRIM(str)
    k = 1
    DO WHILE (str(k:lstr) /= ' ')
       IF (str(k:k) == ' ') THEN
          str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
          k = k-1
       END IF
       k = k+1
    END DO
  END SUBROUTINE RemoveAllSpaces
  !
  SUBROUTINE Replace (ca,cb,str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Replace ALL appearances of character set ca in string str by character set cb
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),       INTENT(in) :: ca
    CHARACTER (LEN=LEN(ca)), INTENT(in) :: cb                ! LEN(ca) must be LEN(cb)
    CHARACTER (LEN=*),    INTENT(inout) :: str
    INTEGER                             :: j,lca
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    lca = LEN(ca)
    DO j=1,LEN_TRIM(str)-lca+1
       IF (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    END DO
  END SUBROUTINE Replace
  !
  SUBROUTINE CheckVariables(Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check variable names: check for conflict with function names
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), DIMENSION(:), INTENT(inout) :: Var       ! Array with variable names
    INTEGER                             :: m
    INTEGER                             :: i, j

    m=size(Var)

    DO i=1,m
       CALL RemoveAllSpaces(Var(i))
    END DO

    DO i=1,m
       IF(SCAN(TRIM(Var(i)),'+-*/^(),') > 0) THEN
          WRITE(*,*) '*** Error: var ', TRIM(Var(i)), ' must not include the following charaters: +-*/^(),'
             STOP
       END IF
    END DO

    DO i=1,m
       DO j=cPow+1,VarBegin-1
          IF(TRIM(Var(i)) == TRIM(Funcs(j))) THEN
             WRITE(*,*) '*** Error: var ', TRIM(Var(i)), ' must not be a math function name!'
             STOP
          END IF
       END DO
    END DO

  END SUBROUTINE CheckVariables
  !
  SUBROUTINE CheckVarExpressions (Varname, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check expression variable names: check for individuality and conflict with function and variable names
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), DIMENSION(:), INTENT(inout) :: Varname   ! Array with expression variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(inout) :: Var       ! Array with variable names
    INTEGER                             :: n, m
    INTEGER                             :: i, j

    n=size(Varname)
    m=size(Var)

    DO i=1,n
       CALL RemoveAllSpaces(Varname(i))
    END DO

    DO i=1,m
       CALL RemoveAllSpaces(Var(i))
    END DO

    DO i=1,n
       IF(SCAN(TRIM(Varname(i)),'+-*/^(),') > 0) THEN
          WRITE(*,*) '*** Error: varname ', TRIM(Varname(i)), ' must not include the following charaters: +-*/^(),'
             STOP
       END IF
    END DO

    ! check for individuality
    DO i=1,n-1
       DO j=i+1,n
          IF(TRIM(Varname(i)) == TRIM(Varname(j))) THEN
             WRITE(*,*) '*** Error: varname ', TRIM(Varname(i)), ' is defined more than once!'
             STOP
          END IF
       END DO
    END DO

    ! check for conflict with function names
    DO i=1,n
       DO j=cPow+1,VarBegin-1
          IF(TRIM(Varname(i)) == TRIM(Funcs(j))) THEN
             WRITE(*,*) '*** Error: varname ', TRIM(Varname(i)), ' must not be a math function name!'
             STOP
          END IF
       END DO
    END DO

    ! check for conflict with variable names
    DO i=1,n
       DO j=1,m
          IF(TRIM(Varname(i)) == TRIM(Funcs(j))) THEN
             WRITE(*,*) '*** Error: expression varname ', TRIM(Varname(i)), ' must not be a variable name!'
             STOP
          END IF
       END DO
    END DO 
    
  END SUBROUTINE CheckVarExpressions
  !
  SUBROUTINE InsertVarExpressions (Varname, Varexpr, str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Replace ALL appearances of variable names in string str by variable expressions
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),    INTENT(inout) :: Varname(:)
    CHARACTER (LEN=*),    INTENT(inout) :: Varexpr(:)
    CHARACTER (LEN=*),    INTENT(inout) :: str(:)
    INTEGER                             :: n, m
    INTEGER                             :: i, j, k, l
    INTEGER                             :: b, e
    CHARACTER (LEN=length)              :: ch
    !----- -------- --------- --------- --------- --------- --------- --------- -------

    IF (size(varname) /= size(varexpr) ) THEN
        WRITE(*,*) '*** Error: size(varname)= ', size(Varname), ' is not equal to size(varexpr)= ', size(Varexpr)
        STOP
    END IF

    n=size(Varname)
    m=size(str)

    ! remove spaces from varnames and varexpressions
    DO i=1,n
       CALL RemoveAllSpaces(Varname(i))
       CALL RemoveAllSpaces(Varexpr(i))
    END DO

    DO i=1,m
       CALL RemoveAllSpaces(str(i)) 
    END DO

    DO i=1,m
       b=1
       step1: DO
          IF(b>LEN_TRIM(str(i))) EXIT step1
          IF(SCAN(str(i)(b:b),'+-*/^(),') == 0) THEN
             e=b
             step2: DO
                IF(e>LEN_TRIM(str(i))+1) EXIT step1
                IF( (SCAN(str(i)(e:e),'+-*/^(),') > 0) .or. (e == LEN_TRIM(str(i))+1) ) THEN
                   DO j=1,n
                      IF(TRIM(Varname(j)) == str(i)(b:e-1)) THEN
                         k=LEN_TRIM(str(i))
                         l=LEN_TRIM(Varexpr(j))
                         IF( k+l-e+b+2 > LEN(str(i)) ) THEN
                             WRITE(*,*) '*** Error: function string length is short - there might be a cross referece!'
                             STOP
                         END IF
                         str(i)(b+l+2:k+l+2-e+b)=str(i)(e:k)
                         str(i)(b:b+l+1)='('//TRIM(Varexpr(j))//')'
                         e=b
                         EXIT
                      END IF
                      IF(j==n) THEN
                         b=e+1
                         EXIT step2
                      END IF
                   END DO
                ELSE
                   e=e+1
                END IF
             END DO step2
          ELSE
             b=b+1
          END IF 
       END DO step1
    END DO

  END SUBROUTINE InsertVarExpressions
  !
  SUBROUTINE Compile (Comp, i, F, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Compile i-th function string F into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(INOUT) :: Comp
    INTEGER,                         INTENT(in) :: i         ! Function identifier
    CHARACTER (LEN=*),               INTENT(in) :: F         ! Function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER                                     :: istat
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (ASSOCIATED(Comp(i)%ByteCode)) DEALLOCATE ( Comp(i)%ByteCode, &
                                                   Comp(i)%Immed,    &
                                                   Comp(i)%NArg,     &
                                                   Comp(i)%Var ) ! lam81

    Comp(i)%ByteCodeSize = 0
    Comp(i)%ImmedSize    = 0
    Comp(i)%StackSize    = 0
    Comp(i)%StackPtr     = 0
    CALL CompileSubstr (Comp,i,F,1,LEN_TRIM(F),Var)               ! Compile string to determine size
    ALLOCATE ( Comp(i)%ByteCode(Comp(i)%ByteCodeSize), & 
               Comp(i)%Immed(Comp(i)%ImmedSize),       &
               Comp(i)%NArg(Comp(i)%ByteCodeSize),     & ! lam81
               Comp(i)%Var(SIZE(Var)),                 & ! lam81
               STAT = istat                            )
    IF (istat /= 0) THEN
       WRITE(*,*) '*** Parser error: Memmory allocation for byte code failed'
       STOP
    ELSE
       Comp(i)%ByteCodeSize = 0
       Comp(i)%ImmedSize    = 0
       Comp(i)%StackSize    = 0
       Comp(i)%StackPtr     = 0
       Comp(i)%Var          = 0
       CALL CompileSubstr (Comp,i,F,1,LEN_TRIM(F),Var)            ! Compile string into bytecode
    END IF
    !
  END SUBROUTINE Compile
  !
  SUBROUTINE AddCompiledByte (Comp, i, b, narg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Add compiled byte to bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(INOUT) :: Comp
    INTEGER,     INTENT(in) :: i                             ! Function identifier  
    INTEGER(is), INTENT(in) :: b                             ! Value of byte to be added
    INTEGER,     INTENT(in), OPTIONAL :: narg                ! Number of function arguments
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    Comp(i)%ByteCodeSize = Comp(i)%ByteCodeSize + 1
    IF (ASSOCIATED(Comp(i)%ByteCode)) THEN
        Comp(i)%ByteCode(Comp(i)%ByteCodeSize) = b
        IF (PRESENT(narg)) THEN
            Comp(i)%NArg(Comp(i)%ByteCodeSize) = narg
        END IF
    END IF
  END SUBROUTINE AddCompiledByte
  !
  FUNCTION MathItemIndex (Comp, i, F, Var) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return math item index, if item is real number, enter it into Comp-structure
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(INOUT) :: Comp
    INTEGER,                         INTENT(in) :: i         ! Function identifier  
    CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n         ! Byte value of math item
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    IF (SCAN(F(1:1),'0123456789.') > 0) THEN                 ! Check for begin of a number
       Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
       IF (ASSOCIATED(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = RealNum (F)
       n = cImmed
    ELSE                                                     ! Check for a variable
       n = VariableIndex (F, Var)
       IF (n > 0) THEN
          IF (ASSOCIATED(Comp(i)%Var)) Comp(i)%Var(n) = 1    ! lam81
          n = VarBegin+n-1
       END IF
    END IF
  END FUNCTION MathItemIndex
  !
  FUNCTION CompletelyEnclosed (F, b, e) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(in) :: F                       ! Function substring
    INTEGER,           INTENT(in) :: b,e                     ! First and last pos. of substring
    LOGICAL                       :: res
    INTEGER                       :: j,k
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    res=.false.
    IF (F(b:b) == '(' .AND. F(e:e) == ')') THEN
       k = 0
       DO j=b+1,e-1
          IF     (F(j:j) == '(') THEN
             k = k+1
          ELSEIF (F(j:j) == ')') THEN
             k = k-1
          END IF
          IF (k < 0) EXIT
       END DO
       IF (k == 0) res=.true.                                ! All opened parenthesis closed
    END IF
  END FUNCTION CompletelyEnclosed
  !
  RECURSIVE SUBROUTINE CompileSubstr (Comp, i, F, b, e, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Compile i-th function string F into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    TYPE(tComp), DIMENSION(:), POINTER, INTENT(INOUT) :: Comp
    INTEGER,                         INTENT(in) :: i         ! Function identifier  
    CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
    INTEGER,                         INTENT(in) :: b,e       ! Begin and end position substring
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n
    INTEGER                                     :: b2,j,k,io,pc
    CHARACTER (LEN=*),                PARAMETER :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                                            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    INTEGER                                     :: ArgCnt    ! Argument Counter    ! lam81
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for special cases of substring
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF     (F(b:b) == '+') THEN                              ! Case 1: F(b:e) = '+...'
!      WRITE(*,*)'1. F(b:e) = "+..."'
       CALL CompileSubstr (Comp, i, F, b+1, e, Var)
       RETURN
    ELSEIF (CompletelyEnclosed (F, b, e)) THEN               ! Case 2: F(b:e) = '(...)'
!      WRITE(*,*)'2. F(b:e) = "(...)"'
       CALL CompileSubstr (Comp, i, F, b+1, e-1, Var)
       RETURN
    ELSEIF (SCAN(F(b:b),calpha) > 0) THEN        
       n = MathFunctionIndex (F(b:e))
       IF (n > 0) THEN
          b2 = b+INDEX(F(b:e),'(')-1
          IF (CompletelyEnclosed(F, b2, e)) THEN             ! Case 3: F(b:e) = 'fcn(...)'
!            WRITE(*,*)'3. F(b:e) = "fcn(...)"'
             ArgCnt=1                                         ! lam81
             k=b2+1                                           ! lam81
             pc=1                                             ! lam81
             DO WHILE (k < e)                                 ! lam81
                IF (F(k:k) == '(') THEN                       ! lam81
                    pc = pc + 1                               ! lam81
                ELSE IF (F(k:k) == ')') THEN                  ! lam81
                    pc = pc - 1                               ! lam81
                ELSE IF (F(k:k) == ',' .AND. pc==1) THEN      ! lam81
                    CALL CompileSubstr(Comp, i, F, b2+1, k-1, Var)  ! lam81
                    b2=k                                      ! lam81
                    ArgCnt = ArgCnt + 1                       ! lam81
                END IF                                        ! lam81
                k=k+1                                         ! lam81
             END DO                                           ! lam81
             CALL CompileSubstr(Comp, i, F, b2+1, e-1, Var)
             CALL AddCompiledByte (Comp, i, n, ArgCnt)        ! lam81
             RETURN
          END IF
       END IF
    ELSEIF (F(b:b) == '-') THEN
       IF (CompletelyEnclosed (F, b+1, e)) THEN              ! Case 4: F(b:e) = '-(...)'
!         WRITE(*,*)'4. F(b:e) = "-(...)"'
          CALL CompileSubstr (Comp, i, F, b+2, e-1, Var)
          CALL AddCompiledByte (Comp, i, cNeg)
          RETURN
       ELSEIF (SCAN(F(b+1:b+1),calpha) > 0) THEN
          n = MathFunctionIndex (F(b+1:e))
          IF (n > 0) THEN
             b2 = b+INDEX(F(b+1:e),'(')
             IF (CompletelyEnclosed(F, b2, e)) THEN          ! Case 5: F(b:e) = '-fcn(...)'
!               WRITE(*,*)'5. F(b:e) = "-fcn(...)"'
                ArgCnt=1                                         ! lam81
                k=b2+1                                           ! lam81
                pc=1                                             ! lam81
                DO WHILE (k < e)                                 ! lam81
                   IF (F(k:k) == '(') THEN                       ! lam81
                       pc = pc + 1                               ! lam81
                   ELSE IF (F(k:k) == ')') THEN                  ! lam81
                       pc = pc - 1                               ! lam81
                   ELSE IF (F(k:k) == ',' .AND. pc==1) THEN      ! lam81
                       CALL CompileSubstr(Comp, i, F, b2+1, k-1, Var)  ! lam81
                       b2=k                                      ! lam81
                       ArgCnt = ArgCnt + 1                       ! lam81
                   END IF                                        ! lam81
                   k=k+1                                         ! lam81
                END DO                                           ! lam81
                CALL CompileSubstr(Comp, i, F, b2+1, e-1, Var)
                CALL AddCompiledByte (Comp, i, n, ArgCnt)              ! lam81
                CALL AddCompiledByte (Comp, i, cNeg)
                RETURN
             END IF
          END IF
       ENDIF
    END IF
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    !----- -------- --------- --------- --------- --------- --------- --------- -------    
    DO io=cAdd,cPow                                          ! Increasing priority +-*/^
       k = 0
       DO j=e,b,-1
          IF     (F(j:j) == ')') THEN
             k = k+1
          ELSEIF (F(j:j) == '(') THEN
             k = k-1
          END IF
          IF (k == 0 .AND. F(j:j) == Ops(io) .AND. IsBinaryOp (j, F)) THEN
             IF (ANY(F(j:j) == Ops(cMul:cPow)) .AND. F(b:b) == '-') THEN ! Case 6: F(b:e) = '-...Op...' with Op > -
!               WRITE(*,*)'6. F(b:e) = "-...Op..." with Op > -'
                CALL CompileSubstr (Comp, i, F, b+1, e, Var)
                CALL AddCompiledByte (Comp, i, cNeg)
                RETURN                 
             ELSE                                                        ! Case 7: F(b:e) = '...BinOp...'
!               WRITE(*,*)'7. Binary operator',F(j:j)
                CALL CompileSubstr (Comp, i, F, b, j-1, Var)
                CALL CompileSubstr (Comp, i, F, j+1, e, Var)
                CALL AddCompiledByte (Comp, i, OperatorIndex(Ops(io)))
                Comp(i)%StackPtr = Comp(i)%StackPtr - 1
                RETURN
             END IF
          END IF
       END DO
    END DO
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for remaining items, i.e. variables or explicit numbers
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    b2 = b
    IF (F(b:b) == '-') b2 = b2+1
    n = MathItemIndex(Comp, i, F(b2:e), Var)
!   WRITE(*,*)'8. AddCompiledByte ',n
    CALL AddCompiledByte (Comp, i, n)
    Comp(i)%StackPtr = Comp(i)%StackPtr + 1
    IF (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
    IF (b2 > b) CALL AddCompiledByte (Comp, i, cNeg)
  END SUBROUTINE CompileSubstr
  !
  FUNCTION IsBinaryOp (j, F) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check if operator F(j:j) in string F is binary operator
    ! Special cases already covered elsewhere:              (that is corrected in v1.1)
    ! - operator character F(j:j) is first character of string (j=1)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,           INTENT(in) :: j                       ! Position of Operator
    CHARACTER (LEN=*), INTENT(in) :: F                       ! String
    LOGICAL                       :: res                     ! Result
    INTEGER                       :: k
    LOGICAL                       :: Dflag,Pflag
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    res=.true.
    IF (F(j:j) == '+' .OR. F(j:j) == '-') THEN               ! Plus or minus sign:
       IF (j == 1) THEN                                      ! - leading unary operator ?
          res = .false.
       ELSEIF (SCAN(F(j-1:j-1),'+-*/^(') > 0) THEN           ! - other unary operator ?
          res = .false.
       ELSEIF (SCAN(F(j+1:j+1),'0123456789') > 0 .AND. &     ! - in exponent of real number ?
               SCAN(F(j-1:j-1),'eEdD')       > 0) THEN
          Dflag=.false.; Pflag=.false.
          k = j-1
          DO WHILE (k > 1)                                   !   step to the left in mantissa 
             k = k-1
             IF     (SCAN(F(k:k),'0123456789') > 0) THEN
                Dflag=.true.
             ELSEIF (F(k:k) == '.') THEN
                IF (Pflag) THEN
                   EXIT                                      !   * EXIT: 2nd appearance of '.'
                ELSE
                   Pflag=.true.                              !   * mark 1st appearance of '.'
                ENDIF
             ELSE
                EXIT                                         !   * all other characters
             END IF
          END DO
          IF (Dflag .AND. (k == 1 .OR. SCAN(F(k:k),'+-*/^(') > 0)) res = .false.
       END IF
    END IF
  END FUNCTION IsBinaryOp
  !
  FUNCTION RealNum (str, ibegin, inext, error) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),  INTENT(in) :: str                    ! String
    REAL(dp)                       :: res                    ! Real number
    INTEGER, OPTIONAL, INTENT(out) :: ibegin,              & ! Start position of real number
                                      inext                  ! 1st character after real number
    LOGICAL, OPTIONAL, INTENT(out) :: error                  ! Error flag
    INTEGER                        :: ib,in,istat
    LOGICAL                        :: Bflag,               & ! .T. at begin of number in str
                                      InMan,               & ! .T. in mantissa of number
                                      Pflag,               & ! .T. after 1st '.' encountered
                                      Eflag,               & ! .T. at exponent identifier 'eEdD'
                                      InExp,               & ! .T. in exponent of number
                                      DInMan,              & ! .T. if at least 1 digit in mant.
                                      DInExp,              & ! .T. if at least 1 digit in exp.
                                      err                    ! Local error flag
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    Bflag=.true.; InMan=.false.; Pflag=.false.; Eflag=.false.; InExp=.false.
    DInMan=.false.; DInExp=.false.
    ib   = 1
    in   = 1
    DO WHILE (in <= LEN_TRIM(str))
       SELECT CASE (str(in:in))
       CASE (' ')                                            ! Only leading blanks permitted
          ib = ib+1
          IF (InMan .OR. Eflag .OR. InExp) EXIT
       CASE ('+','-')                                        ! Permitted only
          IF     (Bflag) THEN           
             InMan=.true.; Bflag=.false.                     ! - at beginning of mantissa
          ELSEIF (Eflag) THEN               
             InExp=.true.; Eflag=.false.                     ! - at beginning of exponent
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE ('0':'9')                                        ! Mark
          IF     (Bflag) THEN           
             InMan=.true.; Bflag=.false.                     ! - beginning of mantissa
          ELSEIF (Eflag) THEN               
             InExp=.true.; Eflag=.false.                     ! - beginning of exponent
          ENDIF
          IF (InMan) DInMan=.true.                           ! Mantissa contains digit
          IF (InExp) DInExp=.true.                           ! Exponent contains digit
       CASE ('.')
          IF     (Bflag) THEN
             Pflag=.true.                                    ! - mark 1st appearance of '.'
             InMan=.true.; Bflag=.false.                     !   mark beginning of mantissa
          ELSEIF (InMan .AND..NOT.Pflag) THEN
             Pflag=.true.                                    ! - mark 1st appearance of '.'
          ELSE
             EXIT                                            ! - otherwise STOP
          END IF
       CASE ('e','E','d','D')                                ! Permitted only
          IF (InMan) THEN
             Eflag=.true.; InMan=.false.                     ! - following mantissa
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE DEFAULT
          EXIT                                               ! STOP at all other characters
       END SELECT
       in = in+1
    END DO
    err = (ib > in-1) .OR. (.NOT.DInMan) .OR. ((Eflag.OR.InExp).AND..NOT.DInExp)
    IF (err) THEN
       res = zero
    ELSE
       READ(str(ib:in-1),*,IOSTAT=istat) res
       err = istat /= 0
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
    IF (PRESENT(error))  error  = err
  END FUNCTION RealNum
  !  
  SUBROUTINE LowCase (str1, str2)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Transform upper case letters in str1 into lower case letters, result is str2
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),  INTENT(in) :: str1
    CHARACTER (LEN=*), INTENT(out) :: str2
    INTEGER                        :: j,k
    CHARACTER (LEN=*),   PARAMETER :: lc = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    str2 = str1
    DO j=1,LEN_TRIM(str1)
       k = INDEX(uc,str1(j:j))
       IF (k > 0) str2(j:j) = lc(k:k)
    END DO
  END SUBROUTINE LowCase
  !
  ! function values and derivatives
  function dotf(fi,fj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3)
    real(dp)             :: res

    res = fi(1)*fj(1)+fi(2)*fj(2)+fi(3)*fj(3)

  end function dotf
  !
  function dotd(fi,fj,di,dj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),di(3),dj(3)
    real(dp)             :: res

    res = di(1)*fj(1)+fi(1)*dj(1)+di(2)*fj(2)+fi(2)*dj(2)+di(3)*fj(3)+fi(3)*dj(3)

  end function dotd
  !
  function crossxf(fi,fj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3)
    real(dp)             :: res

    res = fi(2)*fj(3)-fi(3)*fj(2)

  end function crossxf
  !
  function crossxd(fi,fj,di,dj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),di(3),dj(3)
    real(dp)             :: res

    res = di(2)*fj(3)+fi(2)*dj(3)-di(3)*fj(2)-fi(3)*dj(2)

  end function crossxd
  !
  function crossyf(fi,fj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3)
    real(dp)             :: res

    res = fi(3)*fj(1)-fi(1)*fj(3)
  
  end function crossyf
  !
  function crossyd(fi,fj,di,dj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),di(3),dj(3)
    real(dp)             :: res

    res = di(3)*fj(1)+fi(3)*dj(1)-di(1)*fj(3)-fi(1)*dj(3)

  end function crossyd
  !
  function crosszf(fi,fj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3)
    real(dp)             :: res

    res = fi(1)*fj(2)-fi(2)*fj(1)
  
  end function crosszf
  !
  function crosszd(fi,fj,di,dj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),di(3),dj(3)
    real(dp)             :: res

    res = di(1)*fj(2)+fi(1)*dj(2)-di(2)*fj(1)-fi(2)*dj(1)

  end function crosszd
  !
  function normf(fi) result(res)

    implicit none
    real(dp), intent(in) :: fi(3)
    real(dp)             :: res

    res = sqrt(fi(1)*fi(1)+fi(2)*fi(2)+fi(3)*fi(3))

  end function normf
  !
  function normd(fi,di) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),di(3)
    real(dp)             :: res
  
    res = (fi(1)*di(1)+fi(2)*di(2)+fi(3)*di(3))/normf(fi)

  end function normd
  !
  function disf(fi,fj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3)
    real(dp)             :: res

    res = normf(fi-fj)

  end function disf
  !
  function disd(fi,fj,di,dj) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),di(3),dj(3)
    real(dp)             :: res

    res = normd(fi-fj,di-dj)
  
  end function disd
  !
  function angf(fi,fj,fk) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),fk(3)
    real(dp)             :: res

    real(dp)             :: fij(3),fkj(3)

    fij=fi-fj
    fkj=fk-fj

    res = dotf(fij,fkj)/normf(fij)/normf(fkj) 

    if(res > one) then
       res = one
    else if(res < -one) then
       res = -one
    end if

    res = acos(res)

  end function angf
  !
  function angd(fi,fj,fk,di,dj,dk) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),fk(3),di(3),dj(3),dk(3)
    real(dp)             :: res

    real(dp)             :: fij(3),fkj(3),dij(3),dkj(3)
    real(dp)             :: normfij,normfkj,normdij,normdkj

    fij=fi-fj
    fkj=fk-fj
    dij=di-dj
    dkj=dk-dj

    normfij = normf(fij)
    normfkj = normf(fkj)
    normdij = normd(fij,dij)
    normdkj = normd(fkj,dkj)

    res = -(dotd(fij,fkj,dij,dkj)*normfij*normfkj-dotf(fij,fkj)*(normdij*normfkj+normfij*normdkj)) / &
           (sin(angf(fi,fj,fk))*normfij*normfij*normfkj*normfkj)

  end function angd
  !
  function dihf(fi,fj,fk,fl) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),fk(3),fl(3)
    real(dp)             :: res

    real(dp)             :: fij(3),fkj(3),fkl(3)
    real(dp)             :: fmj(3),fnk(3)

    fij=fi-fj
    fkj=fk-fj
    fkl=fk-fl

    fmj(1)=crossxf(fij,fkj)
    fmj(2)=crossyf(fij,fkj)
    fmj(3)=crosszf(fij,fkj)

    fnk(1)=crossxf(fkj,fkl)
    fnk(2)=crossyf(fkj,fkl)
    fnk(3)=crosszf(fkj,fkl)

    res = dotf(fmj,fnk)/normf(fmj)/normf(fnk)

    if(res > one) then
       res = one
    else if(res < -one) then
       res = -one
    end if

    res = sign(one,dotf(fij,fnk))*acos(res)
    
  end function dihf
  !
  function dihd(fi,fj,fk,fl,di,dj,dk,dl) result(res)

    implicit none
    real(dp), intent(in) :: fi(3),fj(3),fk(3),fl(3),di(3),dj(3),dk(3),dl(3)
    real(dp)             :: res

    real(dp)             :: fij(3),fkj(3),fkl(3)
    real(dp)             :: fmj(3),fnk(3)
    real(dp)             :: normfkj,normfkj2,normfmj2,normfnk2
    real(dp)             :: normfkj_normfmj2,normfkj_normfnk2
    real(dp)             :: fij_dot_fkj_normfkj2,fkj_dot_fkl_normfkj2

    fij=fi-fj
    fkj=fk-fj
    fkl=fk-fl

    fmj(1)=crossxf(fij,fkj)
    fmj(2)=crossyf(fij,fkj)
    fmj(3)=crosszf(fij,fkj)

    fnk(1)=crossxf(fkj,fkl)
    fnk(2)=crossyf(fkj,fkl)
    fnk(3)=crosszf(fkj,fkl)

    normfkj=normf(fkj)
    normfkj2=normfkj**2
    normfmj2=normf(fmj)**2
    normfnk2=normf(fnk)**2

    normfkj_normfmj2=normfkj/normfmj2
    normfkj_normfnk2=normfkj/normfnk2
    fij_dot_fkj_normfkj2=dotf(fij,fkj)/normfkj2
    fkj_dot_fkl_normfkj2=dotf(fkj,fkl)/normfkj2
    
    res = dotf(normfkj_normfmj2*fmj,di+(fij_dot_fkj_normfkj2-one)*dj-fij_dot_fkj_normfkj2*dk) - &
          dotf(normfkj_normfnk2*fnk,dl+(fkj_dot_fkl_normfkj2-one)*dk-fkj_dot_fkl_normfkj2*dj)

  end function dihd    
  !
END MODULE fparser_mod
