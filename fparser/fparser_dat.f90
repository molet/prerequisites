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

module fparser_dat_mod

  use prec_mod
  use constants_mod

IMPLICIT NONE

  INTEGER(is),                              PARAMETER :: cImmed   = 1,          &
                                                         ! operators
                                                         cNeg     = 2,          &
                                                         cAdd     = 3,          &
                                                         cSub     = 4,          &
                                                         cMul     = 5,          &
                                                         cDiv     = 6,          &
                                                         cPow     = 7,          &
                                                         ! functions with single argument
                                                         cAbs     = 8,          &
                                                         cExp     = 9,          &
                                                         cLog10   = 10,         &
                                                         cLog     = 11,         &
                                                         cSqrt    = 12,         &
                                                         cSinh    = 13,         &
                                                         cCosh    = 14,         &
                                                         cTanh    = 15,         &
                                                         cSin     = 16,         &
                                                         cCos     = 17,         &
                                                         cTan     = 18,         &
                                                         cAsinh   = 19,         &
                                                         cAcosh   = 20,         &
                                                         cAtanh   = 21,         &
                                                         cAsin    = 22,         &
                                                         cAcos    = 23,         &
                                                         cAtan    = 24,         &
                                                         cErfc    = 25,         &
                                                         cErf     = 26,         &
                                                         cSgn     = 27,         &
                                                         ! functions with two arguments
                                                         cSign    = 28,         &
                                                         cLogb    = 29,         &
                                                         ! functions with arbitrary number of arguments
                                                         cSum     = 30,         &
                                                         cAve     = 31,         &
                                                         cMean    = 32,         &
                                                         cAmean   = 33,         &
                                                         cProd    = 34,         &
                                                         cGmean   = 35,         &
                                                         cHmean   = 36,         &
                                                         cMin     = 37,         &
                                                         cMax     = 38,         &
                                                         ! spatial functions
                                                         cWrap    = 39,         &
                                                         cDot     = 40,         &
                                                         cCrossx  = 41,         &
                                                         cCrossy  = 42,         &
                                                         cCrossz  = 43,         &
                                                         cNorm    = 44,         &
                                                         cDis     = 45,         &
                                                         cAng     = 46,         &
                                                         cDih     = 47,         &
                                                         VarBegin = 48
  CHARACTER (LEN=1), DIMENSION(cAdd:cPow),  PARAMETER :: Ops      = (/ '+',     &
                                                                       '-',     &
                                                                       '*',     &
                                                                       '/',     &
                                                                       '^' /)
  CHARACTER (LEN=6), DIMENSION(cPow+1:VarBegin-1), PARAMETER :: Funcs  = (/ 'abs   ', &
                                                                            'exp   ', &
                                                                            'log10 ', &
                                                                            'log   ', &
                                                                            'sqrt  ', &
                                                                            'sinh  ', &
                                                                            'cosh  ', &
                                                                            'tanh  ', &
                                                                            'sin   ', &
                                                                            'cos   ', &
                                                                            'tan   ', &
                                                                            'asinh ', &
                                                                            'acosh ', &
                                                                            'atanh ', &
                                                                            'asin  ', &
                                                                            'acos  ', &
                                                                            'atan  ', &
                                                                            'erfc  ', &
                                                                            'erf   ', &
                                                                            'sgn   ', &
                                                                            'sign  ', &
                                                                            'logb  ', &
                                                                            'sum   ', &
                                                                            'ave   ', &
                                                                            'mean  ', &
                                                                            'amean ', &
                                                                            'prod  ', &
                                                                            'gmean ', &
                                                                            'hmean ', &
                                                                            'min   ', &
                                                                            'max   ', &
                                                                            'wrap  ', &
                                                                            'dot   ', &
                                                                            'crossx', &
                                                                            'crossy', &
                                                                            'crossz', &
                                                                            'norm  ', &
                                                                            'dis   ', &
                                                                            'ang   ', &
                                                                            'dih   ' /)

  INTEGER :: INDX
  INTEGER, DIMENSION(cPow+1:VarBegin-1), PARAMETER :: NArgFuncs = (/ (1, INDX=1,cSgn-cAbs+1), &
                                                                     (2, INDX=1,cLogb-cSign+1), &
                                                                     (0, INDX=1,cMax-cSum+1), &
                                                                     2,6,6,6,6,3,6,9,12 /)

  CHARACTER (LEN=35), DIMENSION(17), PARAMETER :: EvalfErrText = (/ 'Division by zero                   ', &
                                                                    'Argument of LOG10 <= 0             ', &
                                                                    'Argument of LOG <= 0               ', &
                                                                    'Argument of SQRT < 0               ', &
                                                                    'Argument of ACOSH < 1              ', &
                                                                    'Argument of ATANH < -1             ', &
                                                                    'Argument of ATANH > 1              ', &
                                                                    'Argument of ASIN < -1              ', &
                                                                    'Argument of ASIN > 1               ', &
                                                                    'Argument of ACOS < -1              ', &
                                                                    'Argument of ACOS > 1               ', &
                                                                    'Base of LOGB <= 0                  ', &
                                                                    'Base of LOGB = 1                   ', &
                                                                    'Argument of LOGB <= 0              ', &
                                                                    'Product of arguments of GMEAN < 0  ', &
                                                                    'At least one argument of HMEAN = 0 ', &
                                                                    'Denominator in HMEAN = 0           ' /)

  CHARACTER (LEN=35), DIMENSION(14), PARAMETER :: EvaldErrText = (/ 'ABS is not differentiable at 0     ', &
                                                                    'SQRT is not differentiable at 0    ', &
                                                                    'ACOSH is not differentiable at -1  ', &
                                                                    'ACOSH is not differentiable at 1   ', &
                                                                    'ATANH is not differentiable at -1  ', &
                                                                    'ATANH is not differentiable at 1   ', &
                                                                    'ASIN is not differentiable at -1   ', &
                                                                    'ASIN is not differentiable at 1    ', &
                                                                    'ACOS is not differentiable at -1   ', &
                                                                    'ACOS is not differentiable at 1    ', &
                                                                    'NORM is not differentiable at 0    ', &
                                                                    'DIS is not differentiable at 0     ', &
                                                                    'ANG is not differentiable at -pi   ', &
                                                                    'ANG is not differentiable at pi    ' /)
 


  TYPE tComp
     INTEGER(is), DIMENSION(:), POINTER :: ByteCode => null()
     INTEGER                            :: ByteCodeSize
     REAL(dp),    DIMENSION(:), POINTER :: Immed => null()
     INTEGER                            :: ImmedSize
     INTEGER,     DIMENSION(:), POINTER :: NArg => null()   ! lam81
     INTEGER                            :: StackSize
     INTEGER                            :: StackPtr
     INTEGER,     DIMENSION(:), POINTER :: Var => null()
  END TYPE tComp
! TYPE (tComp),  DIMENSION(:),  POINTER :: Comp              ! Bytecode
  INTEGER,   DIMENSION(:),  ALLOCATABLE :: ipos              ! Associates function strings

end module fparser_dat_mod
