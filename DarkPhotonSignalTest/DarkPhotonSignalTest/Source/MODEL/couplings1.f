ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'

      INCLUDE 'model_functions.inc'
      GC_1 = -(MDL_EEM*MDL_COMPLEXI)
      GC_3 = MDL_COMPLEXI*MDL_GEDP
      END
