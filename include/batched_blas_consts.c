
#include "bblas_types.h"

#include <assert.h>
#include <stdbool.h>

bblas_enum_t bblas_diag_const(char lapack_char)
{
  if('N'==lapack_char){
    return BblasNonUnit;
  }
  if('n'==lapack_char){
    return BblasNonUnit;
  }
  if('U'==lapack_char){
    return BblasUnit;
  }
  if('u'==lapack_char){
    return BblasUnit;
  }
  return BblasInvalid;
}

bblas_enum_t bblas_info_const(char lapack_char)
{

  if('A'==lapack_char){
    return BblasErrorsReportAll;
  }
  if('a'==lapack_char){
    return BblasErrorsReportAll;
  }

  if('G'==lapack_char){
    return BblasErrorsReportGroup;
  }
  if('g'==lapack_char){
    return BblasErrorsReportGroup;
  }

  if('O'==lapack_char){
    return BblasErrorsReportAny;
  }
  if('o'==lapack_char){
    return BblasErrorsReportAny;
  }

  if('N'==lapack_char){
    return BblasErrorsReportNone;
  }
  if('n'==lapack_char){
    return BblasErrorsReportNone;
  }

  return BblasInvalid;
}

bblas_enum_t bblas_direct_const(char lapack_char)
{
  if('F'==lapack_char){
    return BblasForward;
  }
  if('f'==lapack_char){
    return BblasForward;
  }

  if('B'==lapack_char){
    return BblasBackward;
  }
  if('b'==lapack_char){
    return BblasBackward;
  }

  return BblasInvalid;

}
bblas_enum_t bblas_norm_const(char lapack_char)
{
  if('O'==lapack_char){
    return BblasOneNorm;
  }
  if('o'==lapack_char){
    return BblasOneNorm;
  }
  if('1'==lapack_char){
    return BblasOneNorm;
  }

  if('2'==lapack_char){
    return BblasTwoNorm;
  }
  if('F'==lapack_char){
    return BblasFrobeniusNorm;
  }
  if('f'==lapack_char){
    return BblasFrobeniusNorm;
  }
  if('E'==lapack_char){
    return BblasFrobeniusNorm;
  }
  if('e'==lapack_char){
    return BblasFrobeniusNorm;
  }

  if('I'==lapack_char){
    return BblasInfNorm;
  }
  if('i'==lapack_char){
    return BblasInfNorm;
  }

  if('M'==lapack_char){
    return BblasMaxNorm;
  }
  if('m'==lapack_char){
    return BblasMaxNorm;
  }

  return BblasInvalid;

}

bblas_enum_t bblas_side_const(char lapack_char)
{
  if('L'==lapack_char){
    return BblasLeft;
  }
  if('l'==lapack_char){
    return BblasLeft;
  }

  if('R'==lapack_char){
    return BblasRight;
  }
  if('r'==lapack_char){
    return BblasRight;
  }

  return BblasInvalid;
}

bblas_enum_t bblas_storev_const(char lapack_char)
{
  if('C'==lapack_char){
    return BblasColumnwise;
  }
  if('c'==lapack_char){
    return BblasColumnwise;
  }

  if('R'==lapack_char){
    return BblasRowwise;
  }
  if('r'==lapack_char){
    return BblasRowwise;
  }

  return BblasInvalid;
}

bblas_enum_t bblas_trans_const(char lapack_char)
{
  if('N'==lapack_char){
    return BblasNoTrans;
  }
  if('n'==lapack_char){
    return BblasNoTrans;
  }

  if('T'==lapack_char){
    return BblasTrans;
  }
  if('t'==lapack_char){
    return BblasTrans;
  }

  if('C'==lapack_char){
    return BblasConjTrans;
  }
  if('c'==lapack_char){
    return BblasConjTrans;
  }

  return BblasInvalid;
}

bblas_enum_t bblas_uplo_const(char lapack_char)
{

  if('U'==lapack_char){
    return BblasUpper;
  }
  if('u'==lapack_char){
    return BblasUpper;
  }

  if('L'==lapack_char){
    return BblasLower;
  }
  if('l'==lapack_char){
    return BblasLower;
  }

  return BblasGeneral;

}


