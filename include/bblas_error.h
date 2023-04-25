#ifndef ICL_BBLAS_ERROR_H
#define ICL_BBLAS_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include "bblas_types.h"
//#include "bblas_error.h"

#define bblas_error(msg) printf("BBLAS ERROR %s \n",msg);
static inline void bblas_set_info(int error_flag, int *info,
                                  int batch_count, int code) {
    int i;
    switch (error_flag) {
    case BblasErrorsReportAll :
        for( i=0; i < batch_count; i++) {
            info[i] = code;
        }
        break;
    case BblasErrorsReportGroup :
        info[0] = code;
        break;
    case BblasErrorsReportAny :
        info[0] = code;
        break;
    default :
        bblas_error("illegal value of info");
        info[0] = -1;
    }
}

#endif
