#include "BayesFM.h"
#include <stdlib.h>  // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

static R_NativePrimitiveArgType befa_t[] = {
    INTSXP,  INTSXP,  INTSXP,  INTSXP,  REALSXP, INTSXP,  LGLSXP,  INTSXP,
    REALSXP, LGLSXP,  INTSXP,  INTSXP,  INTSXP,  INTSXP,  LGLSXP,  REALSXP,
    INTSXP,  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, LGLSXP,  INTSXP,  REALSXP, INTSXP,  LGLSXP
    };

static R_NativePrimitiveArgType simnfacprior_t[] = {
    INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, LGLSXP
    };

static const R_FortranMethodDef ForMethods[] = {
    {"befa",         (DL_FUNC) &F77_NAME(befa),         31, befa_t        },
    {"simnfacprior", (DL_FUNC) &F77_NAME(simnfacprior),  8, simnfacprior_t},
    {NULL, NULL, 0}
    };

void attribute_visible R_init_BayesFM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, ForMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
