#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define PHY_API_IMPLEMENTATION
#include <phy.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


SEXP C_mk_shift(SEXP, SEXP);
SEXP C_mk_shift_backtrack(SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_mk_shift, 2),
    CALLDEF(C_mk_shift_backtrack, 2),
    {NULL, NULL, 0}
};


void attribute_visible R_init_mk(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
