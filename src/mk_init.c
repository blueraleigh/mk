#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define PHY_API_IMPLEMENTATION
#include <phy.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


SEXP C_mk_model(SEXP, SEXP);
SEXP C_mk_loglikelihood(SEXP, SEXP);
SEXP C_mk_marginal_asr(SEXP, SEXP);
SEXP C_mk2_model(SEXP, SEXP);
SEXP C_mk2_loglikelihood(SEXP, SEXP);
SEXP C_mk2_marginal_asr(SEXP, SEXP);
SEXP C_mk_shift(SEXP, SEXP);
SEXP C_mk_shift_backtrack(SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_mk_model, 2),
    CALLDEF(C_mk_loglikelihood, 2),
    CALLDEF(C_mk_marginal_asr, 2),
    CALLDEF(C_mk2_model, 2),
    CALLDEF(C_mk2_loglikelihood, 3),
    CALLDEF(C_mk2_marginal_asr, 3),
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
