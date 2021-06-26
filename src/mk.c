#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <phy.h>

#define twotothe256 \
    115792089237316195423570985008687907853269984665640564039457584007913129639936.0

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

#define SCLK(i, j) mk->h[(j) + (i) * mk->k]
#define DCLK(i, j) mk->g[(j) + (i) * mk->k]
#define UCLK(i, j) mk->f[(j) + (i) * mk->k]

struct mk {

    /* number of character states */
    int k;

    /* number of nodes */
    int n;

    /* rate of switching */
    double rate;

    /* downpass scaling factors */
    int *lzd;

    /* uppass scaling factors */
    int *lzu;

    /* downpass conditional likelihood arrays (at stem) */
    double *h;

    /* downpass conditional likelihood arrays (at node) */
    double *g;

    /* uppass conditional likelihood arrays (at node) */
    double *f;

    struct phy *phy;
};


static double pij(int i, int j, int k, double rate, double len)
{
    double D = exp(-rate*k*len);
    return (i != j) ? (1 - D) / (double)k : D + (1 - D) / (double)k;
}


static void branch_downpass(struct phy_node *node, struct mk *mk)
{
    int i;
    int j;
    int u = phy_node_index(node);
    double len = phy_node_brlen(node);

    for (i = 0; i < mk->k; ++i)
    {
        SCLK(u, i) = 0;
        for (j = 0; j < mk->k; ++j)
        {
            SCLK(u, i) += pij(i, j, mk->k, mk->rate, len) * DCLK(u, j);
        }
    }
}


static void node_downpass(struct phy_node *node, struct mk *mk)
{
    int i;
    int v;
    int u = phy_node_index(node);

    int scale = 1;

    struct phy_node *d;

    for (i = 0; i < mk->k; ++i)
    {
        DCLK(u, i) = 1;
        for (d = phy_node_lfdesc(node); d != 0; d = phy_node_next(d))
        {
            v = phy_node_index(d);
            DCLK(u, i) *= SCLK(v, i);
            mk->lzd[u] += (i == 0) ? mk->lzd[v] : 0;
        }
    }

    for (i = 0; scale && (i < mk->k); ++i)
    {
        scale = (DCLK(u, i) < minlikelihood) && (DCLK(u, i) > minusminlikelihood);
    }

    if (scale)
    {
        mk->lzd[u] += 1;
        for (i = 0; i < mk->k; ++i)
        {
            DCLK(u, i) *= twotothe256;
        }
    }
}


static void node_uppass(struct phy_node *node, struct mk *mk)
{
    int i;
    int j;
    int u;
    int v;
    int w;
    int scale = 1;

    double len;
    double sclk;

    struct phy_node *sib;
    struct phy_node *anc = phy_node_anc(node);

    v = phy_node_index(node);

    if (!anc)
    {
        for (i = 0; i < mk->k; ++i)
            UCLK(v, i) = 1 / (double)(mk->k);
        return;
    }

    len = phy_node_brlen(node);
    u = phy_node_index(anc);

    for (i = 0; i < mk->k; ++i)
    {
        UCLK(v, i) = 0;
        for (j = 0; j < mk->k; ++j)
        {
            sclk = 1;
            for (sib = phy_node_next(node); sib != 0; sib = phy_node_next(sib))
            {
                w = phy_node_index(sib);
                sclk *= SCLK(w, j);
                mk->lzu[v] += (i == 0) ? mk->lzu[u] + mk->lzd[w] : 0;
            }
            for (sib = phy_node_prev(node); sib != 0; sib = phy_node_prev(sib))
            {
                w = phy_node_index(sib);
                sclk *= SCLK(w, j);
                mk->lzu[v] += (i == 0) ? mk->lzu[u] + mk->lzd[w] : 0;
            }
            UCLK(v, i) += UCLK(u, j) * sclk * pij(j, i, mk->k, mk->rate, len);
        }
    }

    scale = 1;
    for (i = 0; scale && (i < mk->k); ++i)
        scale = (UCLK(v, i) < minlikelihood) && (UCLK(v, i) > minusminlikelihood);

    if (scale)
    {
        for (i = 0; i < mk->k; ++i)
            UCLK(v, i) *= twotothe256;
        mk->lzu[v] += 1;
    }
}


static double mk_loglikelihood(struct phy_node *node, struct mk *mk)
{
    int i;
    int u = phy_node_index(node);

    double lk = 0;

    for (i = 0; i < mk->k; ++i)
    {
        lk += DCLK(u, i);
    }
    return log(lk) - log(mk->k) + mk->lzd[u] * log(minlikelihood);
}


static double asr_compute(
    int j,
    struct phy_node *node,
    struct mk *mk)
{
    int i;
    int u;
    int v;
    int w;
    double len = phy_node_brlen(node);
    double lk = 0;

    struct phy_node *anc;
    struct phy_node *sib;

    anc = phy_node_anc(node);

    if (anc) {

        v = phy_node_index(node);
        u = phy_node_index(anc);

        for (i = 0; i < mk->k; ++i)
        {
            for (sib = phy_node_next(node); sib != 0; sib = phy_node_next(sib))
            {
                w = phy_node_index(sib);
                
                lk += UCLK(u, i) * SCLK(w, i) * pij(i, j, mk->k, mk->rate, len) * DCLK(v, j);
            }
            for (sib = phy_node_prev(node); sib != 0; sib = phy_node_prev(sib))
            {
                w = phy_node_index(sib);
                
                lk += UCLK(u, i) * SCLK(w, i) * pij(i, j, mk->k, mk->rate, len) * DCLK(v, j);
            }
        }

    } else {

        return DCLK(phy_node_index(node), j);

    }

    return lk;
}


static void asr_normalize(int r, double *out)
{
    int i;
    double norm = 0;

    for (i = 0; i < r; ++i)
        norm += out[i];

    for (i = 0; i < r; ++i)
        out[i] /= norm;
}


static void mk_asr_compute(double *asr, struct phy_node *node, struct mk *mk)
{
    int i;
    for (i = 0; i < mk->k; ++i)
        asr[i] = asr_compute(i, node, mk);
    asr_normalize(mk->k, asr);
}


void C_mk_model_free(SEXP model)
{
    struct mk *mk = (struct mk *)R_ExternalPtrAddr(model);
    free(mk->lzd);
    free(mk->lzu);
    free(mk->h);
    free(mk->g);
    free(mk->f);
    free(mk);
    R_ClearExternalPtr(model);
}


SEXP C_mk_model(SEXP x, SEXP rtree)
{
    struct mk *mk = calloc(1, sizeof(*mk));
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    mk->k = INTEGER(getAttrib(x, R_DimSymbol))[0];
    mk->n = phy_nnode(phy);
    mk->phy = phy;

    mk->lzd = calloc(mk->n, sizeof(int));
    mk->lzu = calloc(mk->n, sizeof(int));
    mk->h = calloc(mk->n * mk->k, sizeof(double));
    mk->g = calloc(mk->n * mk->k, sizeof(double));
    mk->f = calloc(mk->n * mk->k, sizeof(double));
    memcpy(mk->g, REAL(x), mk->k * phy_ntip(phy) * sizeof(double));

    SEXP exptr = PROTECT(R_MakeExternalPtr(mk, R_NilValue, R_NilValue));

    R_RegisterCFinalizer(exptr, C_mk_model_free);

    UNPROTECT(1);
    return exptr;
}


SEXP C_mk_loglikelihood(SEXP rate, SEXP model)
{
    struct phy_node *node;
    struct phy_cursor *cursor;
    struct mk *mk = (struct mk *)R_ExternalPtrAddr(model);

    mk->rate = REAL(rate)[0];
    cursor = phy_cursor_prepare(mk->phy, phy_root(mk->phy), ALL_NODES, 
        POSTORDER);

    memset(mk->lzd, 0, mk->n * sizeof(int));

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        if (phy_node_istip(node))
        {
            branch_downpass(node, mk);
        }
        else
        {
            node_downpass(node, mk);
            branch_downpass(node, mk);
        }
    }
    return ScalarReal(mk_loglikelihood(phy_root(mk->phy), mk));
}


SEXP C_mk_marginal_asr(SEXP rate, SEXP model)
{
    int i;
    struct phy_node *node;
    struct phy_cursor *cursor;
    struct mk *mk = (struct mk *)R_ExternalPtrAddr(model);

    mk->rate = REAL(rate)[0];
    cursor = phy_cursor_prepare(mk->phy, phy_root(mk->phy), ALL_NODES, 
        POSTORDER);


    memset(mk->lzd, 0, mk->n * sizeof(int));
    memset(mk->lzu, 0, mk->n * sizeof(int));

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        if (phy_node_istip(node))
        {
            branch_downpass(node, mk);
        }
        else
        {
            node_downpass(node, mk);
            branch_downpass(node, mk);
        }
    }

    SEXP ASR = PROTECT(allocMatrix(REALSXP, mk->k, mk->n));
    double *asr = REAL(ASR);

    cursor = phy_cursor_prepare(mk->phy, phy_root(mk->phy), ALL_NODES, 
        PREORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        node_uppass(node, mk);
        mk_asr_compute(asr + phy_node_index(node) * mk->k, node, mk);
    }

    UNPROTECT(1);
    return ASR;
}
