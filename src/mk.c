#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <phy.h>

#define twotothe256 \
    115792089237316195423570985008687907853269984665640564039457584007913129639936.0

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

#define SCLK(i, j) mk->h[(i) + (j) * mk->n]
#define DCLK(i, j) mk->g[(i) + (j) * mk->n]

struct mk {

    /* number of character states */
    int k;

    /* number of nodes */
    int n;

    /* rate of switching */
    double rate;

    /* downpass scaling factors */
    int *lzd;

    /* downpass conditional likelihood arrays (at stem) */
    double *h;

    /* downpass conditional likelihood arrays (at node) */
    double *g;

    void *mem;

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
        for (d = phy_node_lfdesc(node); d; d = phy_node_next(d))
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


static double mk_loglikelihood(struct phy_node *node, struct mk *mk)
{
    int i;
    int u = phy_node_index(node);

    double lk = 0;

    for (i = 0; i < mk->k; ++i)
    {
        lk += DCLK(u, i);
        Rprintf("%g %d\n", DCLK(u, i), mk->lzd[u]);
    }
    return log(lk) - log(mk->k) + mk->lzd[u] * log(minlikelihood);
}


void C_mk_model_free(SEXP model)
{
    struct mk *mk = (struct mk *)R_ExternalPtrAddr(model);
    free(mk->mem);
    free(mk);
    R_ClearExternalPtr(model);
}


SEXP C_mk_model(SEXP x, SEXP rtree)
{
    struct mk *mk = calloc(1, sizeof(*mk));
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    mk->k = INTEGER(getAttrib(x, R_DimSymbol))[1];
    mk->n = phy_nnode(phy);
    mk->phy = phy;

    mk->mem = malloc(mk->n * sizeof(int) + 2 * mk->k * mk->n * sizeof(double));

    mk->lzd = (int *)(mk->mem);
    mk->h = (double *)(mk->lzd + mk->n);
    mk->g = mk->h + mk->k * mk->n;

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
