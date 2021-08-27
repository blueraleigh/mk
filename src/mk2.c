#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <phy.h>

#define twotothe256 \
    115792089237316195423570985008687907853269984665640564039457584007913129639936.0

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

#define SCLK(i, j) mk->h[(j) + (i) * 2]
#define DCLK(i, j) mk->g[(j) + (i) * 2]
#define UCLK(i, j) mk->f[(j) + (i) * 2]

struct mk {

    /* number of nodes */
    int n;

    /* rate of switching */
    double tau;

    /* asymmetry in rates */
    double eps;

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


static double pij(int i, int j, double tau, double eps, double len)
{
    double D = exp(-tau*len);
    double pi[2] = {
        1 / (1+eps), 
        eps / (1+eps)
    };
    return (i != j) ? (1-D) * pi[j] : pi[i] + D * pi[1-i];
}


static void branch_downpass(struct phy_node *node, struct mk *mk)
{
    int i;
    int j;
    int u = phy_node_index(node);
    double len = phy_node_brlen(node);

    for (i = 0; i < 2; ++i)
    {
        SCLK(u, i) = 0;
        for (j = 0; j < 2; ++j)
        {
            SCLK(u, i) += pij(i, j, mk->tau, mk->eps, len) * DCLK(u, j);
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

    for (i = 0; i < 2; ++i)
    {
        DCLK(u, i) = 1;
        for (d = phy_node_lfdesc(node); d != 0; d = phy_node_next(d))
        {
            v = phy_node_index(d);
            DCLK(u, i) *= SCLK(v, i);
            mk->lzd[u] += (i == 0) ? mk->lzd[v] : 0;
        }
    }

    for (i = 0; scale && (i < 2); ++i)
    {
        scale = (DCLK(u, i) < minlikelihood) && (DCLK(u, i) > minusminlikelihood);
    }

    if (scale)
    {
        mk->lzd[u] += 1;
        for (i = 0; i < 2; ++i)
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
        for (i = 0; i < 2; ++i)
            UCLK(v, i) = (i == 0) ? 1 / (1+mk->eps) : mk->eps / (1+mk->eps);
        return;
    }

    len = phy_node_brlen(node);
    u = phy_node_index(anc);

    for (i = 0; i < 2; ++i)
    {
        UCLK(v, i) = 0;
        for (j = 0; j < 2; ++j)
        {
            sclk = 1;
            for (sib = phy_node_next(node); sib != 0; sib = phy_node_next(sib))
            {
                w = phy_node_index(sib);
                sclk *= SCLK(w, j);
                mk->lzu[v] += (i == 0 && j == 0) ? mk->lzu[u] + mk->lzd[w] : 0;
            }
            for (sib = phy_node_prev(node); sib != 0; sib = phy_node_prev(sib))
            {
                w = phy_node_index(sib);
                sclk *= SCLK(w, j);
                mk->lzu[v] += (i == 0 && j == 0) ? mk->lzu[u] + mk->lzd[w] : 0;
            }
            UCLK(v, i) += UCLK(u, j) * sclk * pij(j, i, mk->tau, mk->eps, len);
        }
    }

    scale = 1;
    for (i = 0; scale && (i < 2); ++i)
        scale = (UCLK(v, i) < minlikelihood) && (UCLK(v, i) > minusminlikelihood);

    if (scale)
    {
        for (i = 0; i < 2; ++i)
            UCLK(v, i) *= twotothe256;
        mk->lzu[v] += 1;
    }
}


static double mk_loglikelihood(struct phy_node *node, struct mk *mk)
{
    int i;
    int u = phy_node_index(node);

    double lk = 0;

    for (i = 0; i < 2; ++i)
    {
        lk += (i == 0) ? DCLK(u, i) / (1+mk->eps) : DCLK(u, i) * (mk->eps / (1+mk->eps));
    }
    return log(lk) + mk->lzd[u] * log(minlikelihood);
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
    double sclk;

    struct phy_node *anc;
    struct phy_node *sib;

    anc = phy_node_anc(node);

    if (anc) {

        v = phy_node_index(node);
        u = phy_node_index(anc);

        for (i = 0; i < 2; ++i)
        {
            sclk = 1;
            for (sib = phy_node_next(node); sib != 0; sib = phy_node_next(sib))
            {
                w = phy_node_index(sib);
                sclk *= SCLK(w, i);
            }
            for (sib = phy_node_prev(node); sib != 0; sib = phy_node_prev(sib))
            {
                w = phy_node_index(sib);
                sclk *= SCLK(w, i);
            }
            lk += UCLK(u, i) * sclk * pij(i, j, mk->tau, mk->eps, len) * DCLK(v, j);
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
    for (i = 0; i < 2; ++i)
        asr[i] = asr_compute(i, node, mk);
    asr_normalize(2, asr);
}


void C_mk2_model_free(SEXP model)
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


SEXP C_mk2_model(SEXP x, SEXP rtree)
{
    struct mk *mk = calloc(1, sizeof(*mk));
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    mk->n = phy_nnode(phy);
    mk->phy = phy;

    mk->lzd = calloc(mk->n, sizeof(int));
    mk->lzu = calloc(mk->n, sizeof(int));
    mk->h = calloc(mk->n * 2, sizeof(double));
    mk->g = calloc(mk->n * 2, sizeof(double));
    mk->f = calloc(mk->n * 2, sizeof(double));
    memcpy(mk->g, REAL(x), 2 * phy_ntip(phy) * sizeof(double));

    SEXP exptr = PROTECT(R_MakeExternalPtr(mk, R_NilValue, R_NilValue));

    R_RegisterCFinalizer(exptr, C_mk2_model_free);

    UNPROTECT(1);
    return exptr;
}


SEXP C_mk2_loglikelihood(SEXP tau, SEXP eps, SEXP model)
{
    struct phy_node *node;
    struct phy_cursor *cursor;
    struct mk *mk = (struct mk *)R_ExternalPtrAddr(model);

    mk->tau = REAL(tau)[0];
    mk->eps = REAL(eps)[0];
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


SEXP C_mk2_marginal_asr(SEXP tau, SEXP eps, SEXP model)
{
    int i;
    struct phy_node *node;
    struct phy_cursor *cursor;
    struct mk *mk = (struct mk *)R_ExternalPtrAddr(model);

    mk->tau = REAL(tau)[0];
    mk->eps = REAL(eps)[0];
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

    SEXP ASR = PROTECT(allocMatrix(REALSXP, 2, mk->n));
    double *asr = REAL(ASR);

    cursor = phy_cursor_prepare(mk->phy, phy_root(mk->phy), ALL_NODES, 
        PREORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        node_uppass(node, mk);
        mk_asr_compute(asr + phy_node_index(node) * 2, node, mk);
    }

    UNPROTECT(1);
    return ASR;
}
