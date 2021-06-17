#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <phy.h>


/*
** Given a character with k discrete states and set X of terminal node
** character states related by phylogeny T, we assume that the character
** has evolved on T to produce X under a fully symmetric Markov
** process of character state evolution with rate r. Assuming a rate of
** evolution slow enough that the expected number of character state
** changes on any given branch of T is typically much less than 1, the
** log probability of X is approximately
**
** log P(x) ~ N * log(r) - r * V + const
**
** where the constant includes terms that do not depend on r. Here, V is the
** sum of all branch lengths in T and N is the minimum number of character 
** state changes needed to explain X, which is readily computed via Fitch 
** parsimony.
** 
** The maximum likelihood estimate of r is then simply N / V.
**
** This struct tracks these two quantities (N and V) needed to compute
** the likelihood.
*/
struct mk {
    // minimum number of character state changes needed
    // to explain data, computed via Fitch's downpass algorithm
    int n;
    // summed branch length over which changes occur.
    // ML estimate of rate is then just n / v.
    double v;
};


struct dp {
    /* number of shifted processes */
    int n_proc;

    /* only filled entries are valid */
    int filled;

    /* Fitch parsimony downpass state set */
    int x;

    /* log likelihood of shifted processes */
    double logp;

    /* log likelihood of root process and shifted processes */
    double score;

    /* sufficient statistics log for the root process */
    struct mk mk;

    /* index of item in the dp array of left descendant that was
    ** used to form this entry */
    int j;

    /* index of item in the dp array of right descendant that was
    ** used to form this entry */
    int k;
};


static void mk_push(int n, double v, struct mk *mk)
{
    mk->n += n;
    mk->v += v;
}


static void mk_merge(struct mk *a, struct mk *b, struct mk *c)
{
    c->n = a->n + b->n;
    c->v = a->v + b->v;
}


static double mk_loglikelihood(struct mk *mk)
{
    double rate = mk->n / mk->v;
    return rate > 0 ? mk->n * log(rate) - rate * mk->v : 0;
}


static void downpass(struct phy *phy, int *n_edge)
{
    int i;
    int j;
    int k;
    int x;

    double score;

    struct phy_node *lf;
    struct phy_node *rt;
    struct phy_node *node;
    struct phy_cursor *cursor;

    struct dp *vec;
    struct dp *l_vec;
    struct dp *r_vec;

    struct mk mk;
    struct mk empty = {0, 0};

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        lf = phy_node_lfdesc(node);
        rt = phy_node_rtdesc(node);

        vec = (struct dp *)phy_node_data(node);
        l_vec = (struct dp *)phy_node_data(lf);
        r_vec = (struct dp *)phy_node_data(rt);

        mk_merge(&(l_vec[0].mk), &(r_vec[0].mk), &(vec[0].mk));

        x = l_vec[0].x & r_vec[0].x;
        if (!x)
        {
            x = l_vec[0].x | r_vec[0].x;
            mk_push(1, phy_node_brlen(lf) + phy_node_brlen(rt), &(vec[0].mk));
        }
        else
            mk_push(0, phy_node_brlen(lf) + phy_node_brlen(rt), &(vec[0].mk));
        vec[0].n_proc = 1;
        vec[0].x = x;
        vec[0].logp = mk_loglikelihood(&(vec[0].mk));
        vec[0].score = vec[0].logp;
        vec[0].filled = 1;
        vec[0].j = -1;
        vec[0].k = -1;

        for (j = 0; j <= n_edge[phy_node_index(lf)]; ++j)
        {
            for (k = 0; k <= n_edge[phy_node_index(rt)]; ++k)
            {
                if (!l_vec[j].filled || !r_vec[k].filled)
                    continue;

                i = j + k + 2;

                if (j == 0 && k == 0)
                    mk_merge(&empty, &empty, &mk);
                else if (j == 0 && k > 0)
                    mk_merge(&empty, &(r_vec[k].mk), &mk);
                else if (j > 0 && k == 0)
                    mk_merge(&(l_vec[j].mk), &empty, &mk);
                else
                    mk_merge(&(l_vec[j].mk), &(r_vec[k].mk), &mk);

                x = l_vec[j].x & r_vec[k].x;
                if (!x)
                {
                    x = l_vec[j].x | r_vec[k].x;
                    mk_push(1, phy_node_brlen(lf) + phy_node_brlen(rt), &mk);
                }
                else
                    mk_push(0, phy_node_brlen(lf) + phy_node_brlen(rt), &mk);

                score = l_vec[j].logp + r_vec[k].logp + mk_loglikelihood(&mk);

                if (!vec[i].filled || score > vec[i].score)
                {
                    vec[i].n_proc = l_vec[j].n_proc + r_vec[k].n_proc;
                    vec[i].filled = 1;
                    vec[i].score = score;
                    vec[i].j = j;
                    vec[i].k = k;
                    vec[i].logp = l_vec[j].logp + r_vec[k].logp;
                    vec[i].x = x;
                    vec[i].mk = mk;
                }
            }
        }
    }
}


static void backtrack(int index, struct phy_node *node, struct phy *phy,
    double aic_w, double bg_rate, double *rate, int *i, int *shift)
{
    struct dp *vec = (struct dp *)phy_node_data(node);

    if (vec[index].j != -1 && vec[index].k != -1)
    {
        rate[phy_node_index(node)] += aic_w * bg_rate;
        backtrack(vec[index].j, phy_node_lfdesc(node), phy, aic_w,
            bg_rate, rate, i, shift);
        backtrack(vec[index].k, phy_node_rtdesc(node), phy, aic_w,
            bg_rate, rate, i, shift);
    }
    else
    {
        if (phy_node_istip(node))
        {
            rate[phy_node_index(node)] += aic_w * bg_rate;
        }
        else
        {
            struct phy_cursor *cursor;

            if (i && shift && node != phy_root(phy))
                shift[(*i)++] = phy_node_index(node)+1;

            // the shifted rate applies on all branches descended from this
            // node, not on the branch subtending this node (which receives
            // the background rate)
            rate[phy_node_index(node)] += aic_w * bg_rate;
            cursor = phy_cursor_prepare(phy, node, ALL_NODES, PREORDER);
            phy_cursor_step(cursor);
            while ((node = phy_cursor_step(cursor)) != 0)
            {
                rate[phy_node_index(node)] += aic_w * (vec[index].mk.n / vec[index].mk.v);
            }
        }
    }
}


static void dp_init(int *x, int *n_edge, struct phy *phy)
{
    int i;
    struct phy_node *node;
    struct phy_cursor *cursor;
    struct dp *dp;

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        if (phy_node_istip(node))
        {
            n_edge[phy_node_index(node)] = 0;//1;
        }
        else
        {
            n_edge[phy_node_index(node)] =
                n_edge[phy_node_index(phy_node_lfdesc(node))] +
                n_edge[phy_node_index(phy_node_rtdesc(node))] + 2;
        }
    }

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        dp = calloc(n_edge[phy_node_index(node)] + 1, sizeof(*dp));
        phy_node_add_data(node, (void *)dp, &free);
    }

    for (i = 0; i < phy_ntip(phy); ++i)
    {
        node = phy_node_get(phy, i);
        dp = (struct dp *)phy_node_data(node);
        dp[0].x = 1<<(x[i]-1);
        dp[0].n_proc = 1;
        dp[0].filled = 1;
        dp[0].j = -1;
        dp[0].k = -1;
    }
}


static double aic(double lnL, int k)
{
    return -2*lnL + 2*k;
}


static double bic(double lnL, int k, int n)
{
    return -2*lnL + log(n)*k;
}


SEXP C_mk_shift(SEXP x, SEXP rtree)
{
    int i;
    int n;
    int best = 0;
    double w = 0;
    double score = R_NegInf;
    struct dp *dp;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int n_edge[phy_nnode(phy)];
    double rate[phy_nnode(phy)];

    dp_init(INTEGER(x), n_edge, phy);
    downpass(phy, n_edge);
    dp = (struct dp *)phy_node_data(phy_root(phy));

    double aic_score[n_edge[phy_node_index(phy_root(phy))]+1];
    double aic_w[n_edge[phy_node_index(phy_root(phy))]+1];
    double min_aic_score = aic_score[0] = aic(dp[0].score, 2*dp[0].n_proc);

    for (i = 1; i <= n_edge[phy_node_index(phy_root(phy))]; ++i)
    {
        if (dp[i].filled)
        {
            aic_score[i] = aic(dp[i].score, 2*dp[i].n_proc);
            if (aic_score[i] < min_aic_score)
            {
                min_aic_score = aic_score[i];
                best = i;
            }
        }
    }

    for (i = 0; i <= n_edge[phy_node_index(phy_root(phy))]; ++i)
    {
        aic_w[i] = 0;
        if (dp[i].filled)
            w += aic_w[i] = exp(-0.5 * (aic_score[i] - min_aic_score));
    }

    memset(rate, 0, phy_nnode(phy) * sizeof(double));

    for (i = 0; i <= n_edge[phy_node_index(phy_root(phy))]; ++i)
    {
        if (dp[i].filled)
        {
            aic_w[i] /= w;
            backtrack(i, phy_root(phy), phy, aic_w[i], dp[i].mk.n / dp[i].mk.v,
                rate, 0, 0);
        }
    }

    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, phy_nnode(phy)));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, phy_nnode(phy)));

    memcpy(REAL(VECTOR_ELT(ans, 0)), rate, phy_nnode(phy) * sizeof(double));
    memcpy(REAL(VECTOR_ELT(ans, 1)), aic_w, phy_nnode(phy) * sizeof(double));

    UNPROTECT(1);
    return ans;
}


SEXP C_mk_shift_backtrack(SEXP index, SEXP rtree)
{
    int i = *INTEGER(index);
    int j = 0;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    struct dp *dp = (struct dp *)phy_node_data(phy_root(phy));
    int shift[phy_nnode(phy)];
    SEXP rate = PROTECT(allocVector(REALSXP, phy_nnode(phy)));
    memset(REAL(rate), 0, phy_nnode(phy) * sizeof(double));
    if (dp && dp[i].filled) {
        backtrack(i, phy_root(phy), phy, 1, dp[i].mk.n / dp[i].mk.v,
            REAL(rate), &j, shift);
    }

    SEXP rates = PROTECT(allocVector(REALSXP, j+1));
    SEXP shifts = PROTECT(allocVector(INTSXP, j+1));
    INTEGER(shifts)[0] = phy_node_index(phy_root(phy)) + 1;
    REAL(rates)[0] = REAL(rate)[phy_node_index(phy_root(phy))];
    for (i = 1; i < j+1; ++i)
    {
        REAL(rates)[i] = REAL(rate)[shift[i-1]-1];
        INTEGER(shifts)[i] = shift[i-1];
    }

    setAttrib(rates, R_NamesSymbol, shifts);
    setAttrib(rate, install("rate"), rates);

    UNPROTECT(3);
    return rate;
}

