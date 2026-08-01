/* Exact computation of
 *   F(n)  = min over sign matrices of max_x |sum_{i<j} a_ij x_i x_j|
 *   F+(n) = min over sign matrices of max_x  sum_{i<j} a_ij x_i x_j
 * by enumeration of switching classes.
 *
 * Every switching class has a unique representative in which all edges
 * at vertex 0 are +1.  The remaining C(n-1,2) edge signs are enumerated
 * in Gray-code order, updating all 2^(n-1) energies incrementally
 * (x_0 = +1 fixed, since the energy is even in x).
 *
 * Usage: exact_fn n [fixbits fixval]
 *   fixbits/fixval: enumerate only classes whose top `fixbits` edge bits
 *   equal fixval (for splitting across processes).  Default 0 0.
 *
 * Output lines are machine-readable key=value.
 * Corruption controls: periodic from-scratch recomputation of one energy,
 * plus a final from-scratch audit of the reported minima.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n, R, ER, K;
static int8_t **sgn;         /* sgn[e][x] = x_i x_j for reduced edge e */
static int ei[64], ej[64];   /* endpoints (1-based vertices) of edge e */

static long long energy_from_scratch(uint64_t mask, int x)
{
    /* all edges at vertex 0 are +1; bit b of x = sign of vertex b+1 (1 -> -1) */
    long long q = 0;
    int xi, xj, i, j, e = 0;
    for (i = 1; i < n; i++) {
        xi = (x >> (i - 1)) & 1 ? -1 : 1;
        q += xi;                       /* edge (0,i), sign +1, x_0 = +1 */
    }
    for (i = 1; i < n; i++)
        for (j = i + 1; j < n; j++, e++) {
            xi = (x >> (i - 1)) & 1 ? -1 : 1;
            xj = (x >> (j - 1)) & 1 ? -1 : 1;
            q += ((mask >> e) & 1 ? -1 : 1) * xi * xj;
        }
    return q;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s n [fixbits fixval]\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    int fixbits = argc > 3 ? atoi(argv[2]) : 0;
    uint64_t fixval = argc > 3 ? strtoull(argv[3], 0, 10) : 0;
    R = n - 1; ER = R * (R - 1) / 2; K = 1 << R;
    if (ER > 62) { fprintf(stderr, "n too large\n"); return 1; }

    int e = 0;
    for (int i = 1; i < n; i++)
        for (int j = i + 1; j < n; j++, e++) { ei[e] = i; ej[e] = j; }

    sgn = malloc(ER * sizeof(int8_t *));
    for (e = 0; e < ER; e++) {
        sgn[e] = malloc(K);
        for (int x = 0; x < K; x++) {
            int b = ((x >> (ei[e] - 1)) ^ (x >> (ej[e] - 1))) & 1;
            sgn[e][x] = b ? -1 : 1;
        }
    }

    int8_t *Q = malloc(K);
    /* initial class: low bits all zero, top fixbits = fixval */
    uint64_t mask0 = fixval << (ER - fixbits);
    for (int x = 0; x < K; x++) {
        int s = 1;                      /* x_0 */
        for (int i = 1; i < n; i++) s += (x >> (i - 1)) & 1 ? -1 : 1;
        long long q = (long long)(s * s - n) / 2;   /* all-plus energy */
        for (e = 0; e < ER; e++)
            if ((mask0 >> e) & 1) q -= 2 * sgn[e][x];
        Q[x] = (int8_t)q;
    }

    int best_abs = 127, best_plus = 127, best_osc = 1000;
    uint64_t best_abs_mask = mask0, best_plus_mask = mask0, best_osc_mask = mask0;
    long long n_abs_ties = 0, n_plus_ties = 0, n_osc_ties = 0;

    uint64_t nlow = 1ULL << (ER - fixbits);
    uint64_t cur = mask0;

    for (uint64_t c = 0; ; c++) {
        if (c) {   /* Gray-code flip */
            int fe = __builtin_ctzll(c);
            uint64_t bit = 1ULL << fe;
            int8_t d = (cur & bit) ? 2 : -2;   /* leaving -1 -> +2, entering -1 -> -2 */
            cur ^= bit;
            int8_t *s = sgn[fe];
            for (int x = 0; x < K; x++) Q[x] += d * s[x];
        }
        int ma = 0, mp = -127, mn = 127;
        for (int x = 0; x < K; x++) {
            int q = Q[x];
            int a = q < 0 ? -q : q;
            if (a > ma) ma = a;
            if (q > mp) mp = q;
            if (q < mn) mn = q;
        }
        if (ma < best_abs)  { best_abs = ma;  best_abs_mask = cur;  n_abs_ties = 1; }
        else if (ma == best_abs) n_abs_ties++;
        if (mp < best_plus) { best_plus = mp; best_plus_mask = cur; n_plus_ties = 1; }
        else if (mp == best_plus) n_plus_ties++;
        if (mp - mn < best_osc) { best_osc = mp - mn; best_osc_mask = cur; n_osc_ties = 1; }
        else if (mp - mn == best_osc) n_osc_ties++;

        if ((c & 0xFFFFF) == 0x12345) {   /* corruption control: spot recompute */
            int x = (int)(c % K);
            if (Q[x] != energy_from_scratch(cur, x)) {
                fprintf(stderr, "CORRUPTION at c=%llu\n", (unsigned long long)c);
                return 2;
            }
        }
        if (c + 1 == nlow) break;
    }

    /* final audit: recompute the reported minima from scratch */
    int audit_abs = 0, audit_plus = -1000000;
    for (int x = 0; x < K; x++) {
        long long q = energy_from_scratch(best_abs_mask, x);
        long long a = q < 0 ? -q : q;
        if (a > audit_abs) audit_abs = (int)a;
    }
    for (int x = 0; x < K; x++) {
        long long q = energy_from_scratch(best_plus_mask, x);
        if (q > audit_plus) audit_plus = (int)q;
    }
    if (audit_abs != best_abs || audit_plus != best_plus) {
        fprintf(stderr, "AUDIT FAILED abs %d vs %d, plus %d vs %d\n",
                audit_abs, best_abs, audit_plus, best_plus);
        return 3;
    }

    printf("n=%d fixbits=%d fixval=%llu classes=%llu\n", n, fixbits,
           (unsigned long long)fixval, (unsigned long long)nlow);
    printf("F_abs=%d witness_mask=%llu ties=%lld\n",
           best_abs, (unsigned long long)best_abs_mask, n_abs_ties);
    printf("F_plus=%d witness_mask=%llu ties=%lld\n",
           best_plus, (unsigned long long)best_plus_mask, n_plus_ties);
    printf("Osc=%d witness_mask=%llu ties=%lld\n",
           best_osc, (unsigned long long)best_osc_mask, n_osc_ties);
    printf("audit=PASSED\n");
    return 0;
}
