/* Upper bounds on F(n) for n beyond exact-enumeration range.
 *
 * Starts from the exact optimal signing at n = 10 (from exact_fn) and
 * grows one vertex at a time.  At each level every extension of every
 * beam member is evaluated exactly (max over all 2^(n-1) sign vectors,
 * Gray-code walk), and the best `BEAM` extensions are kept.
 *
 * Output: certified upper bounds F(n) <= best_M for n = 11..MAXN,
 * with the witness signings printed as edge masks for reverification.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAXN 16
#define BEAM 6

static int8_t Acur[BEAM][MAXN][MAXN];
static int8_t Anext[BEAM][MAXN][MAXN];

/* exact max |S| via Gray walk; A is k x k, zero diagonal */
static long long max_abs(const int8_t A[MAXN][MAXN], int k)
{
    static int x[MAXN];
    static long long t[MAXN];
    for (int i = 0; i < k; i++) x[i] = 1;
    for (int i = 0; i < k; i++) {
        t[i] = 0;
        for (int j = 0; j < k; j++) if (j != i) t[i] += A[i][j];
    }
    long long S = 0;
    for (int i = 0; i < k; i++) S += t[i];
    S /= 2;
    long long best = llabs(S);
    uint64_t steps = 1ULL << (k - 1);
    for (uint64_t c = 1; c < steps; c++) {
        int f = __builtin_ctzll(c) + 1;
        S -= 2LL * x[f] * t[f];
        int xf = x[f]; x[f] = -xf;
        for (int j = 0; j < k; j++) if (j != f) t[j] -= 2LL * xf * A[j][f];
        if (llabs(S) > best) best = llabs(S);
    }
    return best;
}

int main(void)
{
    /* exact F(10) witness: vertex 0 all-plus, reduced edge mask below */
    uint64_t mask = 229890161ULL;
    int n = 10;
    memset(Acur, 0, sizeof Acur);
    for (int j = 1; j < n; j++) Acur[0][0][j] = Acur[0][j][0] = 1;
    int e = 0;
    for (int i = 1; i < n; i++)
        for (int j = i + 1; j < n; j++, e++) {
            int8_t s = (mask >> e) & 1 ? -1 : 1;
            Acur[0][i][j] = Acur[0][j][i] = s;
        }
    long long m0 = max_abs(Acur[0], 10);
    printf("n=10 start_max=%lld (expect 13)\n", m0);
    if (m0 != 13) { fprintf(stderr, "bad start witness\n"); return 2; }
    int nbeam = 1;

    for (n = 11; n <= MAXN; n++) {
        long long bestM[BEAM];
        uint64_t bestExt[BEAM];
        int bestSrc[BEAM];
        for (int b = 0; b < BEAM; b++) bestM[b] = 1LL << 60;

        uint64_t nex = 1ULL << (n - 1);
        for (int src = 0; src < nbeam; src++) {
            for (uint64_t ext = 0; ext < nex; ext++) {
                static int8_t T[MAXN][MAXN];
                memcpy(T, Acur[src], sizeof T);
                for (int j = 0; j < n - 1; j++) {
                    int8_t s = (ext >> j) & 1 ? -1 : 1;
                    T[j][n - 1] = T[n - 1][j] = s;
                }
                long long m = max_abs(T, n);
                /* insert into beam (keep distinct values simple: allow dups) */
                if (m < bestM[BEAM - 1]) {
                    int p = BEAM - 1;
                    while (p > 0 && bestM[p - 1] > m) {
                        bestM[p] = bestM[p - 1]; bestExt[p] = bestExt[p - 1];
                        bestSrc[p] = bestSrc[p - 1]; p--;
                    }
                    bestM[p] = m; bestExt[p] = ext; bestSrc[p] = src;
                }
            }
        }
        int newbeam = nbeam * 2 > BEAM ? BEAM : nbeam * 2;
        if (newbeam > BEAM) newbeam = BEAM;
        for (int b = 0; b < newbeam; b++) {
            memcpy(Anext[b], Acur[bestSrc[b]], sizeof Anext[b]);
            for (int j = 0; j < n - 1; j++) {
                int8_t s = (bestExt[b] >> j) & 1 ? -1 : 1;
                Anext[b][j][n - 1] = Anext[b][n - 1][j] = s;
            }
        }
        memcpy(Acur, Anext, sizeof Acur);
        nbeam = newbeam;
        printf("n=%d upper_bound=%lld ratio=%.6f ext_mask=%llu\n",
               n, bestM[0], (double)bestM[0] / pow(n, 1.5),
               (unsigned long long)bestExt[0]);
        fflush(stdout);
    }
    /* print the final best signing for reverification */
    printf("final_matrix_rows:\n");
    for (int i = 0; i < MAXN; i++) {
        for (int j = 0; j < MAXN; j++) putchar(Acur[0][i][j] > 0 ? '+' : (i == j ? '0' : '-'));
        putchar('\n');
    }
    return 0;
}
