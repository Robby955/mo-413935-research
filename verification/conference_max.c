/* Exact Boolean maxima of symmetric conference matrices and their
 * principal submatrices.
 *
 * Constructions:
 *   - Paley conference of order q+1 for prime q = 1 mod 4 (q = 5,13,17,29)
 *   - Petersen two-graph Seidel matrix for order 10 (J - I - 2A)
 *   - GF(25) Paley conference for order 26
 * Every matrix is checked to satisfy C^2 = (n-1) I before use.
 *
 * max_x |sum_{i<j} C_ij x_i x_j| is computed exactly by a Gray-code walk
 * over the 2^(n-1) sign vectors with x_0 = +1, maintaining t = Cx.
 *
 * Usage: conference_max mode
 *   mode = maxima  : exact maxima for orders 6,10,14,18,26,30
 *   mode = chains  : restriction-decay chains inside the order-30 matrix
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static int C[32][32];

static int legendre(int a, int p)
{
    a %= p; if (a < 0) a += p;
    if (a == 0) return 0;
    int r = 1, base = a, e = (p - 1) / 2;
    while (e) { if (e & 1) r = r * base % p; base = base * base % p; e >>= 1; }
    return r == 1 ? 1 : -1;
}

static void build_paley(int q)
{
    int n = q + 1;
    memset(C, 0, sizeof C);
    for (int j = 1; j < n; j++) C[0][j] = C[j][0] = 1;
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            if (i != j) C[i][j] = legendre(i - j, q);
}

/* GF(25) = GF(5)[t]/(t^2-2); chi(z) = z^12 */
typedef struct { int a, b; } f25;
static f25 mul25(f25 x, f25 y)
{
    f25 r;
    r.a = (x.a * y.a + 2 * x.b * y.b) % 5;
    r.b = (x.a * y.b + x.b * y.a) % 5;
    return r;
}
static int chi25(f25 z)
{
    if (z.a == 0 && z.b == 0) return 0;
    f25 r = {1, 0}, base = z;
    int e = 12;
    while (e) { if (e & 1) r = mul25(r, base); base = mul25(base, base); e >>= 1; }
    if (r.b != 0 || (r.a != 1 && r.a != 4)) { fprintf(stderr, "chi25 broken\n"); exit(2); }
    return r.a == 1 ? 1 : -1;
}
static void build_gf25(void)
{
    int n = 26;
    memset(C, 0, sizeof C);
    for (int j = 1; j < n; j++) C[0][j] = C[j][0] = 1;
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            if (i != j) {
                f25 x = {(i - 1) % 5, (i - 1) / 5}, y = {(j - 1) % 5, (j - 1) / 5};
                f25 d = {(x.a - y.a + 5) % 5, (x.b - y.b + 5) % 5};
                C[i][j] = chi25(d);
            }
}

static void build_petersen(void)
{
    int pairs[10][2], k = 0;
    for (int a = 0; a < 5; a++)
        for (int b = a + 1; b < 5; b++) { pairs[k][0] = a; pairs[k][1] = b; k++; }
    memset(C, 0, sizeof C);
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++) {
            if (i == j) continue;
            int adj = (pairs[i][0] != pairs[j][0] && pairs[i][0] != pairs[j][1] &&
                       pairs[i][1] != pairs[j][0] && pairs[i][1] != pairs[j][1]);
            C[i][j] = adj ? -1 : 1;   /* J - I - 2A */
        }
}

static void check_conference(int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            long long s = 0;
            for (int k2 = 0; k2 < n; k2++) s += (long long)C[i][k2] * C[k2][j];
            long long want = i == j ? n - 1 : 0;
            if (s != want) { fprintf(stderr, "not conference at %d %d\n", i, j); exit(2); }
        }
}

static long long g_maxS, g_minS;   /* signed spread of the last boolean_max call */
static uint64_t g_bestx;           /* witness sign vector (bit i = vertex i negative) */

/* exact max |S(x)| over x with x_0=+1, for the principal submatrix on idx[0..k-1] */
static long long boolean_max(const int *idx, int k)
{
    static int x[32];
    static long long t[32];
    for (int i = 0; i < k; i++) x[i] = 1;
    for (int i = 0; i < k; i++) {
        t[i] = 0;
        for (int j = 0; j < k; j++) if (j != i) t[i] += C[idx[i]][idx[j]];
    }
    long long S = 0;
    for (int i = 0; i < k; i++) S += t[i];
    S /= 2;
    long long best = S < 0 ? -S : S;
    g_maxS = S; g_minS = S;
    g_bestx = 0;
    uint64_t steps = 1ULL << (k - 1);
    for (uint64_t c = 1; c < steps; c++) {
        int f = __builtin_ctzll(c) + 1;      /* flip vertex f (never 0) */
        S -= 2LL * x[f] * t[f];
        int xf = x[f];
        x[f] = -xf;
        for (int j = 0; j < k; j++)
            if (j != f) t[j] -= 2LL * xf * C[idx[j]][idx[f]];
        long long a = S < 0 ? -S : S;
        if (a > best) {
            best = a;
            g_bestx = 0;
            for (int i = 1; i < k; i++) if (x[i] < 0) g_bestx |= 1ULL << i;
        }
        if (S > g_maxS) g_maxS = S;
        if (S < g_minS) g_minS = S;
        if ((c & 0xFFFFFF) == 0x54321) {     /* corruption control */
            long long q = 0;
            for (int i = 0; i < k; i++)
                for (int j = i + 1; j < k; j++)
                    q += (long long)C[idx[i]][idx[j]] * x[i] * x[j];
            if (q != S) { fprintf(stderr, "CORRUPTION\n"); exit(2); }
        }
    }
    return best;
}

static void report(const char *name, int n)
{
    check_conference(n);
    int idx[32];
    for (int i = 0; i < n; i++) idx[i] = i;
    long long m = boolean_max(idx, n);
    double ratio = (double)m / (pow(n, 1.5));
    double ceil2 = 0.5 * sqrt((double)(n - 1) / n);
    printf("order=%d source=%s maximum=%lld maxS=%lld minS=%lld "
           "ratio=%.9f spectral_ceiling=%.9f witness_x=%llu\n",
           n, name, m, g_maxS, g_minS, ratio, ceil2,
           (unsigned long long)g_bestx);
}

static uint64_t rng_state = 413935;
static uint64_t rng(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}

int main(int argc, char **argv)
{
    const char *mode = argc > 1 ? argv[1] : "maxima";
    if (!strcmp(mode, "maxima")) {
        build_paley(5);  report("paley_q5", 6);
        build_petersen(); report("petersen", 10);
        build_paley(13); report("paley_q13", 14);
        build_paley(17); report("paley_q17", 18);
        build_gf25();    report("gf25", 26);
        build_paley(29); report("paley_q29", 30);
    } else if (!strcmp(mode, "chains")) {
        build_paley(29);
        check_conference(30);
        for (int chain = 0; chain < 3; chain++) {
            int idx[32];
            for (int i = 0; i < 30; i++) idx[i] = i;
            for (int i = 29; i > 0; i--) {   /* Fisher-Yates */
                int j = rng() % (i + 1);
                int tmp = idx[i]; idx[i] = idx[j]; idx[j] = tmp;
            }
            printf("chain=%d\n", chain);
            for (int k = 5; k <= 30; k++) {
                long long m = boolean_max(idx, k);
                printf("  k=%d max=%lld ratio=%.6f\n", k, m, (double)m / pow(k, 1.5));
            }
        }
    }
    return 0;
}
