/* jglue.c -- exhaustive census of induced Omega-join ("gluing (J)") instances.
 *
 * Input: gentourng output for one order n (upper-triangle bit strings, '1' at
 * (i,j), i<j, means i->j).  For every tournament T and every ordered pair of
 * vertex subsets (A,B) with A not inside B and B not inside A, test
 *
 *     H(T[A u B]) == H(T[A]) + H(T[B]) - 1          (cover form, U = A u B)
 *     H(T)        == H(T[A]) + H(T[B]) - 1          (aux form,  A u B < V)
 *
 * where H is the Hamiltonian path count.  By THM-4114 / inclusion-exclusion,
 * H(A u B) = H(A) + H(B) - H(W) + M, W = A n B, M = mixed pack weight >= 0,
 * so (J) holds iff M = H(W) - 1.  Cover instances are recorded at the order
 * |A u B| at which they occur (we only test U = V, smaller U are covered by
 * the census of smaller orders).
 *
 * Also checked:
 *   (N2) one-vertex overlap, prime p,q: H(A u B) >= p+q+3, and equality only
 *        if A-v and B-v are both transitive (one-apex at the shared vertex);
 *   (N3) W transitive and A-W={x}, B-W={y}: never H = p+q-1 with p,q >= 3;
 *   CONDITIONAL.md quantitative bound on every cover.
 *
 * usage: gentourng n | ./jglue n [aux]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 10
static int n, N;
static uint32_t out_[MAXN];
static uint32_t H[1 << MAXN];
static uint32_t dp[1 << MAXN][MAXN];
static char isprime[1 << 16];

static void sieve(void) {
    int L = 1 << 16;
    memset(isprime, 1, L);
    isprime[0] = isprime[1] = 0;
    for (int i = 2; i * i < L; i++)
        if (isprime[i])
            for (int j = i * i; j < L; j += i) isprime[j] = 0;
}

static void allH(void) {
    H[0] = 1;
    for (int S = 1; S < N; S++) {
        uint32_t tot = 0;
        for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) continue;
            int R = S & ~(1 << v);
            uint32_t c;
            if (!R) c = 1;
            else {
                c = 0;
                for (int u = 0; u < n; u++)
                    if ((R >> u & 1) && (out_[u] >> v & 1)) c += dp[R][u];
            }
            dp[S][v] = c;
            tot += c;
        }
        H[S] = tot;
    }
}

/* statistics */
static long long nT = 0;
static long long cov_all = 0, cov_pp = 0, cov_pp_exact = 0;  /* unordered */
static long long cov_byw[MAXN + 1], cov_pp_byw[MAXN + 1];
static long long aux_pp = 0;
static long long n2_cases = 0, n2_eq = 0, n2_eq_oneapex = 0, n2_viol = 0;
static int n2_min = 1 << 30, n2b_min = 1 << 30;
static long long n3_cases = 0, n3_hits = 0;
static long long cond_viol = 0, cond_tested = 0, cond2_viol = 0;
/* distinct (p,q,N) with p<=q primes >= 11, cover form: min overlap seen */
#define VMAX 4096
static unsigned char minw_pq[VMAX][VMAX / 8]; /* lazily: store min overlap+1 in a map instead */
typedef struct { int p, q, Nv, w, r, exact, ord; uint32_t A, B; char t[64]; } inst_t;
static inst_t *inst = NULL; static int ninst = 0, capinst = 0;

static int find_inst(int p, int q) {
    for (int i = 0; i < ninst; i++) if (inst[i].p == p && inst[i].q == q) return i;
    return -1;
}

static char cur_line[256];

static void record(int p, int q, int Nv, int w, int r, int exact, uint32_t A, uint32_t B) {
    if (p > q) { int t = p; p = q; q = t; uint32_t T = A; A = B; B = T; }
    int i = find_inst(p, q);
    if (i < 0) {
        if (ninst == capinst) { capinst = capinst ? 2 * capinst : 256; inst = realloc(inst, capinst * sizeof(inst_t)); }
        i = ninst++;
        inst[i].p = p; inst[i].q = q; inst[i].Nv = Nv; inst[i].w = 99;
    }
    if (w < inst[i].w) {
        inst[i].w = w; inst[i].r = r; inst[i].exact = exact; inst[i].A = A; inst[i].B = B;
        inst[i].ord = n;
        strncpy(inst[i].t, cur_line, 63); inst[i].t[63] = 0;
    }
}

static int popc(uint32_t x) { return __builtin_popcount(x); }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: jglue n [aux]\n"); return 1; }
    n = atoi(argv[1]); N = 1 << n;
    int doaux = argc > 2;
    sieve();
    char line[256];
    int V = N - 1;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' ) continue;
        int len = strlen(line);
        while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
        if (len != n * (n - 1) / 2) continue;
        strcpy(cur_line, line);
        memset(out_, 0, sizeof out_);
        int k = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++, k++) {
                if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i;
            }
        allH();
        nT++;
        uint32_t HV = H[V];
        /* cover form: A u B = V, both proper, neither contains the other. Enumerate
           ordered pairs with A < B numerically to make them unordered. */
        for (uint32_t A = 1; A < (uint32_t)V; A++) {
            uint32_t hA = H[A];
            uint32_t compA = V & ~A;
            /* B = compA | X, X subset of A, X != A (B proper) ; need A not subset of B: X != A */
            for (uint32_t X = A;; X = (X - 1) & A) {
                uint32_t B = compA | X;
                if (X != A && B > A && (A & ~B)) {
                    uint32_t hB = H[B];
                    uint32_t W = A & B;
                    uint32_t hW = H[W];
                    int w = popc(W);
                    /* CONDITIONAL quantitative bound */
                    {
                        long long a = H[A & ~W], b = H[B & ~W], r = hW, p = hA, q = hB;
                        long long bound = p + q - r + (a - 1) * (q - r) + (b - 1) * (p - r) - (a - 1) * (b - 1);
                        cond_tested++;
                        if ((long long)HV < bound) { cond_viol++; if (getenv("DBG2") && cond_viol <= 3) printf("CV A=%u B=%u p=%lld q=%lld r=%lld a=%lld b=%lld H=%u bound=%lld corrected=%lld T=%s\n", A, B, p, q, r, a, b, HV, bound, p+q-r+(a-1)*(q-r)+(b-1)*(p-r)-r*(a-1)*(b-1), cur_line); }
                        { long long cb = p+q-r+(a-1)*(q-r)+(b-1)*(p-r)-r*(a-1)*(b-1); if ((long long)HV < cb) cond2_viol++; }
                    }
                    if (hA >= 3 && hB >= 3) {
                        int pp = isprime[hA] && isprime[hB];
                        if (w == 1 && pp && !(hA == 3 && hB == 3)) { int d2 = (int)HV - (int)(hA + hB); if (d2 < n2b_min) n2b_min = d2; if (d2 == 7 && getenv("DBG4")) printf("N2B d=%d p=%u q=%u A=%u B=%u T=%s\n", d2, hA, hB, A, B, cur_line); }
                        if (w == 1 && pp) {
                            n2_cases++;
                            int d = (int)HV - (int)(hA + hB);
                            if (d < n2_min) n2_min = d;
                            if (d < 3) n2_viol++;
                            if (d == 3) {
                                n2_eq++;
                                if (H[A & ~W] == 1 && H[B & ~W] == 1) { n2_eq_oneapex++; if (getenv("DBG3")) printf("N2OA %u %u\n", hA < hB ? hA : hB, hA < hB ? hB : hA); }
                                else if (getenv("DBG")) printf("N2EQ p=%u q=%u A=%u B=%u H(A-v)=%u H(B-v)=%u T=%s\n", hA, hB, A, B, H[A&~W], H[B&~W], cur_line);
                            }
                        }
                        if (popc(A & ~W) == 1 && popc(B & ~W) == 1 && hW == 1) {
                            n3_cases++;
                            if (HV == hA + hB - 1) { n3_hits++; if (getenv("DBG")) printf("N3J p=%u q=%u N=%u A=%u B=%u W=%u T=%s\n", hA, hB, HV, A, B, W, cur_line); }
                        }
                        if (HV == hA + hB - 1) {
                            cov_all++; cov_byw[w]++;
                            if (pp && hA >= 11 && hB >= 11) {
                                cov_pp++; cov_pp_byw[w]++;
                                long long M = (long long)HV - hA - hB + hW;
                                int exact = (M == 0);
                                if (exact) cov_pp_exact++;
                                record(hA, hB, HV, w, hW, exact, A, B);
                            }
                        }
                    }
                }
                if (X == 0) break;
            }
        }
        if (doaux) {
            /* aux form: A u B strictly inside V; restrict to prime counts >= 11 */
            static uint32_t L[1 << MAXN]; int nL = 0;
            for (uint32_t S = 1; S < (uint32_t)V; S++) if (H[S] >= 11 && isprime[H[S]]) L[nL++] = S;
            for (int i = 0; i < nL; i++)
                for (int j = i + 1; j < nL; j++) {
                    uint32_t A = L[i], B = L[j];
                    if ((A | B) == (uint32_t)V) continue;
                    if (!(A & ~B) || !(B & ~A)) continue;
                    if (HV == H[A] + H[B] - 1) aux_pp++;
                }
        }
    }
    printf("n=%d tournaments=%lld\n", n, nT);
    printf("cover (J) instances, unordered, H(A),H(B)>=3: %lld   by overlap w:", cov_all);
    for (int w = 0; w <= n; w++) printf(" %d:%lld", w, cov_byw[w]);
    printf("\n");
    printf("cover (J) with both prime >= 11: %lld (exact M=0: %lld)   by w:", cov_pp, cov_pp_exact);
    for (int w = 0; w <= n; w++) printf(" %d:%lld", w, cov_pp_byw[w]);
    printf("\n");
    if (doaux) printf("aux (J) with A u B < V, both prime >= 11: %lld\n", aux_pp);
    printf("N2 one-vertex overlap, prime p,q>=3: cases=%lld min(H-p-q)=%d violations(<3)=%lld equality=%lld equality&one-apex=%lld\n",
           n2_cases, n2_cases ? n2_min : 0, n2_viol, n2_eq, n2_eq_oneapex);
    printf("N2b one-vertex overlap, primes not both 3: min(H-p-q)=%d\n", n2b_min);
    printf("N3 two apexes over a common transitive W: cases=%lld exact joins=%lld\n", n3_cases, n3_hits);
    printf("CONDITIONAL bound: tested=%lld violations=%lld ; corrected bound (intersection weight r(a-1)(b-1)) violations=%lld\n", cond_tested, cond_viol, cond2_viol);
    for (int i = 0; i < ninst; i++)
        printf("PQ p=%d q=%d N=%d minw=%d H(W)=%d exact=%d A=%u B=%u T=%s\n", inst[i].p, inst[i].q, inst[i].Nv,
               inst[i].w, inst[i].r, inst[i].exact, inst[i].A, inst[i].B, inst[i].t);
    return 0;
}
