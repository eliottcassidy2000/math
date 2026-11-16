/* glue.c <nA> <nB> [strongonly]: additive gluings on a shared arc.  For every class A of order nA and B of order nB
   (gentourng), every arc u->v of A and u'->v' of B (identified), every orientation of the cross arcs between A-{u,v} and
   B-{u,v}: H(U), H(A), H(B), the order statistics h_A(uv), h_A(vu), h_B(uv), h_B(vu), the path-splice bound and whether it is
   attained.  Reports counts of H(U) - H(A) - H(B) in {-3,-1,1,3} with bound-attained flags, and prints examples of the
   -1 and +1 cases with p, q >= 3 (A, B as gentourng strings, arcs, cross mask). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static i64 dp[1<<10][10][2];
/* count Hamiltonian paths of the sub-tournament on mask S with vertices adjacency adj (n total), returning total and the
   number with u before v */
static void count(u64 *adj, int n, int S, int u, int v, i64 *tot, i64 *uv) {
    int N = 1 << n; for (int m = 0; m < N; m++) if ((m & S) == m) for (int l = 0; l < n; l++) dp[m][l][0] = dp[m][l][1] = 0;
    for (int x = 0; x < n; x++) if ((S >> x) & 1) dp[1<<x][x][0] = 1;
    for (int m = 1; m < N; m++) { if ((m & S) != m) continue; for (int l = 0; l < n; l++) for (int f = 0; f < 2; f++) { i64 a = dp[m][l][f]; if (!a) continue; u64 nx = adj[l] & S & ~(u64)m; while (nx) { int x = __builtin_ctzll(nx); nx &= nx-1; int f2 = f; if (x == v && ((m >> u) & 1)) f2 = 1; dp[m|(1<<x)][x][f2] += a; } } }
    *tot = 0; *uv = 0; for (int l = 0; l < n; l++) { *tot += dp[S][l][0] + dp[S][l][1]; *uv += dp[S][l][1]; }
}
static void parse(const char *s, int n, u64 *adj) { int k = 0; for (int i = 0; i < n; i++) adj[i] = 0; for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (s[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; } }
static u64 reachfrom(u64 s, u64 *A, int n) { u64 r = s, p = 0; while (r != p) { p = r; for (int x = 0; x < n; x++) if ((r>>x)&1) r |= A[x]; } return r; }
static int strongT(u64 *adj, int n) { u64 radj[12]; for (int x = 0; x < n; x++) { radj[x] = 0; for (int y = 0; y < n; y++) if ((adj[y]>>x)&1) radj[x] |= 1ULL<<y; } u64 full = (1ULL<<n)-1; return reachfrom(1, adj, n) == full && reachfrom(1, radj, n) == full; }
int main(int argc, char **argv) {
    int nA = atoi(argv[1]), nB = atoi(argv[2]); int strongonly = argc > 3;
    char cmd[64]; static char As[20000][64], Bs[20000][64]; int cA = 0, cB = 0;
    sprintf(cmd, "gentourng -q %d", nA); FILE *f = popen(cmd, "r"); while (fscanf(f, "%63s", As[cA]) == 1) cA++; pclose(f);
    sprintf(cmd, "gentourng -q %d", nB); f = popen(cmd, "r"); while (fscanf(f, "%63s", Bs[cB]) == 1) cB++; pclose(f);
    long total = 0, hist[7][2] = {{0}}; long shown = 0;
    for (int ia = 0; ia < cA; ia++) for (int ib = 0; ib < cB; ib++) {
        u64 A[12], B[12]; parse(As[ia], nA, A); parse(Bs[ib], nB, B);
        if (strongonly && (!strongT(A, nA) || !strongT(B, nB))) continue;
        i64 HA, HAuv, HB, HBuv; 
        for (int u = 0; u < nA; u++) for (int v = 0; v < nA; v++) if ((A[u]>>v)&1) for (int u2 = 0; u2 < nB; u2++) for (int v2 = 0; v2 < nB; v2++) if ((B[u2]>>v2)&1) {
            count(A, nA, (1<<nA)-1, u, v, &HA, &HAuv); count(B, nB, (1<<nB)-1, u2, v2, &HB, &HBuv);
            /* build U: A's vertices 0..nA-1; B's others mapped to nA.. */
            int n = nA + nB - 2; int mapB[12]; int idx = nA; for (int x = 0; x < nB; x++) mapB[x] = (x == u2) ? u : (x == v2) ? v : idx++;
            u64 base[12]; for (int i = 0; i < n; i++) base[i] = 0;
            for (int i = 0; i < nA; i++) base[i] = A[i];
            for (int x = 0; x < nB; x++) for (int y = 0; y < nB; y++) if ((B[x]>>y)&1) { if ((x == u2 || x == v2) && (y == u2 || y == v2)) continue; base[mapB[x]] |= 1ULL << mapB[y]; }
            int crossA[12], ca = 0, crossB[12], cb = 0; for (int i = 0; i < nA; i++) if (i != u && i != v) crossA[ca++] = i; for (int i = nA; i < n; i++) crossB[cb++] = i;
            int np = ca * cb;
            for (int bits = 0; bits < (1 << np); bits++) {
                u64 adj[12]; for (int i = 0; i < n; i++) adj[i] = base[i]; int t = 0;
                for (int i = 0; i < ca; i++) for (int j = 0; j < cb; j++) { if ((bits >> t) & 1) adj[crossA[i]] |= 1ULL << crossB[j]; else adj[crossB[j]] |= 1ULL << crossA[i]; t++; }
                i64 HU, HUuv; count(adj, n, (1<<n)-1, u, v, &HU, &HUuv); total++;
                i64 diff = HU - HA - HB; if (diff < -3 || diff > 3 || (diff % 2 == 0)) continue;
                i64 bound = HAuv * HBuv + (HA - HAuv) * (HB - HBuv); int att = (bound == HU);
                hist[diff + 3][att]++;
                if ((diff == -1 || diff == 1) && HA >= 3 && HB >= 3 && shown < 24) { shown++; printf("diff=%lld H(U)=%lld H(A)=%lld H(B)=%lld bound=%lld attained=%d A=%s arc=(%d,%d) B=%s arc=(%d,%d) cross=%d hA(uv)=%lld hB(uv)=%lld strongU=%d\n", diff, HU, HA, HB, bound, att, As[ia], u, v, Bs[ib], u2, v2, bits, HAuv, HBuv, strongT(adj, n)); }
            }
        }
    }
    printf("orders (%d,%d)%s: instances %ld; by diff: -3: %ld+%ld (bound attained / not), -1: %ld+%ld, +1: %ld+%ld, +3: %ld+%ld\n", nA, nB, strongonly ? " strong only" : "", total, hist[0][1], hist[0][0], hist[2][1], hist[2][0], hist[4][1], hist[4][0], hist[6][1], hist[6][0]);
    return 0;
}
