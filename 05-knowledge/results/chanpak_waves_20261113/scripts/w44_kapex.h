/* Wave 44: exact Hamiltonian-path counts for k-apex tapes (k <= 4).
   Tape: spine vertices 0..m-1 with i -> j for i < j (transitive); apexes 0..k-1.
   beat[i] bit a = 1 iff spine vertex i beats apex a;  arc[a][b] = 1 iff apex a -> apex b (a != b).
   A Hamiltonian path visits the apexes in some order pi(1..k) and cuts the spine into k+1 increasing runs
   R0 pi1 R1 pi2 ... pik Rk.  Conditions: R0's last vertex beats pi1; R_j (1<=j<k) starts at a vertex beaten by pi_j
   and ends at a vertex beating pi_(j+1), and if empty needs pi_j -> pi_(j+1); Rk starts at a vertex beaten by pi_k.
   Scan the spine; each vertex joins one run.  Per-run states: R0 {legal, illegal}; middle {empty, legal, illegal};
   Rk {empty, nonempty}.  2 * 3^(k-1) * 2 states per order, k! orders. */
#include <string.h>
typedef unsigned long long ull;
static int K_, M_;
static unsigned beat_[64];
static int arc_[4][4];
static int perm_[24][4], nperm_;
static void kapex_init_perms(int k) {
    int p[4] = {0, 1, 2, 3}; nperm_ = 0;
    /* generate all permutations of 0..k-1 (Heap-free simple recursion via nested loops for k <= 4) */
    int a, b, c, d;
    for (a = 0; a < k; a++) for (b = 0; b < k; b++) for (c = 0; c < k; c++) for (d = 0; d < k; d++) {
        int q[4] = {a, b, c, d}; int ok = 1, seen = 0;
        for (int t = 0; t < k; t++) { if (seen >> q[t] & 1) { ok = 0; break; } seen |= 1 << q[t]; }
        if (!ok) continue;
        int dup = 0; for (int s = 0; s < nperm_; s++) { int same = 1; for (int t = 0; t < k; t++) if (perm_[s][t] != q[t]) same = 0; if (same) dup = 1; }
        if (dup) continue;
        for (int t = 0; t < k; t++) perm_[nperm_][t] = q[t];
        nperm_++;
        if (k < 4) { /* the unused nested loops over-generate; dedupe handles it */ }
    }
    (void)p;
}
/* state encoding: run 0 in {0,1} (0 legal), runs 1..k-1 in {0 empty, 1 legal, 2 illegal}, run k in {0,1} */
static ull kapex_count(void) {
    int k = K_, m = M_;
    int nst = 2 * 2; for (int j = 1; j < k; j++) nst *= 3;
    static ull dp[2][216];
    ull total = 0;
    for (int s = 0; s < nperm_; s++) {
        const int *pi = perm_[s];
        memset(dp, 0, sizeof dp);
        dp[0][0] = 1;      /* all runs empty / R0 legal */
        int cur = 0;
        for (int i = 0; i < m; i++) {
            int nx = cur ^ 1; memset(dp[nx], 0, sizeof(ull) * nst);
            unsigned bi = beat_[i];
            for (int st = 0; st < nst; st++) {
                ull x = dp[cur][st]; if (!x) continue;
                /* decode */
                int r[5]; int t = st; r[0] = t % 2; t /= 2;
                for (int j = 1; j < k; j++) { r[j] = t % 3; t /= 3; }
                r[k] = t % 2;
                for (int j = 0; j <= k; j++) {
                    int q[5]; memcpy(q, r, sizeof q);
                    if (j == 0) {
                        q[0] = (bi >> pi[0] & 1) ? 0 : 1;
                    } else if (j < k) {
                        if (r[j] == 0 && (bi >> pi[j - 1] & 1)) continue;          /* entry: pi_j must beat the vertex */
                        q[j] = (bi >> pi[j] & 1) ? 1 : 2;                         /* last vertex beats pi_(j+1)? */
                    } else {
                        if (r[k] == 0 && (bi >> pi[k - 1] & 1)) continue;          /* entry: pi_k must beat the vertex */
                        q[k] = 1;
                    }
                    int ns = q[k]; for (int jj = k - 1; jj >= 1; jj--) ns = ns * 3 + q[jj]; ns = ns * 2 + q[0];
                    dp[nx][ns] += x;
                }
            }
            cur = nx;
        }
        for (int st = 0; st < nst; st++) {
            ull x = dp[cur][st]; if (!x) continue;
            int r[5]; int t = st; r[0] = t % 2; t /= 2;
            for (int j = 1; j < k; j++) { r[j] = t % 3; t /= 3; }
            r[k] = t % 2;
            if (r[0] != 0) continue;
            int ok = 1;
            for (int j = 1; j < k; j++) {
                if (r[j] == 0) { if (!arc_[pi[j - 1]][pi[j]]) { ok = 0; break; } }
                else if (r[j] != 1) { ok = 0; break; }
            }
            if (ok) total += x;
        }
    }
    return total;
}
