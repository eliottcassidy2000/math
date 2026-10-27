/* audit_lattice.c -- independent audit of THM-4820 (D): the catalytic minimum g on Z[w].  Written from scratch.
 *
 * Two independent computations:
 *
 * (1) g_W(x) = min over binary merge trees with unit leaves, root x, ALL node values of norm <= W, of sum beta+.
 *     Algorithm: Knuth's generalisation of Dijkstra's algorithm to superior grammar functions (Knuth 1977): the
 *     combination g(y) + g(z) + beta(y,z)+ is monotone and >= max(g(y), g(z)), so points can be finalised in order of
 *     cost; each finalised point is combined with every previously finalised point.  (The producer used repeated
 *     relaxation rounds instead.)  Seeds: the six units (mode 6) or only 1, w, w^2 (mode 3; tests (A4)).
 *
 * (2) The window-free tight closure.  Phi = floor(Q/3); a merge is tight iff beta+ + Phi(y) + Phi(z) - Phi(y+z) = 0.
 *     G = least set containing 0 and the units and closed under 'x = y + z tight, y, z != 0, y, z in G'.
 *     Audit lemma: g(x) = Phi(x) iff x in G, and every all-tight tree with root x has all node norms
 *     <= Q(x) + [3 | Q(x)].  So membership of x in G is decided by points of norm <= Q(x) + 1: window-free.
 *     Computed on all points of norm <= N1 (exact for every x with Q(x) + [3|Q(x)] <= N1).
 *
 * Usage: ./audit_lattice W N1 mode(6|3) [dump] [tree a b]      (N1 = 0 skips the closure; tree writes audit_value_tree.txt)
 * Output lines:  P a b Q gW inG     (inG = -1 if the closure was not computed / not exact for that point)
 *                N n gW(n)           for rational n with n^2 <= W
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int Qf(int a, int b) { return a * a - a * b + b * b; }
static int R;          /* coordinate radius: all points of norm <= max(W,N1) have |a|,|b| <= R */
static int SIDE;
static inline int ix(int a, int b) { return (a + R) * SIDE + (b + R); }

typedef struct { int c, p; } HE;
static HE *heap; static long hn = 0, hcap = 0;
static void hpush(int c, int p) { if (hn == hcap) { hcap = hcap ? 2 * hcap : 1 << 16; heap = realloc(heap, sizeof(HE) * hcap); }
    long i = hn++; heap[i].c = c; heap[i].p = p;
    while (i > 0) { long q = (i - 1) / 2; if (heap[q].c <= heap[i].c) break; HE t = heap[q]; heap[q] = heap[i]; heap[i] = t; i = q; } }
static HE hpop(void) { HE top = heap[0]; heap[0] = heap[--hn]; long i = 0;
    for (;;) { long l = 2 * i + 1, r = l + 1, m = i; if (l < hn && heap[l].c < heap[m].c) m = l; if (r < hn && heap[r].c < heap[m].c) m = r;
        if (m == i) break; HE t = heap[m]; heap[m] = heap[i]; heap[i] = t; i = m; } return top; }

int main(int argc, char **argv) {
    int W = atoi(argv[1]), N1 = atoi(argv[2]), mode = atoi(argv[3]), dump = 0, want_tree = 0, tra = 0, trb = 0;
    for (int i = 4; i < argc; i++) { if (!strcmp(argv[i], "dump")) dump = 1;
        if (!strcmp(argv[i], "tree") && i + 2 < argc) { want_tree = 1; tra = atoi(argv[i + 1]); trb = atoi(argv[i + 2]); i += 2; } }
    int M = W > N1 ? W : N1;
    R = 1; while (3 * R * R <= 4 * M + 8) R++; R += 2;   /* Q(a,b) >= 3a^2/4 and >= 3b^2/4 */
    SIDE = 2 * R + 1;
    int S = SIDE * SIDE;
    int *g = malloc(sizeof(int) * S); char *fin = calloc(S, 1);
    int *dya = malloc(sizeof(int) * S), *dyb = malloc(sizeof(int) * S);   /* derivation: g[t] achieved by y + (t - y) */
    const int INF = 1 << 29;
    for (int i = 0; i < S; i++) g[i] = INF;
    int units6[6][2] = {{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,-1}};   /* 1, -1, w, -w, -w^2, w^2 */
    int units3[3][2] = {{1,0},{0,1},{-1,-1}};                       /* 1, w, w^2 */
    if (mode == 6) for (int i = 0; i < 6; i++) { g[ix(units6[i][0], units6[i][1])] = 0; hpush(0, ix(units6[i][0], units6[i][1])); }
    else for (int i = 0; i < 3; i++) { g[ix(units3[i][0], units3[i][1])] = 0; hpush(0, ix(units3[i][0], units3[i][1])); }
    int *fa = malloc(sizeof(int) * S), *fb = malloc(sizeof(int) * S), *fg = malloc(sizeof(int) * S), *fq = malloc(sizeof(int) * S); int nf = 0;
    long combos = 0;
    while (hn) {
        HE e = hpop(); int p = e.p; if (fin[p] || e.c != g[p]) continue;
        fin[p] = 1; int xa = p / SIDE - R, xb = p % SIDE - R, qx = Qf(xa, xb), gx = g[p];
        fa[nf] = xa; fb[nf] = xb; fg[nf] = gx; fq[nf] = qx; nf++;
        for (int j = 0; j < nf; j++) {   /* combine with every finalised point, including itself */
            int sa = xa + fa[j], sb = xb + fb[j]; int qs = Qf(sa, sb); if (qs > W) continue;
            combos++;
            int t = ix(sa, sb); if (fin[t]) continue;
            int beta = qs - qx - fq[j]; int c = gx + fg[j] + (beta > 0 ? beta : 0);
            if (c < g[t]) { g[t] = c; dya[t] = xa; dyb[t] = xb; hpush(c, t); }
        }
    }
    fprintf(stderr, "Knuth-Dijkstra: W=%d mode=%d finalised %d points, %ld combinations\n", W, mode, nf, combos);

    /* tight closure on norms <= N1 */
    signed char *inG = malloc(S); memset(inG, -1, S);
    if (N1 > 0) {
        int *pa = malloc(sizeof(int) * S), *pb = malloc(sizeof(int) * S), np = 0;
        for (int a = -R; a <= R; a++) for (int b = -R; b <= R; b++) if (Qf(a, b) <= N1) { pa[np] = a; pb[np] = b; np++; }
        /* tight splits: for each x, all y != 0, z = x - y != 0, both of norm <= N1, with zero slack */
        long cap = 1 << 20, nt = 0; int *tx = malloc(sizeof(int) * cap), *ty = malloc(sizeof(int) * cap), *tz = malloc(sizeof(int) * cap);
        for (int i = 0; i < np; i++) {
            int xa = pa[i], xb = pb[i], qx = Qf(xa, xb);
            for (int j = 0; j < np; j++) {
                int ya = pa[j], yb = pb[j]; if (ya == 0 && yb == 0) continue;
                int za = xa - ya, zb = xb - yb; if (za == 0 && zb == 0) continue;
                int qy = Qf(ya, yb), qz = Qf(za, zb); if (qz > N1) continue;
                int beta = qx - qy - qz; int slack = (beta > 0 ? beta : 0) + qy / 3 + qz / 3 - qx / 3;
                if (slack < 0) { fprintf(stderr, "NEGATIVE SLACK?!\n"); return 2; }
                if (slack == 0) { if (nt == cap) { cap *= 2; tx = realloc(tx, sizeof(int) * cap); ty = realloc(ty, sizeof(int) * cap); tz = realloc(tz, sizeof(int) * cap); }
                    tx[nt] = ix(xa, xb); ty[nt] = ix(ya, yb); tz[nt] = ix(za, zb); nt++; }
            }
        }
        char *G = calloc(S, 1);
        G[ix(0, 0)] = 1; for (int i = 0; i < 6; i++) G[ix(units6[i][0], units6[i][1])] = 1;
        int rounds = 0, changed = 1;
        while (changed) { changed = 0; rounds++;
            for (long t = 0; t < nt; t++) if (!G[tx[t]] && G[ty[t]] && G[tz[t]]) { G[tx[t]] = 1; changed = 1; } }
        fprintf(stderr, "tight closure: %d points of norm <= %d, %ld tight splits, %d rounds\n", np, N1, nt, rounds);
        for (int i = 0; i < np; i++) { int q = Qf(pa[i], pb[i]); int B = q + (q % 3 == 0 ? 1 : 0);
            if (B <= N1) inG[ix(pa[i], pb[i])] = G[ix(pa[i], pb[i])]; }
    }
    if (want_tree) {   /* derivation DAG of (tra, trb): lines 'D xa xb = ya yb + za zb' (each value once), units are leaves */
        FILE *fo = fopen("audit_value_tree.txt", "w"); char *seen = calloc(S, 1); int *st = malloc(sizeof(int) * S); int sp = 0;
        st[sp++] = ix(tra, trb); seen[ix(tra, trb)] = 1;
        while (sp) { int t = st[--sp]; int ta = t / SIDE - R, tb = t % SIDE - R;
            int isunit = 0; for (int i = 0; i < 6; i++) if (ta == units6[i][0] && tb == units6[i][1]) isunit = 1;
            if (isunit && (mode == 6 || (ta == 1 && tb == 0) || (ta == 0 && tb == 1) || (ta == -1 && tb == -1))) continue;
            int ya = dya[t], yb = dyb[t], za = ta - ya, zb = tb - yb;
            fprintf(fo, "D %d %d = %d %d + %d %d\n", ta, tb, ya, yb, za, zb);
            int c1 = ix(ya, yb), c2 = ix(za, zb); if (!seen[c1]) { seen[c1] = 1; st[sp++] = c1; } if (!seen[c2]) { seen[c2] = 1; st[sp++] = c2; } }
        fclose(fo); fprintf(stderr, "derivation of (%d,%d) written, g = %d\n", tra, trb, g[ix(tra, trb)]); }
    for (int n = 1; n * n <= W; n++) printf("N %d %d\n", n, g[ix(n, 0)]);
    if (dump) for (int a = -R; a <= R; a++) for (int b = -R; b <= R; b++) { int q = Qf(a, b); if (q <= M)
        printf("P %d %d %d %d %d\n", a, b, q, q <= W ? g[ix(a, b)] : -1, inG[ix(a, b)]); }
    return 0;
}
