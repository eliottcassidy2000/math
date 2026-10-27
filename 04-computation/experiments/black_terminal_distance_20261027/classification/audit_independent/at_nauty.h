/* at_nauty.h -- nauty helpers for the auditor (digraph canonical forms, colored automorphism group sizes).
 * Compile with -DMAXN=WORDSIZE and link /opt/homebrew/opt/nauty/lib/libnauty1.a */
#ifndef AT_NAUTY_H
#define AT_NAUTY_H
#include <math.h>
#include "nauty.h"
#include "nautinv.h"
#include "at.h"

/* canonical form of T (as a digraph) into cg[0..n-1] (m = 1 setword per row) */
static void nt_canon(const tour *T, graph *cg) {
    int n = T->n, m = 1;
    graph g[MAXN];
    int lab[MAXN], ptn[MAXN], orb[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(opt);
    statsblk st;
    opt.getcanon = TRUE; opt.defaultptn = TRUE;
    EMPTYGRAPH(g, m, n);
    for (int u = 0; u < n; u++) for (uint32_t c = T->o[u]; c; c &= c - 1) ADDONEARC(g, u, lowbit(c), m);
    densenauty(g, lab, ptn, orb, &opt, &st, m, n, cg);
}
static int nt_same(const graph *a, const graph *b, int n) { for (int i = 0; i < n; i++) if (a[i] != b[i]) return 0; return 1; }

/* canonical code string of T (hex of canonical rows) for hashing/comparison */
static void nt_canon_str(const tour *T, char *buf) {
    graph cg[MAXN]; nt_canon(T, cg);
    int k = 0;
    for (int i = 0; i < T->n; i++) k += sprintf(buf + k, "%llx.", (unsigned long long)cg[i]);
    buf[k] = 0;
}

static int nt_selfconverse(const tour *T) {
    tour C; t_converse(T, &C);
    graph a[MAXN], b[MAXN];
    nt_canon(T, a); nt_canon(&C, b);
    return nt_same(a, b, T->n);
}

/* |Aut(T; colouring)| where col[v] is a small colour id (vertices of equal colour may be exchanged) */
static double nt_autsize_col(const tour *T, const int *col) {
    int n = T->n, m = 1;
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orb[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(opt);
    statsblk st;
    opt.getcanon = FALSE; opt.defaultptn = FALSE;
    EMPTYGRAPH(g, m, n);
    for (int u = 0; u < n; u++) for (uint32_t c = T->o[u]; c; c &= c - 1) ADDONEARC(g, u, lowbit(c), m);
    int k = 0, maxc = 0;
    for (int v = 0; v < n; v++) if (col[v] > maxc) maxc = col[v];
    for (int c = 0; c <= maxc; c++) {
        int start = k;
        for (int v = 0; v < n; v++) if (col[v] == c) lab[k++] = v;
        if (k > start) { for (int i = start; i < k; i++) ptn[i] = 1; ptn[k - 1] = 0; }
    }
    densenauty(g, lab, ptn, orb, &opt, &st, m, n, cg);
    return st.grpsize1 * pow(10.0, st.grpsize2);
}
static double nt_autsize(const tour *T) { int col[MAXV] = {0}; return nt_autsize_col(T, col); }
#endif
