/* Audit E part 2 (2026-10-15): independent re-implementation of the THM-4655 (7) certificate.
   Block parameters are computed by a DIFFERENT method from z1_cone_certificate.py (which uses subset DP):
   enumerate every permutation pi of the block D and every cut position c (0..|D|), keep (pi, c) iff pi[0..c) and pi[c..)
   are both directed paths in T[D]; then (unnormalised, i.e. times h(D)):
      h      = #{pi : pi is a path}                    e~ = #{pi path, last in I}     s~ = #{pi path, first in O}
      sigma~ = #{(pi,c)}                                (X = set of pi[0..c), any c)
      beta~  = #{(pi,c) : c < |D|, pi[c] in O}          gamma~ = #{(pi,c) : 0 < c < |D|, pi[c-1] in I, pi[c] in O}
      alpha~ = #{(pi,c) : c > 0, pi[c-1] in I}
   and the v-parts restrict to v in pi[0..c).  (O = N+(z) cap D, I = N-(z) cap D.)
   For every strong block D (singleton; gentourng -z -c classes of sizes 3..MAXS), every pattern O (2^|D|) and every v:
     (1) all four bilinear conditions (a >= 0, a+b >= 0, a+c >= 0, a+b+c+d >= 0, region P, Q >= 1) at all 16 corners
         (g, s' in {0,1}; Lp in {1,P}; Rp in {1,Q}) for 3N - H >= 0 and for 2H - 3N >= 0, exactly as in THM-4655 (7)(iv);
     (2) the closed forms L4 = 1 + sigma/2 + 2(beta(v)+gamma(v)+alpha(v)) - (bbar+gbar+abar+s) and
         U4 = 1 + sigma/2 + 2(bbar+gbar+abar) - (beta(v)+gamma(v)+alpha(v)+e); minimum slack per size (exact rational);
     (3) the general lemmas used for the other corners: beta - beta(v) <= sigma/2, alpha - alpha(v) <= sigma/2 - 1,
         alpha(v) <= sigma/2, beta(v) <= sigma/2 - 1, and the block inequalities (i): gamma + s <= beta <= sigma - 1,
         e + gamma <= alpha <= sigma - 1, and the end-block bounds (ii): O nonempty => beta >= 1, gamma + s >= 1;
         I nonempty => alpha >= 1, e + gamma >= 1.
   All arithmetic in integers (everything multiplied by h(D) > 0; sigma~ is even).
   Usage: for s in 3 4 5 6 7; do gentourng -q -z -c $s; done | ./z1_cert_audit                                        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static unsigned out[16];
static int perm[16];
/* list of valid (pi, c): store first/last of prefix and suffix and the prefix mask */
typedef struct { unsigned char pf, pl, sf, sl, c; unsigned short mask; } cut_t;
static cut_t *cuts; static long ncuts, capc;
static ll hcount; static int *pathfirst, *pathlast; static long npaths, capp;
static int is_path(int a, int b) { for (int i = a; i + 1 < b; i++) if (!(out[perm[i]] >> perm[i + 1] & 1)) return 0; return 1; }
static void gen(int k, unsigned used) {
    if (k == n) {
        if (is_path(0, n)) {
            if (npaths == capp) { capp = capp ? 2 * capp : 1024; pathfirst = realloc(pathfirst, capp * sizeof(int)); pathlast = realloc(pathlast, capp * sizeof(int)); }
            pathfirst[npaths] = perm[0]; pathlast[npaths] = perm[n - 1]; npaths++;
        }
        for (int c = 0; c <= n; c++) {
            if (!is_path(0, c) || !is_path(c, n)) continue;
            if (ncuts == capc) { capc = capc ? 2 * capc : 4096; cuts = realloc(cuts, capc * sizeof(cut_t)); }
            cut_t t; t.c = c; unsigned m = 0; for (int i = 0; i < c; i++) m |= 1u << perm[i];
            t.mask = m; t.pf = c ? perm[0] : 255; t.pl = c ? perm[c - 1] : 255; t.sf = c < n ? perm[c] : 255; t.sl = c < n ? perm[n - 1] : 255;
            cuts[ncuts++] = t;
        }
        return;
    }
    for (int v = 0; v < n; v++) if (!(used >> v & 1)) { perm[k] = v; gen(k + 1, used | 1u << v); }
}
/* bilinear nonnegativity on [1,inf)^2 with a common positive denominator */
static int bil(ll a2, ll b2, ll c2, ll d2) { /* all doubled to keep sigma/2 integral: we pass 2a,2b,... */
    return a2 >= 0 && a2 + b2 >= 0 && a2 + c2 >= 0 && a2 + b2 + c2 + d2 >= 0;
}
int main(void) {
    char line[256];
    long types = 0, fails_corner = 0, fails_lemma = 0, fails_closed = 0;
    /* min slack per size as fraction num/den */
    ll minLn[16], minLd[16], minUn[16], minUd[16]; long typesz[16];
    for (int i = 0; i < 16; i++) { minLn[i] = minUn[i] = 1LL << 60; minLd[i] = minUd[i] = 1; typesz[i] = 0; }
    int first = 1;
    for (;;) {
        if (first) { first = 0; n = 1; out[0] = 0; }          /* the singleton block first */
        else {
            if (!fgets(line, sizeof line, stdin)) break;
            if (line[0] != '&') continue;
            line[strcspn(line, "\r\n")] = 0;
            n = line[1] - 63; memset(out, 0, sizeof out);
            for (int k = 0; k < n * n; k++) { int b = line[2 + k / 6] - 63; if (b >> (5 - k % 6) & 1) out[k / n] |= 1u << (k % n); }
        }
        ncuts = 0; npaths = 0; gen(0, 0);
        ll h = npaths;
        for (unsigned O = 0; O < (1u << n); O++) {             /* pattern: O = out-neighbours of z in D */
            types++; typesz[n]++;
            ll e = 0, s = 0;
            for (long p = 0; p < npaths; p++) { if (!(O >> pathlast[p] & 1)) e++; if (O >> pathfirst[p] & 1) s++; }
            ll sig = 0, be = 0, ga = 0, al = 0;
            ll sigv[16] = {0}, bev[16] = {0}, gav[16] = {0}, alv[16] = {0};
            for (long q = 0; q < ncuts; q++) {
                cut_t *t = &cuts[q]; int c = t->c;
                int bterm = (c < n) && (O >> t->sf & 1);
                int gterm = (c > 0 && c < n) && !(O >> t->pl & 1) && (O >> t->sf & 1);
                int aterm = (c > 0) && !(O >> t->pl & 1);
                sig++; be += bterm; ga += gterm; al += aterm;
                for (int v = 0; v < n; v++) if (t->mask >> v & 1) { sigv[v]++; bev[v] += bterm; gav[v] += gterm; alv[v] += aterm; }
            }
            /* (i), (ii) */
            if (!(ga + s <= be && be <= sig - h && e + ga <= al && al <= sig - h)) fails_lemma++;
            if (O && !(be >= h && ga + s >= h)) fails_lemma++;
            if (O != (1u << n) - 1 && !(al >= h && e + ga >= h)) fails_lemma++;
            for (int v = 0; v < n; v++) {
                if (2 * sigv[v] != sig) fails_lemma++;                       /* fair coin: v-part of sigma is half */
                ll bb = be - bev[v], gb = ga - gav[v], ab = al - alv[v];
                if (!(2 * bb <= sig && 2 * ab <= sig - 2 * h && 2 * alv[v] <= sig && 2 * bev[v] <= sig - 2 * h)) fails_lemma++;
                /* (1) all corners, everything doubled (x2) and times h:  a = sigma/2 -> 2a = sig */
                for (int g = 0; g <= 1; g++) for (int sp = 0; sp <= 1; sp++) for (int LpP = 0; LpP <= 1; LpP++) for (int RpQ = 0; RpQ <= 1; RpQ++) {
                    /* 3N - H = 2Rp + 2 e s' + (3bv - be) Q + 3gv - ga - g s + P (sig Q/2 + 3 av - al) - Lp   (times h; Lp, Rp in units of h) */
                    ll A2 = sig;
                    ll B2 = 2 * (3 * alv[v] - al) - (LpP ? 2 * h : 0);
                    ll C2 = 2 * (3 * bev[v] - be) + (RpQ ? 4 * h : 0);
                    ll D2 = 2 * (2 * e * sp + 3 * gav[v] - ga - g * s) - (LpP ? 0 : 2 * h) + (RpQ ? 0 : 4 * h);
                    if (!bil(A2, B2, C2, D2)) fails_corner++;
                    /* 2H - 3N = -Rp - e s' + (2be - 3bv) Q + 2ga - 3gv + 2 g s + P (sig Q/2 + 2 al - 3 av) + 2 Lp */
                    B2 = 2 * (2 * al - 3 * alv[v]) + (LpP ? 4 * h : 0);
                    C2 = 2 * (2 * be - 3 * bev[v]) - (RpQ ? 2 * h : 0);
                    D2 = 2 * (-e * sp + 2 * ga - 3 * gav[v] + 2 * g * s) + (LpP ? 0 : 4 * h) - (RpQ ? 0 : 2 * h);
                    if (!bil(A2, B2, C2, D2)) fails_corner++;
                }
                /* (2) closed forms (times 2h) */
                ll L4 = 2 * h + sig + 4 * (bev[v] + gav[v] + alv[v]) - 2 * (bb + gb + ab + s);
                ll U4 = 2 * h + sig + 4 * (bb + gb + ab) - 2 * (bev[v] + gav[v] + alv[v] + e);
                if (L4 < 0 || U4 < 0) fails_closed++;
                ll den = 2 * h;
                if ((__int128)L4 * minLd[n] < (__int128)minLn[n] * den) { minLn[n] = L4; minLd[n] = den; }
                if ((__int128)U4 * minUd[n] < (__int128)minUn[n] * den) { minUn[n] = U4; minUd[n] = den; }
            }
        }
    }
    printf("block types (block, z-pattern) checked: %ld; corner-condition failures: %ld; lemma failures ((i),(ii), fair coin, the four comparisons): %ld; L4/U4 negative: %ld\n",
           types, fails_corner, fails_lemma, fails_closed);
    for (int sz = 1; sz < 16; sz++) if (typesz[sz]) {
        ll a = minLn[sz], b = minLd[sz], g; ll x = a < 0 ? -a : a, y = b; while (y) { ll t = x % y; x = y; y = t; } g = x ? x : 1;
        ll c = minUn[sz], d = minUd[sz]; x = c < 0 ? -c : c; y = d; while (y) { ll t = x % y; x = y; y = t; } ll g2 = x ? x : 1;
        printf("  size %d: %ld (block, pattern) types; min slack L4 = %lld/%lld, min slack U4 = %lld/%lld\n", sz, typesz[sz], a / g, b / g, c / g2, d / g2);
    }
    return 0;
}
