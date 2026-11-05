/* ribbon_levels.c -- breadth-first enumeration of span-w ribbon words by length, deduplicated on the
   observable integer state (export of intrep.py) and the condensation state.  Two prefixes with the
   same reduced state have the same H for every continuation, so merging them loses no value.
   Per length n it records the set of H values of ALL words (for the tiling law) and of the STRONG words
   (accepting condensation states), for H <= N; prefixes with H > N are dropped (insertion lemma).
   usage: ribbon_levels <intrep file> <automaton file> <w> <N> <maxlen> <out prefix> [maxstates]        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#define MAXR 64
#define MAXL 8
#define MAXC 256
static int r, nl, w, nc, na, maxlen;
static int lmask[MAXL];
static int64_t A[MAXL][MAXR][MAXR];
/* state coordinates are int32 (checked for overflow); hash entries are compact */
static int64_t y0v[MAXR], hcoef[MAXR];
static int trans[MAXC][MAXL];
static uint8_t accept[MAXC];
static uint64_t N;
typedef struct { int32_t y[MAXR]; int8_t c; } state_t;
/* open-addressing hash set of states */
typedef struct { state_t *tab; uint8_t *used; size_t cap, n; } hset;
static uint64_t hash_state(const state_t *s) {
    uint64_t h = 1469598103934665603ull;
    for (int k = 0; k < r; k++) { h ^= (uint64_t)s->y[k]; h *= 1099511628211ull; h ^= h >> 29; }
    h ^= (uint64_t)s->c * 0x9E3779B97F4A7C15ull; h ^= h >> 31;
    return h;
}
static void hs_init(hset *H, size_t cap) { H->cap = cap; H->n = 0; H->tab = malloc(cap * sizeof(state_t)); H->used = calloc(cap, 1); if (!H->tab || !H->used) { fprintf(stderr, "alloc fail\n"); exit(3); } }
static void hs_free(hset *H) { free(H->tab); free(H->used); }
static int hs_insert_raw(hset *H, const state_t *s) {
    uint64_t h = hash_state(s) % H->cap;
    while (H->used[h]) {
        if (H->tab[h].c == s->c && memcmp(H->tab[h].y, s->y, sizeof(int32_t) * r) == 0) return 0;
        h = (h + 1) % H->cap;
    }
    H->used[h] = 1; H->tab[h] = *s; H->n++;
    return 1;
}
static void hs_grow(hset *H) {
    hset G; hs_init(&G, H->cap * 2);
    for (size_t i = 0; i < H->cap; i++) if (H->used[i]) hs_insert_raw(&G, &H->tab[i]);
    hs_free(H); *H = G;
}
static int hs_insert(hset *H, const state_t *s) {
    if (H->n * 10 >= H->cap * 7) hs_grow(H);
    return hs_insert_raw(H, s);
}
int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s intrepfile automatonfile w N maxlen outprefix [maxstates]\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("intrep"); return 1; }
    w = atoi(argv[3]); N = strtoull(argv[4], 0, 10); maxlen = atoi(argv[5]);
    size_t maxstates = argc > 7 ? strtoull(argv[7], 0, 10) : 20000000;
    int skip_empty = argc > 8 ? atoi(argv[8]) : 0;   /* 1 = restrict to nonempty letters (the 3-letter alphabet for span 2) */
    if (fscanf(f, "%d %d", &r, &nl) != 2) return 1;
    for (int k = 0; k < r; k++) if (fscanf(f, "%lld", (long long *)&y0v[k]) != 1) return 1;
    for (int k = 0; k < r; k++) if (fscanf(f, "%lld", (long long *)&hcoef[k]) != 1) return 1;
    for (int li = 0; li < nl; li++) {
        if (fscanf(f, "%d", &lmask[li]) != 1) return 1;
        for (int k = 0; k < r; k++) for (int l = 0; l < r; l++) if (fscanf(f, "%lld", (long long *)&A[li][k][l]) != 1) return 1;
    }
    fclose(f);
    f = fopen(argv[2], "r"); if (!f) { perror("automaton"); return 1; }
    if (fscanf(f, "%d %d", &nc, &na) != 2) return 1;
    for (int k = 0; k < na; k++) { int s; if (fscanf(f, "%d", &s) != 1) return 1; accept[s] = 1; }
    for (int s = 0; s < nc; s++) for (int li = 0; li < nl; li++) if (fscanf(f, "%d", &trans[s][li]) != 1) return 1;
    fclose(f);
    uint64_t nodd = N / 2 + 1;
    uint8_t *att_all = calloc(nodd / 8 + 1, 1), *att_strong = calloc(nodd / 8 + 1, 1), *level = calloc(nodd / 8 + 1, 1);
    att_strong[0] |= 1 | 2 | 4;  /* 1, 3, 5: K1, C3, the strong 4-tournament */
    char name[512]; snprintf(name, sizeof name, "%s_levels.txt", argv[6]);
    FILE *o = fopen(name, "w");
    fprintf(o, "# span-%d ribbon words by length, deduplicated on (observable integer state of rank %d, condensation state); H <= %llu\n", w, r, (unsigned long long)N);
    fprintf(o, "# n  distinct_states  values_n  min_n  max_n  solid_lo  solid_hi  misses_inside_[min,max]  strong_values_n  cumulative_odd_misses_all  cumulative_odd_misses_strong\n");
    hset cur, nxt; hs_init(&cur, 1024);
    state_t s0; memset(&s0, 0, sizeof s0); for (int k = 0; k < r; k++) s0.y[k] = (int32_t)y0v[k]; s0.c = 0;
    hs_insert(&cur, &s0);
    clock_t t0 = clock();
    for (int n = 1; n <= maxlen; n++) {
        size_t cap = cur.n * 2 + 4096;
        hs_init(&nxt, cap);
        memset(level, 0, nodd / 8 + 1);
        uint64_t vmin = ~0ull, vmax = 0, nstrong_vals = 0;
        for (size_t i = 0; i < cur.cap; i++) if (cur.used[i]) {
            const state_t *s = &cur.tab[i];
            for (int li = 0; li < nl; li++) {
                int m = lmask[li];
                if (skip_empty && n >= 2 && m == 0) continue;
                if (n >= 2 && (m >> (n - 1))) continue;
                if (n == 1 && m) continue;
                int c2 = (s->c < 0) ? -1 : trans[s->c][li];
                state_t t; t.c = c2;
                int64_t acc64[MAXR]; int ovf = 0;
                for (int k = 0; k < r; k++) { int64_t acc = 0; for (int l = 0; l < r; l++) acc += A[li][k][l] * (int64_t)s->y[l]; acc64[k] = acc; if (acc > 2000000000LL || acc < -2000000000LL) ovf = 1; }
                int64_t H = 0; for (int k = 0; k < r; k++) H += hcoef[k] * acc64[k];
                if (H <= 0 || (uint64_t)H > N) continue;
                if (ovf) { fprintf(stderr, "int32 overflow at n=%d\n", n); return 5; }
                for (int k = 0; k < r; k++) t.y[k] = (int32_t)acc64[k];
                if (!(H & 1)) { fprintf(stderr, "even H\n"); return 2; }
                uint64_t idx = (uint64_t)H >> 1;
                level[idx >> 3] |= 1u << (idx & 7); att_all[idx >> 3] |= 1u << (idx & 7);
                if ((uint64_t)H < vmin) vmin = H; if ((uint64_t)H > vmax) vmax = H;
                if (c2 >= 0 && accept[c2]) att_strong[idx >> 3] |= 1u << (idx & 7);
                /* dead condensation states are kept (c = -1) for the all-words value set */
                if (nxt.n < maxstates) hs_insert(&nxt, &t);
                else { static int warned = 0; if (!warned) { fprintf(stderr, "state cap reached at n=%d\n", n); warned = 1; } }
            }
        }
        /* statistics for this level */
        uint64_t nvals = 0, solid_lo = 0, solid_hi = 0, best = 0, start = 0, inrun = 0, misses_inside = 0;
        for (uint64_t i = 0; i < nodd; i++) {
            int on = (level[i >> 3] >> (i & 7)) & 1;
            if (on) nvals++;
            if (2 * i + 1 >= vmin && 2 * i + 1 <= vmax && !on) misses_inside++;
            if (on && !inrun) { inrun = 1; start = i; }
            if (!on && inrun) { inrun = 0; if (i - start > best) { best = i - start; solid_lo = 2 * start + 1; solid_hi = 2 * (i - 1) + 1; } }
        }
        if (inrun && nodd - start > best) { best = nodd - start; solid_lo = 2 * start + 1; solid_hi = 2 * (nodd - 1) + 1; }
        uint64_t cum_all = 0, cum_strong = 0;
        for (uint64_t i = 0; i < nodd; i++) { if (!((att_all[i >> 3] >> (i & 7)) & 1)) cum_all++; if (!((att_strong[i >> 3] >> (i & 7)) & 1)) cum_strong++; }
        fprintf(o, "%d %zu %llu %llu %llu %llu %llu %llu %llu %llu\n", n, nxt.n, (unsigned long long)nvals, (unsigned long long)vmin, (unsigned long long)vmax,
                (unsigned long long)solid_lo, (unsigned long long)solid_hi, (unsigned long long)misses_inside, (unsigned long long)cum_all, (unsigned long long)cum_strong);
        fflush(o);
        fprintf(stderr, "n=%d states %zu values %llu [%llu,%llu] solid [%llu,%llu] cum misses all %llu strong %llu  (%.0f s)\n", n, nxt.n, (unsigned long long)nvals,
                (unsigned long long)vmin, (unsigned long long)vmax, (unsigned long long)solid_lo, (unsigned long long)solid_hi, (unsigned long long)cum_all, (unsigned long long)cum_strong, (double)(clock() - t0) / CLOCKS_PER_SEC);
        hs_free(&cur); cur = nxt;
        if (cur.n == 0) break;
    }
    fprintf(o, "# odd misses <= N (all words):");
    for (uint64_t i = 0; i < nodd; i++) if (!((att_all[i >> 3] >> (i & 7)) & 1)) fprintf(o, " %llu", (unsigned long long)(2 * i + 1));
    fprintf(o, "\n# odd misses <= N (strong words):");
    for (uint64_t i = 0; i < nodd; i++) if (!((att_strong[i >> 3] >> (i & 7)) & 1)) fprintf(o, " %llu", (unsigned long long)(2 * i + 1));
    fprintf(o, "\n");
    fclose(o);
    return 0;
}
