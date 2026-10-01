// opus S15 (2026-10-01), helper for shaved_tournaments_20261001.py (check T3, n = 9).
// Usage: closing m classfile  -- classfile holds one canonical code per line for the m-vertex classes
// (bit b over pairs i<j in lexicographic order: 1 means i->j).
// For every base class on m vertices (codes in file) and every extension by one vertex,
// test whether the (m+1)-tournament has a non-closing Hamiltonian path (start -> end).
// Prints every all-closing tournament found (as out-neighbour masks).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int has_nonclosing(int n, const uint32_t *out) {
    static __thread uint16_t st[1<<12][12];
    int full = (1<<n)-1;
    memset(st, 0, sizeof(uint16_t)*12*(1<<n));
    uint32_t in[12]={0};
    for (int v=0; v<n; v++) for (int w=0; w<n; w++) if (out[v]>>w&1) in[w]|=1u<<v;
    for (int v=0; v<n; v++) st[1<<v][v] = 1<<v;
    for (int m=1; m<=full; m++) for (int e=0; e<n; e++) {
        uint16_t s = st[m][e]; if (!s) continue;
        uint32_t o = out[e] & ~m;
        while (o) { int w = __builtin_ctz(o); o &= o-1; st[m|1<<w][w] |= s; }
    }
    for (int e=0; e<n; e++) if (st[full][e] & in[e]) return 1;
    return 0;
}
int main(int argc, char **argv) {
    int m = atoi(argv[1]); FILE *f = fopen(argv[2], "r");
    unsigned long long code; int n = m+1; long long tested=0, bad=0;
    unsigned long long *codes = malloc(sizeof(unsigned long long)*300000); int nc=0;
    while (fscanf(f, "%llu", &code)==1) codes[nc++]=code;
    #pragma omp parallel for schedule(dynamic) reduction(+:tested,bad)
    for (int c=0; c<nc; c++) {
        uint32_t base[12]={0}; int b=0;
        for (int i=0;i<m;i++) for (int j=i+1;j<m;j++) { if (codes[c]>>b&1) base[i]|=1u<<j; else base[j]|=1u<<i; b++; }
        for (int S=0; S<(1<<m); S++) {
            uint32_t out[12]; for (int v=0; v<m; v++) out[v]=base[v]; out[m]=0;
            for (int v=0; v<m; v++) { if (S>>v&1) out[m]|=1u<<v; else out[v]|=1u<<m; }
            tested++;
            if (!has_nonclosing(n, out)) {
                bad++;
                #pragma omp critical
                { printf("ALLCLOSING n=%d:", n); for (int v=0; v<n; v++) printf(" %u", out[v]); printf("\n"); }
            }
        }
    }
    printf("n=%d tested=%lld allclosing=%lld\n", n, tested, bad);
    return 0;
}
