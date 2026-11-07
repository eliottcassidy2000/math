/* Exhaustive enumerator for aud_perm.py part B (auditor's own code).
   Input file: "nv ne nc" then ne lines "u v comp parity".  For each assignment of the nc class bits,
   o_e = bit[comp] ^ parity, o_e = 1 means u -> v (points into v).  Counts sinks (out-degree 0) and tests
   acyclicity by Kahn's algorithm.  Prints the sink spectrum {count: #assignments} for all and for acyclic ones. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r");
    int nv, ne, nc;
    if (fscanf(f, "%d %d %d", &nv, &ne, &nc) != 3) return 1;
    int *U = malloc(sizeof(int) * ne), *V = malloc(sizeof(int) * ne), *Cc = malloc(sizeof(int) * ne), *Pp = malloc(sizeof(int) * ne);
    for (int e = 0; e < ne; e++) if (fscanf(f, "%d %d %d %d", &U[e], &V[e], &Cc[e], &Pp[e]) != 4) return 1;
    fclose(f);
    long long all[1024], acy[1024];
    memset(all, 0, sizeof(all)); memset(acy, 0, sizeof(acy));
    int *outdeg = malloc(sizeof(int) * nv), *indeg = malloc(sizeof(int) * nv), *stack = malloc(sizeof(int) * nv);
    int *tail = malloc(sizeof(int) * ne), *head = malloc(sizeof(int) * ne);
    /* adjacency by edge list per vertex */
    int *deg = calloc(nv, sizeof(int));
    for (int e = 0; e < ne; e++) { deg[U[e]]++; deg[V[e]]++; }
    int **inc = malloc(sizeof(int *) * nv);
    for (int v = 0; v < nv; v++) { inc[v] = malloc(sizeof(int) * deg[v]); deg[v] = 0; }
    for (int e = 0; e < ne; e++) { inc[U[e]][deg[U[e]]++] = e; inc[V[e]][deg[V[e]]++] = e; }
    for (long long M = 0; M < (1ll << nc); M++) {
        for (int v = 0; v < nv; v++) { outdeg[v] = 0; indeg[v] = 0; }
        for (int e = 0; e < ne; e++) {
            int o = (int)((M >> Cc[e]) & 1) ^ Pp[e];
            if (o) { tail[e] = U[e]; head[e] = V[e]; } else { tail[e] = V[e]; head[e] = U[e]; }
            outdeg[tail[e]]++; indeg[head[e]]++;
        }
        int s = 0;
        for (int v = 0; v < nv; v++) if (outdeg[v] == 0) s++;
        all[s]++;
        int top = 0, seen = 0;
        for (int v = 0; v < nv; v++) if (indeg[v] == 0) stack[top++] = v;
        while (top > 0) {
            int x = stack[--top]; seen++;
            for (int j = 0; j < deg[x]; j++) { int e = inc[x][j]; if (tail[e] == x) { int h = head[e]; if (--indeg[h] == 0) stack[top++] = h; } }
        }
        if (seen == nv) acy[s]++;
    }
    long long tot_ac = 0; for (int s = 0; s < 1024; s++) tot_ac += acy[s];
    printf("all %lld:", 1ll << nc);
    for (int s = 0; s < 1024; s++) if (all[s]) printf(" %d:%lld", s, all[s]);
    printf("\nacyclic %lld:", tot_ac);
    for (int s = 0; s < 1024; s++) if (acy[s]) printf(" %d:%lld", s, acy[s]);
    printf("\n");
    return 0;
}
