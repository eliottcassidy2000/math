/* extgraph10.c: the one-vertex extension graph between the 191,536 classes of order 9 and the 9,733,056 classes of order 10.
   Canonical forms by nauty's labelg (digraph6).  Streams gentourng -q -z 10, writes the ten deletions of each class as
   digraph6 to a chunk file, canonises the chunk with labelg, and maps each canonical string to its order-9 class by a
   hash table.  Reports contacts (distinct parents per class), children per order-9 class, the transitive tournament's
   children, and the canonical-parent tree (least canonical deletion): max branching, childless classes, mean. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define N9 191536
#define HS (1u << 19)
static char can9[N9][20]; static int htab[HS];
static unsigned hashs(const char *s) { unsigned h = 2166136261u; while (*s) { h ^= (unsigned char)*s++; h *= 16777619u; } return h; }
static int lookup(const char *s) { unsigned i = hashs(s) & (HS - 1); while (htab[i] >= 0) { if (!strcmp(can9[htab[i]], s)) return htab[i]; i = (i + 1) & (HS - 1); } return -1; }
static void encode9(u64 *adj, char *out) {   /* 9 vertices: 81 bits row-major, 14 chars */
    int bit = 0; out[0] = '&'; out[1] = 9 + 63; int p = 2; int acc = 0, nb = 0;
    for (int i = 0; i < 9; i++) for (int j = 0; j < 9; j++) { acc = (acc << 1) | (int)((adj[i] >> j) & 1); nb++; if (nb == 6) { out[p++] = acc + 63; acc = 0; nb = 0; } bit++; }
    if (nb) { acc <<= (6 - nb); out[p++] = acc + 63; }
    out[p] = 0;
}
int main(void) {
    char line[256];
    FILE *g = popen("gentourng -q -z 9 | labelg -q", "r"); int k = 0;
    while (fgets(line, sizeof line, g)) { int len = strlen(line); while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0; if (!len) continue; strcpy(can9[k++], line); }
    pclose(g); if (k != N9) { fprintf(stderr, "order 9: %d classes\n", k); return 1; }
    for (unsigned i = 0; i < HS; i++) htab[i] = -1;
    for (int i = 0; i < N9; i++) { unsigned j = hashs(can9[i]) & (HS - 1); while (htab[j] >= 0) j = (j + 1) & (HS - 1); htab[j] = i; }
    /* transitive M_9 */
    u64 tr[9]; for (int i = 0; i < 9; i++) { tr[i] = 0; for (int j = i + 1; j < 9; j++) tr[i] |= 1ULL << j; }
    char t9[20]; encode9(tr, t9); FILE *tf = fopen("tmp_tr.d6", "w"); fprintf(tf, "%s\n", t9); fclose(tf);
    system("labelg -q tmp_tr.d6 > tmp_tr.can"); tf = fopen("tmp_tr.can", "r"); fgets(line, sizeof line, tf); fclose(tf);
    { int len = strlen(line); while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0; }
    int trans = lookup(line); if (trans < 0) { fprintf(stderr, "transitive not found\n"); return 1; }
    static long children[N9], branch[N9]; long hist[12] = {0}; long classes = 0, contacts = 0;
    const int CH = 250000; static char (*buf)[20]; buf = malloc((size_t)CH * 10 * 20); static int cls[10];
    FILE *s = popen("gentourng -q -z 10", "r"); int done = 0;
    while (!done) {
        FILE *cf = fopen("tmp_chunk.d6", "w"); int nc = 0;
        while (nc < CH) {
            if (!fgets(line, sizeof line, s)) { done = 1; break; }
            int len = strlen(line); while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0; if (len < 2 || line[0] != '&') continue;
            u64 adj[10]; for (int i = 0; i < 10; i++) adj[i] = 0; int bit = 0;
            for (int p = 2; p < len; p++) { int c = line[p] - 63; for (int b = 5; b >= 0; b--, bit++) { if (bit >= 100) break; if ((c >> b) & 1) adj[bit / 10] |= 1ULL << (bit % 10); } }
            for (int v = 0; v < 10; v++) {
                u64 d[9]; int r = 0;
                for (int i = 0; i < 10; i++) { if (i == v) continue; u64 row = 0; int c2 = 0; for (int j = 0; j < 10; j++) { if (j == v) continue; if ((adj[i] >> j) & 1) row |= 1ULL << c2; c2++; } d[r++] = row; }
                char e9[20]; encode9(d, e9); fputs(e9, cf); fputc('\n', cf);
            }
            nc++;
        }
        fclose(cf); if (!nc) break;
        system("labelg -q tmp_chunk.d6 > tmp_chunk.can");
        FILE *rf = fopen("tmp_chunk.can", "r");
        for (int j = 0; j < nc; j++) {
            int nd = 0, mn = -1;
            for (int v = 0; v < 10; v++) {
                if (!fgets(line, sizeof line, rf)) { fprintf(stderr, "short read\n"); return 1; }
                int len = strlen(line); while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0;
                int idx = lookup(line); if (idx < 0) { fprintf(stderr, "unknown deletion %s\n", line); return 1; }
                int seen = 0; for (int q = 0; q < nd; q++) if (cls[q] == idx) { seen = 1; break; }
                if (!seen) { cls[nd++] = idx; children[idx]++; if (mn < 0 || strcmp(can9[idx], can9[mn]) < 0) mn = idx; }
            }
            classes++; contacts += nd; hist[nd]++; branch[mn]++;
        }
        fclose(rf);
        fprintf(stderr, "processed %ld classes\n", classes);
    }
    pclose(s);
    long cmin = -1, cmax = 0, childless = 0, bmax = 0, bchildless = 0; double csum = 0;
    for (int i = 0; i < N9; i++) { if (cmin < 0 || children[i] < cmin) cmin = children[i]; if (children[i] > cmax) cmax = children[i]; csum += children[i]; if (!children[i]) childless++; if (branch[i] > bmax) bmax = branch[i]; if (!branch[i]) bchildless++; }
    printf("order 10: %ld classes; contacts %ld; mean distinct parents %.4f\n", classes, contacts, (double)contacts / classes);
    printf("distinct-parent histogram:"); for (int d = 1; d <= 10; d++) printf(" %d:%ld", d, hist[d]); printf("\n");
    printf("children per order-9 class: min %ld, mean %.3f, max %ld; transitive M_9 has %ld (%s); childless %ld\n", cmin, csum / N9, cmax, children[trans], children[trans] == cmax ? "max" : "not max", childless);
    printf("canonical-parent tree: max branching %ld (transitive %ld), childless order-9 classes %ld of %d (%.1f %%), mean %.3f\n", bmax, branch[trans], bchildless, N9, 100.0 * bchildless / N9, (double)classes / N9);
    return 0;
}
