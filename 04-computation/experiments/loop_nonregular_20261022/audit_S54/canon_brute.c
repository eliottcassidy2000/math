/* Independent auditor code: brute-force canonical form of a tournament restricted to score-respecting relabellings.
 * Every automorphism preserves scores, so with vertices grouped by score (ascending) and all permutations inside each
 * score class enumerated, canon = lexicographically minimal upper-triangle code over these relabellings, and
 * |Aut| = number of relabellings attaining it.  Also reports whether the converse tournament has the same canon
 * (self-converse) -- for the converse the score classes are reversed, which the same procedure handles.
 * Input lines "code ..."; output "canon aut selfconverse".  Usage: canon_brute n                                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, A[32][32], B[32][32];
static int order[32], nclass, cstart[33];
static int perm[32];                 /* perm[position] = vertex */
static char best[600], cur[600];
static long long cnt;
static int have;

static void eval(void) {
    int c = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) cur[c++] = B[perm[i]][perm[j]] ? '1' : '0';
    cur[c] = 0;
    if (!have) { strcpy(best, cur); cnt = 1; have = 1; return; }
    int r = strcmp(cur, best);
    if (r < 0) { strcpy(best, cur); cnt = 1; }
    else if (r == 0) cnt++;
}

/* enumerate permutations within each class by recursive swapping (Heap-like via recursion over positions) */
static void rec(int pos) {
    if (pos == n) { eval(); return; }
    /* class of position pos */
    int k = 0; while (!(cstart[k] <= pos && pos < cstart[k + 1])) k++;
    for (int i = pos; i < cstart[k + 1]; i++) {
        int t = perm[pos]; perm[pos] = perm[i]; perm[i] = t;
        rec(pos + 1);
        t = perm[pos]; perm[pos] = perm[i]; perm[i] = t;
    }
}

static void canon(char *out, long long *aut) {
    /* group vertices by score of B */
    int sc[32];
    for (int v = 0; v < n; v++) { sc[v] = 0; for (int w = 0; w < n; w++) sc[v] += B[v][w]; }
    int idx = 0; nclass = 0;
    for (int s = 0; s < n; s++) {
        int any = 0;
        for (int v = 0; v < n; v++) if (sc[v] == s) { if (!any) { cstart[nclass++] = idx; any = 1; } order[idx++] = v; }
    }
    cstart[nclass] = n;
    for (int i = 0; i < n; i++) perm[i] = order[i];
    have = 0; cnt = 0;
    rec(0);
    strcpy(out, best); *aut = cnt;
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    char buf[8192], code[1024], c1[600], c2[600];
    while (fgets(buf, sizeof buf, stdin)) {
        if (sscanf(buf, "%1023s", code) != 1) continue;
        int c = 0;
        for (int i = 0; i < n; i++) { A[i][i] = 0; for (int j = i + 1; j < n; j++) { A[i][j] = code[c++] == '1'; A[j][i] = !A[i][j]; } }
        long long a1, a2;
        memcpy(B, A, sizeof A); canon(c1, &a1);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) B[i][j] = A[j][i];
        canon(c2, &a2);
        printf("%s %lld %d\n", c1, a1, strcmp(c1, c2) == 0);
        fflush(stdout);
    }
    return 0;
}
