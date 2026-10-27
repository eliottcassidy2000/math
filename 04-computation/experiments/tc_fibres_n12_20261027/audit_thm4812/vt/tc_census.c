// tc_census.c -- independent census of tc = H/|Aut| over the strong tournament classes of one order
// (input: `gentourng -c -q n`), histogram of tc values <= VMAX.  Own H (pull Held-Karp) and |Aut| (tourlib).
#include "../common/tourlib.h"
int main(int argc, char **argv) {
    int n = atoi(argv[1]); long long vmax = atoll(argv[2]);
    uint64_t *F = malloc(sizeof(uint64_t) * ((size_t)1 << n) * n);
    static long long hist[100000]; long long classes = 0, nonstrong = 0, divfail = 0; char line[512];
    long long minTc = -1;
    while (fgets(line, sizeof line, stdin)) {
        if (n > 1 && line[0] != '0' && line[0] != '1') continue;
        uint32_t out[32] = {0}; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; k++; }
        if (!strong_landau(n, out)) { nonstrong++; continue; }
        classes++;
        uint32_t in[32]; make_in(n, out, in);
        uint64_t H = ham_pull(n, in, F); long long A = aut_order(n, out);
        if (H % (uint64_t)A) { divfail++; continue; }
        long long tc = (long long)(H / (uint64_t)A);
        if (minTc < 0 || tc < minTc) minTc = tc;
        if (tc <= vmax) hist[tc]++;
    }
    printf("order %d: strong classes %lld (non-strong in input %lld, H%%Aut failures %lld), min tc %lld; tc <= %lld:", n, classes, nonstrong, divfail, minTc, vmax);
    for (long long t = 1; t <= vmax; t++) if (hist[t]) printf(" %lld:%lld", t, hist[t]);
    printf("\n");
    return 0;
}
