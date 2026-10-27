/* w4free.c -- count tournament classes in which every out-neighbourhood (W4-free) / in-neighbourhood (L4-free) /
 * both (locally transitive) is transitive, i.e. contains no 3-cycle.  Reads gentourng.  cc -O3 w4free.c -o w4free */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
int n; uint32_t out_[32];
static int has_c3(uint32_t S) { for (uint32_t a = S; a; a &= a - 1) { int u = __builtin_ctz(a); uint32_t U = out_[u] & S;
    for (uint32_t b = U; b; b &= b - 1) { int v = __builtin_ctz(b); if (out_[v] & S & ~U & ~(1u << u)) { uint32_t W = out_[v] & S; for (uint32_t c = W; c; c &= c - 1) { int w = __builtin_ctz(c); if (out_[w] >> u & 1) return 1; } } } } return 0; }
int main(int argc, char **argv) { n = atoi(argv[1]); char line[256]; long cls = 0, w4f = 0, lt = 0;
    uint32_t full = (1u << n) - 1;
    while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; cls++;
        for (int i = 0; i < n; i++) out_[i] = 0; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i; k++; }
        int okO = 1, okI = 1;
        for (int v = 0; v < n && okO; v++) if (has_c3(out_[v])) okO = 0;
        if (okO) { w4f++; for (int v = 0; v < n && okI; v++) { uint32_t in = full & ~out_[v] & ~(1u << v); if (has_c3(in)) okI = 0; } if (okI) lt++; } }
    printf("n=%d classes %ld: W4-free (every out-neighbourhood transitive) %ld; locally transitive (both) %ld\n", n, cls, w4f, lt); return 0; }
