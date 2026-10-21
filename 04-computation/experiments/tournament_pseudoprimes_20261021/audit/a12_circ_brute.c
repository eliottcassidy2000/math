/* AUDIT a12 (independent auditor, 2026-10-21): brute force over ALL permutations for every connection set S of every
 * circulant tournament Cay(Z_n,S), n = 3,5,7,9 (and 11 optional, 39.9M perms per set): H, C (directed Ham cycles:
 * permutations with v_1 = 0 closing up), rho0 (paths with v_{n+1-i} = -v_i).  Prints one line per S.
 * Usage: ./a12_circ_brute n */
#include <stdio.h>
#include <stdlib.h>
static int n, inS[32], perm[32];
static int arc(int x, int y) { return inS[((y - x) % n + n) % n]; }
static int next_perm(int *a, int k) { int i = k - 2; while (i >= 0 && a[i] >= a[i + 1]) i--; if (i < 0) return 0;
    int j = k - 1; while (a[j] <= a[i]) j--; int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = k - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; } return 1; }
int main(int argc, char **argv) {
    n = atoi(argv[1]); int m = (n - 1) / 2;
    for (unsigned code = 0; code < (1u << m); code++) {
        for (int s = 0; s < n; s++) inS[s] = 0;
        for (int j = 0; j < m; j++) inS[(code >> j & 1) ? n - (j + 1) : j + 1] = 1;
        for (int i = 0; i < n; i++) perm[i] = i;
        unsigned long long H = 0, C = 0, R = 0;
        do { int ok = 1; for (int i = 0; i + 1 < n; i++) if (!arc(perm[i], perm[i + 1])) { ok = 0; break; }
             if (!ok) continue; H++;
             if (perm[0] == 0 && arc(perm[n - 1], perm[0])) C++;
             int sym = 1; for (int i = 0; i < n; i++) if (perm[n - 1 - i] != (n - perm[i]) % n) { sym = 0; break; }
             R += sym;
        } while (next_perm(perm, n));
        printf("n=%d S={", n); int f = 1; for (int s = 1; s < n; s++) if (inS[s]) { printf(f ? "%d" : ",%d", s); f = 0; }
        printf("} H=%llu C=%llu rho0=%llu\n", H, C, R);
    }
    return 0; }
