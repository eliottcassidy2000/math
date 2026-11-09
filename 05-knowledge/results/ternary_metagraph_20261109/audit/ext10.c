/* ext10.c - independent audit of THM-4950 (4): the set of H_3 values of the STRONG tournaments of order 10.
 * Written from scratch by the auditor; no canonical forms are needed.  Every tournament of order 10 is, after deleting
 * vertex 9, isomorphic to one of the 191,536 order-9 classes of gen_classes.c, so it is isomorphic to a one-vertex
 * extension of a class representative.  For each representative R and each out-set O of the new vertex 9:
 *   - cycles through 9 with vertex set S + {9}: sum over u in S & O, w in S \ O of hp_R[S](u -> w)  (u .. w is a
 *     Hamiltonian path of R[S], closed by w -> 9 -> u);
 *   - Z_x(T) = Z_x(R) + x * sum_S cyc9(S) Z_x(R[V \ S])     (independence polynomial of the odd-cycle graph at x),
 *     where Z_x(R[U]) is precomputed for all U by the recursion on the least vertex.
 * H_3 = Z_3(T), and H = Z_2(T) is checked against an independent order-10 Hamiltonian-path DP on every 2048th strong
 * extension, together with a direct order-10 recomputation of H_3.
 * usage: ext10 acls_9.txt out_prefix        (writes out_prefix.bin: byte v = 1 iff some strong order-10 T has H_3 = v,
 *                                             and out_prefix_wit.txt: 'v base_code O' for the least (index, O) hit)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef unsigned long long u64;
typedef long long i64;
#define VMAX (1 << 17)

static void decode(u64 code, int n, uint32_t *out) {
  int b = 0;
  for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++, b++) {
      if ((code >> b) & 1ULL) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}

/* direct order-n recomputation (n = 10) for the sample checks */
static i64 direct_h3(const uint32_t *out, int n, i64 *H2out, i64 *hpout) {
  int F = 1 << n;
  i64 *g = calloc((size_t)F * n, sizeof(i64));
  i64 *cyc = calloc(F, sizeof(i64));
  for (int m = 0; m < n; m++) {
    memset(g, 0, (size_t)F * n * sizeof(i64));
    g[(1 << m) * n + m] = 1;
    for (int S = 1 << m; S < F; S++) {
      if (!((S >> m) & 1) || (S & ((1 << m) - 1))) continue;
      for (int v = m; v < n; v++) {
        i64 c = g[S * n + v];
        if (!c) continue;
        uint32_t nb = out[v] & ~(uint32_t)S & ~((1u << m) - 1);
        while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; g[(S | (1 << u)) * n + u] += c; }
      }
      if (__builtin_popcount(S) >= 3 && (__builtin_popcount(S) & 1)) {
        i64 t = 0;
        for (int v = m + 1; v < n; v++) if ((out[v] >> m) & 1u) t += g[S * n + v];
        cyc[S] = t;
      }
    }
  }
  i64 *z3 = calloc(F, sizeof(i64)), *z2 = calloc(F, sizeof(i64));
  z3[0] = z2[0] = 1;
  for (int U = 1; U < F; U++) {
    int m = __builtin_ctz(U), rest = U & ~(1 << m);
    z3[U] = z3[rest]; z2[U] = z2[rest];
    for (int T = rest; ; T = (T - 1) & rest) {
      int S = T | (1 << m);
      if (cyc[S]) { z3[U] += 3 * cyc[S] * z3[U & ~S]; z2[U] += 2 * cyc[S] * z2[U & ~S]; }
      if (!T) break;
    }
  }
  i64 h3 = z3[F - 1];
  *H2out = z2[F - 1];
  /* Hamiltonian paths */
  memset(g, 0, (size_t)F * n * sizeof(i64));
  for (int v = 0; v < n; v++) g[(1 << v) * n + v] = 1;
  for (int S = 1; S < F; S++)
    for (int v = 0; v < n; v++) {
      i64 c = g[S * n + v];
      if (!c) continue;
      uint32_t nb = out[v] & ~(uint32_t)S;
      while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; g[(S | (1 << u)) * n + u] += c; }
    }
  i64 hp = 0;
  for (int v = 0; v < n; v++) hp += g[(F - 1) * n + v];
  *hpout = hp;
  free(g); free(cyc); free(z3); free(z2);
  return h3;
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: ext10 acls_9.txt out_prefix\n"); return 1; }
  FILE *fi = fopen(argv[1], "r");
  size_t cap = 1 << 18, nr = 0;
  u64 *reps = malloc(cap * sizeof(u64));
  char line[512];
  while (fgets(line, sizeof line, fi)) { u64 c; if (sscanf(line, "%llu", &c) == 1) { if (nr == cap) { cap *= 2; reps = realloc(reps, cap * sizeof(u64)); } reps[nr++] = c; } }
  fclose(fi);
  fprintf(stderr, "read %zu order-9 representatives\n", nr);
  static unsigned char seen[VMAX];
  static u64 wit[VMAX];
  for (int v = 0; v < VMAX; v++) wit[v] = ~0ULL;
  long long nstrong = 0, nchecks = 0, nbad = 0;
  i64 maxh3 = 0, minh3 = (i64)1 << 60;
  const int n = 9, F = 1 << 9, FULL = F - 1;
  int evenS[256]; int ne = 0;
  for (int S = 1; S < F; S++) if (!(__builtin_popcount(S) & 1)) evenS[ne++] = S;
#pragma omp parallel reduction(+:nstrong,nchecks,nbad) reduction(max:maxh3) reduction(min:minh3)
  {
    i64 (*hp)[9][9] = malloc(sizeof(i64) * F * 81);
    i64 (*f)[9] = malloc(sizeof(i64) * F * 9);
    i64 *cyc = malloc(sizeof(i64) * F), *z3 = malloc(sizeof(i64) * F), *z2 = malloc(sizeof(i64) * F);
    unsigned char *lseen = calloc(VMAX, 1);
    u64 *lwit = malloc(sizeof(u64) * VMAX);
    for (int v = 0; v < VMAX; v++) lwit[v] = ~0ULL;
#pragma omp for schedule(dynamic, 64)
    for (long r = 0; r < (long)nr; r++) {
      uint32_t out[10], in[10];
      decode(reps[r], 9, out);
      for (int v = 0; v < 9; v++) in[v] = 0;
      for (int v = 0; v < 9; v++) for (int u = 0; u < 9; u++) if ((out[v] >> u) & 1u) in[u] |= 1u << v;
      /* Hamiltonian path counts hp[S][u][w] of R[S] */
      for (int u = 0; u < 9; u++) {
        memset(f, 0, sizeof(i64) * F * 9);
        f[1 << u][u] = 1;
        for (int S = 1; S < F; S++) {
          if (!((S >> u) & 1)) continue;
          for (int w = 0; w < 9; w++) {
            i64 c = f[S][w];
            if (!c) continue;
            uint32_t nb = out[w] & ~(uint32_t)S;
            while (nb) { int x = __builtin_ctz(nb); nb &= nb - 1; f[S | (1 << x)][x] += c; }
          }
        }
        for (int S = 0; S < F; S++) for (int w = 0; w < 9; w++) hp[S][u][w] = ((S >> u) & 1) ? f[S][w] : 0;
      }
      for (int S = 0; S < F; S++) {
        cyc[S] = 0;
        int sz = __builtin_popcount(S);
        if (sz >= 3 && (sz & 1)) {
          int m = __builtin_ctz(S);
          for (int w = 0; w < 9; w++) if (((S >> w) & 1) && ((out[w] >> m) & 1u)) cyc[S] += hp[S][m][w];
        }
      }
      z3[0] = z2[0] = 1;
      for (int U = 1; U < F; U++) {
        int m = __builtin_ctz(U), rest = U & ~(1 << m);
        z3[U] = z3[rest]; z2[U] = z2[rest];
        for (int T = rest; ; T = (T - 1) & rest) {
          int S = T | (1 << m);
          if (cyc[S]) { z3[U] += 3 * cyc[S] * z3[U & ~S]; z2[U] += 2 * cyc[S] * z2[U & ~S]; }
          if (!T) break;
        }
      }
      for (int O = 1; O < FULL; O++) {
        int I = FULL & ~O;
        /* strong: everything reachable from 9 (via O) and 9 reachable from everything (via I) */
        uint32_t a = (uint32_t)O, last = 0;
        while (a != last) { last = a; for (int v = 0; v < 9; v++) if ((a >> v) & 1u) a |= out[v]; }
        if (a != (uint32_t)FULL) continue;
        uint32_t b = (uint32_t)I; last = 0;
        while (b != last) { last = b; for (int v = 0; v < 9; v++) if ((b >> v) & 1u) b |= in[v]; }
        if (b != (uint32_t)FULL) continue;
        i64 h3 = z3[FULL], h2 = z2[FULL];
        for (int e = 0; e < ne; e++) {
          int S = evenS[e];
          int SO = S & O, SI = S & I;
          if (!SO || !SI) continue;
          i64 c = 0;
          for (int uu = SO; uu; uu &= uu - 1) {
            int u = __builtin_ctz(uu);
            for (int ww = SI; ww; ww &= ww - 1) c += hp[S][u][__builtin_ctz(ww)];
          }
          if (c) { h3 += 3 * c * z3[FULL & ~S]; h2 += 2 * c * z2[FULL & ~S]; }
        }
        nstrong++;
        if (h3 >= VMAX || h3 < 1) { fprintf(stderr, "value out of range %lld\n", h3); exit(2); }
        lseen[h3] = 1;
        u64 key = ((u64)r << 9) | (u64)O;
        if (key < lwit[h3]) lwit[h3] = key;
        if (h3 > maxh3) maxh3 = h3;
        if (h3 < minh3) minh3 = h3;
        if ((nstrong & 2047) == 0) {
          uint32_t t[10];
          for (int v = 0; v < 9; v++) t[v] = out[v];
          t[9] = 0;
          for (int v = 0; v < 9; v++) { if ((O >> v) & 1) t[9] |= 1u << v; else t[v] |= 1u << 9; }
          i64 H2d, hpd;
          i64 h3d = direct_h3(t, 10, &H2d, &hpd);
          nchecks++;
          if (h3d != h3 || H2d != h2 || hpd != h2) nbad++;
        }
      }
    }
#pragma omp critical
    {
      for (int v = 0; v < VMAX; v++) { if (lseen[v]) seen[v] = 1; if (lwit[v] < wit[v]) wit[v] = lwit[v]; }
    }
    free(hp); free(f); free(cyc); free(z3); free(z2); free(lseen); free(lwit);
  }
  long distinct = 0;
  for (int v = 0; v < VMAX; v++) distinct += seen[v];
  char fn[1024];
  snprintf(fn, sizeof fn, "%s.bin", argv[2]);
  FILE *fo = fopen(fn, "wb"); fwrite(seen, 1, VMAX, fo); fclose(fo);
  snprintf(fn, sizeof fn, "%s_wit.txt", argv[2]);
  fo = fopen(fn, "w");
  for (int v = 0; v < VMAX; v++) if (seen[v]) fprintf(fo, "%d %llu %llu\n", v, reps[wit[v] >> 9], wit[v] & 511ULL);
  fclose(fo);
  printf("order-10 strong one-vertex extensions processed (with repeats): %lld\n", nstrong);
  printf("distinct strong order-10 H_3 values: %ld; min %lld; max %lld\n", distinct, minh3, maxh3);
  printf("sample checks against a direct order-10 recomputation (H_3, OCF H, Hamiltonian paths): %lld, failures %lld\n", nchecks, nbad);
  return 0;
}
