/* line_census.c - opus-2026-11-09-S87.
 * The tiling line metagraph (THM-643 frame) aggregated by node invariants.
 * Tiling t: bit for every non-path arc (i, j), 0 <= i < j <= n-1, j >= i + 2; bit set means the arc is reversed (j -> i);
 * path arcs i -> i+1 are fixed.  tau(t) = t xor 1 (all tiles flipped); a LINE is a tau-orbit {t, tau t}.
 * sigma = grid transpose: tile (i, j) -> (n-1-j, n-1-i) (0-based), downstairs the converse; a line is BLUE iff it lies in
 * Fix(sigma) (sigma t = t; then sigma(tau t) = tau t since sigma fixes the all-ones vector), else BLACK.
 * Node invariants come from the class census (alpha_<n>.txt written by alpha_census.c, keyed by kind-pasteur's canonical
 * code): H, strong, self-converse, alpha_1..alpha_4, number of strong components.  The packing number
 * nu = max k with alpha_k > 0 (max number of vertex-disjoint odd cycles = of disjoint cyclic triangles).
 * Output: one record per merged class pair and colour: codeA codeB blue count   (codeA <= codeB), i.e. the merged line
 * multigraph, from which any classification of the nodes can be aggregated.
 * usage: line_census n alpha_n.txt SLICE NSLICES > pairs_n_SLICE.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 11
static int nn; static uint16_t A[MAXN]; static u64 best; static int have;
static u64 code_of(const int *pos){ int inv[MAXN]; for (int v=0; v<nn; v++) inv[pos[v]] = v; u64 c = 0; int b = 0;
  for (int i=0;i<nn;i++) for (int j=i+1;j<nn;j++){ if (A[inv[i]] >> inv[j] & 1) c |= 1ULL << b; b++; } return c; }
static int ncells(const int *cell){ int m=-1; for (int v=0; v<nn; v++) if (cell[v]>m) m=cell[v]; return m+1; }
static void refine(int *cell){ for (;;){ int k = ncells(cell); int sig[MAXN][MAXN+1];
    for (int v=0; v<nn; v++){ sig[v][0] = cell[v]; for (int c=0;c<k;c++) sig[v][c+1]=0; for (int u=0; u<nn; u++) if (A[v] >> u & 1) sig[v][cell[u]+1]++; }
    int ord[MAXN]; for (int v=0; v<nn; v++) ord[v]=v;
    for (int i=1;i<nn;i++){ int x=ord[i], j=i-1; while (j>=0 && memcmp(sig[ord[j]], sig[x], (k+1)*sizeof(int))>0){ ord[j+1]=ord[j]; j--; } ord[j+1]=x; }
    int nc[MAXN]; int c = 0; nc[ord[0]] = 0;
    for (int i=1;i<nn;i++){ if (memcmp(sig[ord[i]], sig[ord[i-1]], (k+1)*sizeof(int))) c++; nc[ord[i]] = c; }
    int changed = (c+1 != k); memcpy(cell, nc, sizeof(int)*nn); if (!changed) return; } }
static void search(int *cell){ refine(cell); int k = ncells(cell);
  if (k == nn){ u64 c = code_of(cell); if (!have || c < best){ best = c; have = 1; } return; }
  int cnt[MAXN]={0}; for (int v=0; v<nn; v++) cnt[cell[v]]++; int target=-1; for (int c=0;c<k;c++) if (cnt[c]>1){ target=c; break; }
  for (int v=0; v<nn; v++) if (cell[v]==target){ int c2[MAXN]; for (int u=0; u<nn; u++){ c2[u] = cell[u] > target ? cell[u]+1 : cell[u]; if (cell[u]==target && u!=v) c2[u] = target+1; } search(c2); } }
static u64 canon(const uint16_t *adj, int n){ nn = n; memcpy(A, adj, sizeof(uint16_t)*n); have = 0; int cell[MAXN]; for (int v=0; v<n; v++) cell[v]=0; search(cell); return best; }
/* hash map code -> index */
typedef struct { u64 *key; int *val; size_t cap; } hmap;
static void hm_put(hmap *h, u64 k, int v){ size_t i = (k*0x9E3779B97F4A7C15ULL) >> 17; i &= h->cap-1; while (h->key[i]!=~0ULL){ if (h->key[i]==k){ h->val[i]=v; return; } i=(i+1)&(h->cap-1); } h->key[i]=k; h->val[i]=v; }
static int hm_get(hmap *h, u64 k){ size_t i = (k*0x9E3779B97F4A7C15ULL) >> 17; i &= h->cap-1; while (h->key[i]!=~0ULL){ if (h->key[i]==k) return h->val[i]; i=(i+1)&(h->cap-1); } return -1; }
/* pair counter: key = (ia * NC + ib) * 2 + blue */
typedef struct { u64 *key; u64 *cnt; size_t cap, used; } pmap;
static void pm_add(pmap *h, u64 k){ size_t i = (k*0x9E3779B97F4A7C15ULL) >> 13; i &= h->cap-1; while (h->key[i]!=~0ULL){ if (h->key[i]==k){ h->cnt[i]++; return; } i=(i+1)&(h->cap-1); } h->key[i]=k; h->cnt[i]=1; h->used++; if (h->used*2 > h->cap){ fprintf(stderr,"pair map full\n"); exit(5); } }
int main(int argc, char **argv){
  int n = atoi(argv[1]); const char *af = argv[2]; int slice = atoi(argv[3]), ns = atoi(argv[4]);
  /* class table */
  hmap hm; hm.cap = 1u<<20; while (hm.cap < 4u*300000u) hm.cap <<= 1; hm.key = malloc(hm.cap*sizeof(u64)); hm.val = malloc(hm.cap*sizeof(int)); memset(hm.key,0xff,hm.cap*sizeof(u64));
  FILE *f = fopen(af, "r"); static char line[4096]; int NC = 0; u64 *codes = malloc(sizeof(u64)*300000);
  while (fgets(line, sizeof line, f)){ unsigned long long c; if (sscanf(line, "%llu", &c)==1){ codes[NC] = c; hm_put(&hm, c, NC); NC++; } }
  fclose(f);
  /* tiles */
  int ti[64], tj[64], m = 0; int idx[MAXN][MAXN];
  for (int i=0;i<n;i++) for (int j=i+2;j<n;j++){ ti[m]=i; tj[m]=j; idx[i][j]=m; m++; }
  int sg[64]; for (int b=0;b<m;b++){ int i2 = n-1-tj[b], j2 = n-1-ti[b]; sg[b] = idx[i2][j2]; }
  u64 all = (m==64) ? ~0ULL : ((1ULL<<m)-1);
  pmap pm; pm.cap = 1u<<22; pm.key = malloc(pm.cap*sizeof(u64)); pm.cnt = malloc(pm.cap*sizeof(u64)); memset(pm.key,0xff,pm.cap*sizeof(u64)); pm.used = 0;
  u64 half = 1ULL << (m-1); long nlines = 0, nblue = 0;
  for (u64 t = slice; t < half; t += ns){          /* lines: t with top bit 0, partner t ^ all */
    u64 ends[2] = { t, t ^ all }; int ci[2];
    for (int e=0;e<2;e++){ uint16_t adj[MAXN]; memset(adj,0,sizeof adj);
      for (int i=0;i+1<n;i++) adj[i] |= 1u<<(i+1);
      for (int b=0;b<m;b++){ int i=ti[b], j=tj[b]; if (ends[e]>>b & 1) adj[j] |= 1u<<i; else adj[i] |= 1u<<j; }
      u64 c = canon(adj, n); ci[e] = hm_get(&hm, c); if (ci[e] < 0){ fprintf(stderr, "unknown class\n"); return 3; } }
    /* blue: sigma t == t */
    u64 st = 0; for (int b=0;b<m;b++) if (t>>b & 1) st |= 1ULL << sg[b];
    int blue = (st == t);
    int a = ci[0], bb = ci[1]; if (a > bb){ int x=a; a=bb; bb=x; }
    pm_add(&pm, (((u64)a*NC + bb) << 1) | blue);
    nlines++; nblue += blue;
  }
  for (size_t i=0;i<pm.cap;i++) if (pm.key[i]!=~0ULL){ u64 k = pm.key[i]; int blue = k & 1; k >>= 1; int a = k / NC, b = k % NC;
    printf("%llu %llu %d %llu\n", codes[a], codes[b], blue, pm.cnt[i]); }
  fprintf(stderr, "n=%d slice %d/%d: lines %ld, blue %ld, merged pairs %zu\n", n, slice, ns, nlines, nblue, pm.used);
  return 0;
}
