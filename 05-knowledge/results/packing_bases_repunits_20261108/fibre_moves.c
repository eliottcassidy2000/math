/* fibre_moves.c - opus-2026-11-08-S86.
 * H-fibres under H-preserving elementary moves. Nodes: all isomorphism classes of tournaments of orders 1..NMAX
 * (from kind-pasteur's classgen.c output files cls_<n>.txt). Edges (undirected):
 *   (D) T -- T - v      when H(T - v) = H(T)        (H-preserving vertex deletion / insertion)
 *   (F) T -- T^(uv)     when H(T^(uv)) = H(T)       (H-preserving single-arc reversal, same order)
 * Question: does every class lie in the same component as a first creator of its value (a class of least order with
 * that H)?  Equivalently: are the first creators of h enough to generate every tournament with H = h by H-preserving
 * single-vertex and single-arc moves?
 * Canonical forms: classgen.c's individualisation-refinement (kind-pasteur-2026-11-07-S16), reproduced.
 * usage: fibre_moves NMAX dir   (reads dir/cls_<n>.txt)
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
static void decode(u64 c, int n, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*n); int b=0;
  for (int i=0;i<n;i++) for (int j=i+1;j<n;j++){ if (c>>b & 1) adj[i] |= 1u<<j; else adj[j] |= 1u<<i; b++; } }
/* per-order sorted arrays of (code, H) for lookup */
typedef struct { u64 code; u64 H; int id; } rec;
static rec *R[MAXN+1]; static size_t NR[MAXN+1]; static int base[MAXN+2];
static int cmp(const void *a, const void *b){ u64 x = ((const rec*)a)->code, y = ((const rec*)b)->code; return x < y ? -1 : x > y; }
static rec *find(int n, u64 code){ size_t lo = 0, hi = NR[n]; while (lo < hi){ size_t m = (lo+hi)/2; if (R[n][m].code < code) lo = m+1; else hi = m; } return (lo < NR[n] && R[n][lo].code == code) ? &R[n][lo] : NULL; }
static int *par;
static int fnd(int x){ while (par[x] != x){ par[x] = par[par[x]]; x = par[x]; } return x; }
static void uni(int a, int b){ a = fnd(a); b = fnd(b); if (a != b) par[a] = b; }
int main(int argc, char **argv){
  int N = atoi(argv[1]); const char *dir = argv[2];
  int total = 0; static char line[1<<16];
  /* first pass: records */
  for (int n=1;n<=N;n++){ char fn[1024]; snprintf(fn, sizeof fn, "%s/cls_%d.txt", dir, n); FILE *f = fopen(fn,"r"); if (!f){ fprintf(stderr,"missing %s\n", fn); return 1; }
    size_t cap = 1024; R[n] = malloc(cap*sizeof(rec)); NR[n] = 0;
    while (fgets(line, sizeof line, f)){ u64 code, H; if (sscanf(line, "%llu %llu", &code, &H) != 2) continue; if (NR[n]==cap){ cap*=2; R[n]=realloc(R[n],cap*sizeof(rec)); } R[n][NR[n]].code = code; R[n][NR[n]].H = H; NR[n]++; }
    fclose(f); qsort(R[n], NR[n], sizeof(rec), cmp); base[n] = total; for (size_t i=0;i<NR[n];i++) R[n][i].id = total + (int)i; total += (int)NR[n]; }
  par = malloc(total*sizeof(int)); for (int i=0;i<total;i++) par[i]=i;
  long dedges = 0, fedges = 0;
  for (int n=1;n<=N;n++){
    for (size_t r=0;r<NR[n];r++){ u64 code = R[n][r].code, H = R[n][r].H; int id = R[n][r].id; uint16_t adj[MAXN]; decode(code, n, adj);
      /* deletions */
      if (n > 1) for (int x=0;x<n;x++){ uint16_t d[MAXN]; int k=0; for (int v=0;v<n;v++) if (v!=x){ uint16_t m=0; int kk=0; for (int u=0;u<n;u++) if (u!=x){ if (adj[v]>>u&1) m|=1u<<kk; kk++; } d[k++]=m; }
        rec *q = find(n-1, canon(d, n-1)); if (!q){ fprintf(stderr,"lookup fail\n"); return 2; } if (q->H == H){ uni(id, q->id); dedges++; } }
      /* single-arc reversals */
      for (int i=0;i<n;i++) for (int j=i+1;j<n;j++){ uint16_t a2[MAXN]; memcpy(a2, adj, sizeof(uint16_t)*n);
        if (a2[i]>>j&1){ a2[i]&=~(1u<<j); a2[j]|=1u<<i; } else { a2[j]&=~(1u<<i); a2[i]|=1u<<j; }
        rec *q = find(n, canon(a2, n)); if (!q){ fprintf(stderr,"lookup fail flip\n"); return 2; } if (q->H == H){ uni(id, q->id); fedges++; } }
    }
    fprintf(stderr, "order %d done (%zu classes)\n", n, NR[n]);
  }
  /* first-creator order per value */
  #define HMAX 400000
  static int firstn[HMAX]; for (int h=0;h<HMAX;h++) firstn[h]=0;
  for (int n=1;n<=N;n++) for (size_t r=0;r<NR[n];r++){ u64 H = R[n][r].H; if (H < HMAX && !firstn[H]) firstn[H] = n; }
  /* components: does each contain a first creator of its H? */
  static char hasfc[1<<22]; memset(hasfc, 0, sizeof hasfc);
  for (int n=1;n<=N;n++) for (size_t r=0;r<NR[n];r++) if (firstn[R[n][r].H] == n) hasfc[fnd(R[n][r].id)] = 1;
  long orphans = 0; long per_order[MAXN+1] = {0};
  for (int n=1;n<=N;n++) for (size_t r=0;r<NR[n];r++) if (!hasfc[fnd(R[n][r].id)]){ orphans++; per_order[n]++; }
  /* number of components per value */
  printf("nodes %d, H-preserving deletion edges %ld, H-preserving flip edges %ld\n", total, dedges, fedges);
  printf("classes NOT connected to a first creator of their value: %ld\n", orphans);
  for (int n=1;n<=N;n++) printf("  order %d: %ld of %zu\n", n, per_order[n], NR[n]);
  /* list values whose fibre (orders <= N) is disconnected */
  int *compcount = calloc(HMAX, sizeof(int));
  static char seen[1<<22];
  for (int n=1;n<=N;n++) for (size_t r=0;r<NR[n];r++){ int c = fnd(R[n][r].id); if (!seen[c]){ seen[c] = 1; u64 H = R[n][r].H; if (H < HMAX) compcount[H]++; } }
  int printed = 0; printf("values with more than one H-preserving component among orders <= %d (value: components):", N);
  for (int h=1; h<HMAX; h+=2) if (compcount[h] > 1){ if (printed < 60) printf(" %d:%d", h, compcount[h]); printed++; }
  printf("\n  (%d values)\n", printed);
  if (argc > 3){ u64 H0 = strtoull(argv[3], 0, 10);              /* optional: list the components of one value */
    printf("classes with H = %llu (order, canonical code, component representative):\n", H0);
    for (int n=1;n<=N;n++) for (size_t r=0;r<NR[n];r++) if (R[n][r].H == H0) printf("  %d %llu comp %d\n", n, R[n][r].code, fnd(R[n][r].id)); }
  return 0;
}
