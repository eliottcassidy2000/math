/* order10.c (derived from classgen.c) - kind-pasteur-2026-11-07-S16. Isomorphism classes of tournaments of order 1..NMAX by vertex extension
 * and canonical forms (individualization-refinement with ordered colour refinement), with, for every class:
 *   code (canonical upper-triangle bit code), H (Hamiltonian paths), c3, strong, SC (self-converse), aut (|Aut|),
 *   and the canonical codes and H values of the n vertex-deleted subtournaments.
 * usage: classgen NMAX outprefix      writes outprefix_<n>.txt for n = 1..NMAX
 * line: code H c3 strong SC aut ccode | del_code_1:del_H_1 ... (ccode = canonical code of the converse)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 12
static int N;
/* ---------- canonical form ---------- */
typedef struct { int cell[MAXN]; } part;
static int nn;                 /* current order */
static uint16_t A[MAXN];       /* out-neighbour masks of the tournament being canonicalised */
static u64 best; static int have; static long leaves_equal;
static u64 code_of(const int *pos){ /* pos[v] = canonical position of v */
  int inv[MAXN]; for (int v=0; v<nn; v++) inv[pos[v]] = v;
  u64 c = 0; int b = 0;
  for (int i=0; i<nn; i++) for (int j=i+1; j<nn; j++){ if (A[inv[i]] >> inv[j] & 1) c |= 1ULL << b; b++; }
  return c;
}
static int ncells(const int *cell){ int m=-1; for (int v=0; v<nn; v++) if (cell[v]>m) m=cell[v]; return m+1; }
/* ordered colour refinement: new colour = (old colour, out-counts per old colour), ordered lexicographically */
static void refine(int *cell){
  for (;;){
    int k = ncells(cell);
    int sig[MAXN][MAXN+1];
    for (int v=0; v<nn; v++){ sig[v][0] = cell[v]; for (int c=0;c<k;c++) sig[v][c+1]=0;
      for (int u=0; u<nn; u++) if (A[v] >> u & 1) sig[v][cell[u]+1]++; }
    int ord[MAXN]; for (int v=0; v<nn; v++) ord[v]=v;
    for (int i=1;i<nn;i++){ int x=ord[i], j=i-1; while (j>=0 && memcmp(sig[ord[j]], sig[x], (k+1)*sizeof(int))>0){ ord[j+1]=ord[j]; j--; } ord[j+1]=x; }
    int nc[MAXN]; int c = 0; nc[ord[0]] = 0;
    for (int i=1;i<nn;i++){ if (memcmp(sig[ord[i]], sig[ord[i-1]], (k+1)*sizeof(int))) c++; nc[ord[i]] = c; }
    int changed = (c+1 != k);
    memcpy(cell, nc, sizeof(int)*nn);
    if (!changed) return;
  }
}
static void search(int *cell){
  refine(cell);
  int k = ncells(cell);
  if (k == nn){ u64 c = code_of(cell); if (!have || c < best){ best = c; have = 1; leaves_equal = 1; } else if (c == best) leaves_equal++; return; }
  /* first non-singleton cell */
  int cnt[MAXN]={0}; for (int v=0; v<nn; v++) cnt[cell[v]]++;
  int target=-1; for (int c=0;c<k;c++) if (cnt[c]>1){ target=c; break; }
  for (int v=0; v<nn; v++) if (cell[v]==target){
    int c2[MAXN]; for (int u=0; u<nn; u++){ c2[u] = cell[u] > target ? cell[u]+1 : cell[u]; if (cell[u]==target && u!=v) c2[u] = target+1; }
    search(c2);
  }
}
static u64 canon(const uint16_t *adj, int n, long *aut){
  nn = n; memcpy(A, adj, sizeof(uint16_t)*n); have = 0;
  int cell[MAXN]; for (int v=0; v<n; v++) cell[v]=0;
  search(cell);
  if (aut) *aut = leaves_equal;   /* number of leaves with the minimal code = |Aut| (no pruning) */
  return best;
}
static void decode(u64 c, int n, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*n); int b=0;
  for (int i=0;i<n;i++) for (int j=i+1;j<n;j++){ if (c>>b & 1) adj[i] |= 1u<<j; else adj[j] |= 1u<<i; b++; } }
/* ---------- invariants ---------- */
static u64 hpaths(const uint16_t *adj, int n){
  static u64 f[1<<MAXN][MAXN]; int F = 1<<n;
  for (int S=0;S<F;S++) for (int v=0;v<n;v++) f[S][v]=0;
  for (int v=0;v<n;v++) f[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<n;v++) if (f[S][v]){ u64 x=f[S][v]; uint16_t out = adj[v] & ~S; while(out){ int u=__builtin_ctz(out); out&=out-1; f[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<n;v++) h+=f[F-1][v]; return h;
}
static int c3count(const uint16_t *adj, int n){ int s=n*(n-1)*(n-2)/6; for(int v=0;v<n;v++){ int d=__builtin_popcount(adj[v]); s-=d*(d-1)/2; } return s; }
static int strongq(const uint16_t *adj, int n){ if(n<=1) return 1; uint16_t full=(1u<<n)-1, r=1, prev=0;
  while(r!=prev){ prev=r; for(int v=0;v<n;v++) if(r>>v&1) r|=adj[v]; } if(r!=full) return 0;
  uint16_t in[MAXN]={0}; for(int v=0;v<n;v++) for(int u=0;u<n;u++) if(adj[v]>>u&1) in[u]|=1u<<v;
  r=1; prev=0; while(r!=prev){ prev=r; for(int v=0;v<n;v++) if(r>>v&1) r|=in[v]; } return r==full; }


int main(int argc, char **argv){
  /* reads lines "W code10 x H status" ; prints canonical code of the 11-tournament, H, status */
  FILE *fp = fopen(argv[1], "r"); char line[512];
  while (fgets(line, sizeof line, fp)){
    char tag[8], st[16]; unsigned long long code; unsigned x; int h;
    if (sscanf(line, "%7s %llu %u %d %15s", tag, &code, &x, &h, st) != 5 || tag[0] != 'W') continue;
    uint16_t adj[MAXN]; decode(code, 10, adj); adj[10] = 0;
    for (int v=0; v<10; v++){ if (x>>v & 1) adj[v] |= 1u<<10; else adj[10] |= 1u<<v; }
    long aut; u64 c = canon(adj, 11, &aut);
    printf("%llu %d %s %ld\n", c, h, st, aut);
  }
  return 0;
}
