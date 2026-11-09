/* ternary_order10.c - opus-2026-11-09-S87.
 * Base-3 odd-cycle packing counts H_3(T) = I(Omega(T), 3) of every STRONG tournament class of order 10.
 * Generation: one-vertex extensions of the order-9 classes (cls_9.txt from kind-pasteur's classgen.c), canonical form
 * as in classgen.c / order10.c (kind-pasteur-2026-11-07-S16), deduplicated per slice.  Every order-10 class arises
 * from some order-9 class by adding a vertex, so the union over all slices covers every class; duplicates across
 * slices only repeat a value, which is harmless because only the SET of attained H_3 values is recorded.
 * For each new strong class: odd cycles by subset DP, the packing vector alpha_k by subset DP, H_3 = sum alpha_k 3^k,
 * and the check sum alpha_k 2^k = H (Hamiltonian paths, independent DP).
 * usage: ternary_order10 cls_9.txt SLICE NSLICES out_prefix   (writes out_prefix_SLICE.bin: one byte per value < 2^17)
 * Witness mode: ternary_order10 cls_9.txt SLICE NSLICES out_prefix target.bin  re-runs the same deterministic slice,
 * records for each value the canonical code (bit b of the 45-bit code, pairs i < j in order, set iff i -> j) of the
 * first class attaining it, stops once every value of target.bin (that slice's own bitmap) has a witness, and writes
 * out_prefix_SLICE.wit (lines 'H_3 code').  verify_witness10.py re-checks the witnesses independently.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 11
#define VMAX (1 << 17)
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
static int strongq(const uint16_t *adj, int n){ uint16_t full=(1u<<n)-1, r=1, prev=0;
  while(r!=prev){ prev=r; for(int v=0;v<n;v++) if(r>>v&1) r|=adj[v]; } if(r!=full) return 0;
  r=1; prev=0; while(r!=prev){ prev=r; for(int v=0;v<n;v++) if(!(r>>v&1) && (adj[v]&r)) r|=1u<<v; } return r==full; }
static u64 cyc[1<<10], dpc[1<<10][10], P[1<<10][5];
static void packing(const uint16_t *adj, int n, u64 *alpha){
  int F = 1 << n;
  for (int S=0;S<F;S++){ cyc[S]=0; for (int v=0;v<n;v++) dpc[S][v]=0; for (int k=0;k<5;k++) P[S][k]=0; }
  for (int s=0;s<n;s++) dpc[1<<s][s]=1;
  for (int S=1;S<F;S++){ int s=__builtin_ctz(S);
    for (int v=0;v<n;v++){ u64 c=dpc[S][v]; if(!c) continue; uint16_t nb = adj[v] & ~S & ~((1u<<(s+1))-1);
      while(nb){ int u=__builtin_ctz(nb); nb&=nb-1; dpc[S|1<<u][u]+=c; } }
    int k=__builtin_popcount(S); if (k>=3 && (k&1)){ u64 t=0; for (int v=0;v<n;v++) if (dpc[S][v] && (adj[v]>>s&1)) t+=dpc[S][v]; cyc[S]=t; } }
  P[0][0]=1;
  for (int U=1;U<F;U++){ int m=__builtin_ctz(U); int rest=U&~(1<<m);
    for (int T=rest;;T=(T-1)&rest){ int S=T|(1<<m); if (cyc[S]) for (int k=0;k<4;k++) if (P[U&~S][k]) P[U][k+1]+=cyc[S]*P[U&~S][k]; if(!T) break; } }
  for (int k=0;k<5;k++){ alpha[k]=0; for (int U=0;U<F;U++) alpha[k]+=P[U][k]; }
}
static u64 hpaths(const uint16_t *adj, int n){ static u64 f[1<<10][10]; int F=1<<n;
  for (int S=0;S<F;S++) for (int v=0;v<n;v++) f[S][v]=0; for (int v=0;v<n;v++) f[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<n;v++) if (f[S][v]){ u64 x=f[S][v]; uint16_t o=adj[v]&~S; while(o){ int u=__builtin_ctz(o); o&=o-1; f[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<n;v++) h+=f[F-1][v]; return h; }
typedef struct { u64 *key; size_t cap, cnt; } hset;
static int hs_insert(hset *h, u64 k){ size_t i = (k*0x9E3779B97F4A7C15ULL) >> 20; i &= h->cap-1;
  while (h->key[i]!=~0ULL){ if (h->key[i]==k) return 0; i=(i+1)&(h->cap-1); } h->key[i]=k; h->cnt++; return 1; }
int main(int argc, char **argv){
  if (argc < 5){ fprintf(stderr, "usage: ternary_order10 cls_9.txt SLICE NSLICES out_prefix\n"); return 1; }
  int slice = atoi(argv[2]), ns = atoi(argv[3]);
  FILE *f9 = fopen(argv[1], "r"); static char line[1<<16]; size_t nr=0, cap=0; u64 *reps=NULL;
  while (fgets(line, sizeof line, f9)){ unsigned long long c; if (sscanf(line, "%llu", &c)==1){ if(nr==cap){cap=cap?2*cap:1024; reps=realloc(reps,cap*sizeof(u64));} reps[nr++]=c; } }
  fclose(f9);
  hset hs; hs.cap = 1ULL<<23; hs.key = malloc(hs.cap*sizeof(u64)); memset(hs.key,0xff,hs.cap*sizeof(u64)); hs.cnt=0;
  static unsigned char seen[VMAX]; long strong_new=0, checks=0, bad=0, alpha3pos=0; u64 maxh3=0;
  int wmode = argc > 5; static unsigned char target[VMAX]; static u64 wit[VMAX]; long need = 0, got = 0;
  if (wmode){ FILE *ft = fopen(argv[5], "rb"); if (!ft || fread(target, 1, VMAX, ft) != VMAX){ fprintf(stderr, "bad target\n"); return 4; } fclose(ft);
    for (int v=0; v<VMAX; v++) need += target[v] != 0; }
  for (size_t r = slice; r < nr; r += ns){ uint16_t adj[MAXN]; decode(reps[r], 9, adj);
    for (int Pt=1; Pt<(1<<9)-1; Pt++){                         /* strong extensions need an in- and an out-neighbour */
      uint16_t a2[MAXN]; for (int v=0;v<9;v++) a2[v]=adj[v]; a2[9]=0;
      for (int v=0;v<9;v++){ if (Pt>>v&1) a2[9] |= 1u<<v; else a2[v] |= 1u<<9; }
      if (!strongq(a2, 10)) continue;
      u64 c = canon(a2, 10);
      if (hs.cnt > hs.cap/2){ memset(hs.key,0xff,hs.cap*sizeof(u64)); hs.cnt=0; }   /* bounded memory: reset (only repeats work) */
      if (!hs_insert(&hs, c)) continue;
      strong_new++;
      uint16_t t[MAXN]; decode(c, 10, t); u64 al[5]; packing(t, 10, al);
      u64 h3 = 0, p3 = 1; for (int k=0;k<5;k++){ h3 += al[k]*p3; p3 *= 3; }
      if (al[3]) alpha3pos++;
      if (h3 >= VMAX){ fprintf(stderr, "H_3 too large: %llu\n", h3); return 2; }
      if (wmode && target[h3] && !seen[h3]){ wit[h3] = c; got++; }
      seen[h3] = 1; if (h3 > maxh3) maxh3 = h3;
      if (wmode && got == need) goto done;
      if ((strong_new & 255) == 0){ u64 h2 = 0, p2 = 1; for (int k=0;k<5;k++){ h2 += al[k]*p2; p2 *= 2; } checks++; if (h2 != hpaths(t, 10)) bad++; }
    }
  }
  done:;
  char fn[512];
  if (wmode){ snprintf(fn, sizeof fn, "%s_%d.wit", argv[4], slice); FILE *fw = fopen(fn, "w");
    for (int v=0; v<VMAX; v++) if (target[v] && seen[v]) fprintf(fw, "%d %llu\n", v, wit[v]); fclose(fw);
    fprintf(stderr, "slice %d/%d witness mode: %ld of %ld target values witnessed after %ld classes\n", slice, ns, got, need, strong_new);
    return got == need ? 0 : 5; }
  snprintf(fn, sizeof fn, "%s_%d.bin", argv[4], slice); FILE *fo = fopen(fn, "wb"); fwrite(seen, 1, VMAX, fo); fclose(fo);
  fprintf(stderr, "slice %d/%d: %zu order-9 reps, strong classes processed %ld, OCF checks %ld (bad %ld), alpha_3 > 0: %ld, max H_3 %llu\n",
          slice, ns, (nr - slice + ns - 1)/ns, strong_new, checks, bad, alpha3pos, maxh3);
  return bad ? 3 : 0;
}
