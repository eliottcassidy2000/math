// antiinv_census.c - all labelled tournaments on n = 2m + 1 vertices that are anti-invariant under the involution psi
// fixing vertex 0 and swapping 2i-1 <-> 2i (i = 1..m), normalised so that 0 -> 2i-1 for every i (swapping the names
// inside a pair is a relabelling).  Free arcs: d_i = [2i-1 -> 2i] (m bits), and for i < j the two bits
// [a_i -> a_j], [a_i -> psi(a_j)] (a_i = 2i-1); the other arcs are forced by psi.  m^2 bits in all.
// For each tournament: t = number of cyclic invariant triples {0, a_i, psi(a_i)} (= number of d_i = 1), the numbers
// F_k of invariant odd cycles of length 2k+1 (half-path DP), F = sum F_k, and for n <= 9 the residue of H mod 4.
// Output: for each t: number of tournaments, number with F odd, histogram of the parity vector (F_1..F_m), and the
// number of violations of H = 1 + 2F (mod 4).
// Check of the quotient formula (Theorem Q of the note): let Z be the tournament on the pairs with
//   i -> j  iff  [a_i -> a_j] + [a_i -> psi(a_j)] is odd,
// and D_1 the set of pairs with a_i -> psi(a_i).  Then F_k = number of directed paths on k vertices of Z that end in
// D_1 (mod 2), for every k.  The program counts the violations.
// Usage: antiinv_census m
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(int argc,char**argv){ int m=atoi(argv[1]); int n=2*m+1; int nb=m*m;
  static u64 cnt[8][2]; static u64 pat[8][64]; u64 viol=0, total=0, violQ=0; int doH=(n<=9);
  static unsigned dph[1<<11][11];
  for(u64 bits=0; bits<(1ULL<<nb); bits++){
    unsigned out[16]; memset(out,0,sizeof out); int b=0;
    for(int i=1;i<=m;i++){ out[0]|=1u<<(2*i-1); out[2*i]|=1u<<0; }
    int t=0;
    for(int i=1;i<=m;i++){ int d=(bits>>b)&1; b++; t+=d; if(d) out[2*i-1]|=1u<<(2*i); else out[2*i]|=1u<<(2*i-1); }
    for(int i=1;i<=m;i++) for(int j=i+1;j<=m;j++){
      int x=(bits>>b)&1; b++; int ai=2*i-1, aj=2*j-1, pi=2*i, pj=2*j;
      if(x){ out[ai]|=1u<<aj; out[pj]|=1u<<pi; } else { out[aj]|=1u<<ai; out[pi]|=1u<<pj; }
      int y=(bits>>b)&1; b++;
      if(y){ out[ai]|=1u<<pj; out[aj]|=1u<<pi; } else { out[pj]|=1u<<ai; out[pi]|=1u<<aj; } }
    // half-path DP over pairs: vertex index v in 1..2m, pair (v-1)/2
    static unsigned dp[32][11]; memset(dp,0,sizeof dp);
    for(int v=1;v<=2*m;v++) if(out[0]>>v&1) dp[1<<((v-1)/2)][v]=1;
    unsigned Fk[8]; memset(Fk,0,sizeof Fk);
    for(int mask=1;mask<(1<<m);mask++){ int k=__builtin_popcount(mask);
      for(int v=1;v<=2*m;v++){ unsigned c=dp[mask][v]; if(!c) continue;
        int pv=((v-1)^1)+1; if(out[v]>>pv&1) Fk[k]+=c;
        for(int w=1;w<=2*m;w++){ int pw=(w-1)/2; if(mask>>pw&1) continue; if(out[v]>>w&1) dp[mask|(1<<pw)][w]+=c; } } }
    unsigned F=0; int pv=0; for(int k=1;k<=m;k++){ F+=Fk[k]; if(Fk[k]&1) pv|=1<<(k-1); }
    { // quotient tournament Z and D_1
      unsigned zout[8]; memset(zout,0,sizeof zout); unsigned D1=0;
      for(int i=1;i<=m;i++){ if(out[2*i-1]>>(2*i)&1) D1|=1u<<(i-1);
        for(int j=1;j<=m;j++){ if(i==j) continue; int x=(out[2*i-1]>>(2*j-1))&1, y=(out[2*i-1]>>(2*j))&1; if(x^y) zout[i-1]|=1u<<(j-1); } }
      for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ int a=(zout[i]>>j)&1, bb=(zout[j]>>i)&1; if(a==bb){ fprintf(stderr,"Z is not a tournament\n"); return 1; } }
      static unsigned zp[32][8]; memset(zp,0,sizeof zp); unsigned Qk[8]; memset(Qk,0,sizeof Qk);
      for(int v=0;v<m;v++) zp[1<<v][v]=1;
      for(int mask=1;mask<(1<<m);mask++){ int k=__builtin_popcount(mask);
        for(int v=0;v<m;v++){ unsigned c=zp[mask][v]; if(!c) continue; if(D1>>v&1) Qk[k]+=c;
          unsigned nbm=zout[v]&~mask; while(nbm){ int w=__builtin_ctz(nbm); nbm&=nbm-1; zp[mask|(1<<w)][w]+=c; } } }
      for(int k=1;k<=m;k++) if((Qk[k]^Fk[k])&1){ violQ++; break; } }
    cnt[t][F&1]++; pat[t][pv]++; total++;
    if(doH){ int full=(1<<n)-1; for(int S=1;S<=full;S++) for(int v=0;v<n;v++) dph[S][v]=0;
      for(int v=0;v<n;v++) dph[1<<v][v]=1;
      for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ unsigned c=dph[S][v]; if(!c) continue; unsigned nbm=out[v]&~S;
        while(nbm){ int w=__builtin_ctz(nbm); nbm&=nbm-1; dph[S|(1<<w)][w]+=c; } }
      unsigned H=0; for(int v=0;v<n;v++) H+=dph[full][v];
      if(((H-1-2*F)&3)!=0) viol++; }
  }
  printf("n = %d (m = %d): %llu normalised anti-invariant labelled tournaments; violations of H = 1 + 2F (mod 4): %llu%s; violations of the quotient formula: %llu\n", n, m, total, viol, doH?"":" (H not computed)", violQ);
  for(int t=0;t<=m;t++){ printf("  t = %d cyclic invariant triples: %llu tournaments, F odd in %llu, F even in %llu; parity vectors (F_1..F_m, F_1 first):", t, cnt[t][0]+cnt[t][1], cnt[t][1], cnt[t][0]);
    for(int p=0;p<(1<<m);p++) if(pat[t][p]){ printf(" "); for(int k=0;k<m;k++) printf("%d",(p>>k)&1); printf(":%llu", pat[t][p]); }
    printf("\n"); }
  return 0; }
