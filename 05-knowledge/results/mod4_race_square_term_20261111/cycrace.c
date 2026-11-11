// cycrace.c - the parity race by cycle length.  For each tournament class of order n: the numbers c_k of directed
// cycles of each odd length k = 3, 5, ..., and the labelled weight n!/|Aut|.  Accumulates, for every subset K of the odd
// lengths, the signed labelled sum of (-1)^(sum_{k in K} c_k).   Usage: cycrace n  < gentourng output
// Output: one line per subset K (bitmask over lengths 3,5,7,9,11): mask, signed sum, total weight.
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 11
static int n; static u64 p[1<<MAXN][MAXN]; static int outm[MAXN], sc[MAXN];
static int perm[MAXN], used; static u64 autcount;
static void autrec(int k){
  if(k==n){ autcount++; return; }
  for(int c=0;c<n;c++){ if(used>>c&1) continue; if(sc[c]!=sc[k]) continue; int ok=1;
    for(int j=0;j<k;j++){ int a=(outm[k]>>j)&1, b=(outm[c]>>perm[j])&1; if(a!=b){ok=0;break;} }
    if(!ok) continue; perm[k]=c; used|=1<<c; autrec(k+1); used&=~(1<<c); }
}
int main(int argc,char**argv){ n=atoi(argv[1]); int full=(1<<n)-1; char line[256];
  u64 fact=1; for(int i=2;i<=n;i++) fact*=i;
  int nl=0; for(int k=3;k<=n;k+=2) nl++;
  static i64 sgn[64]; u64 wtot=0, nt=0; static u64 cnt[MAXN+1];
  u64 chk=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++) outm[i]=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1') outm[i]|=1<<j; else outm[j]|=1<<i; }
    nt++;
    // cycles counted by their least vertex v: paths from v inside {v..n-1} ending at w with w->v
    for(int L=0;L<=n;L++) cnt[L]=0;
    for(int v=0;v<n;v++){
      int allowed=full&~((1<<v)-1); // vertices >= v
      // p[S][w] = number of paths from v with vertex set S (v in S), ending at w
      for(int S=1<<v;S<=full;S++){ if((S&~allowed)||!(S>>v&1)) continue;
        for(int w=v;w<n;w++){ if(!(S>>w&1)){p[S][w]=0;continue;}
          if(S==(1<<v)){ p[S][w]=(w==v); continue; }
          if(w==v){p[S][w]=0;continue;}
          int R=S&~(1<<w); u64 a=0; for(int u=v;u<n;u++) if((R>>u&1)&&(outm[u]>>w&1)) a+=p[R][u]; p[S][w]=a; }
        int sz=__builtin_popcount(S); if(sz>=3){ u64 c=0; for(int w=v+1;w<n;w++) if((S>>w&1)&&(outm[w]>>v&1)) c+=p[S][w]; cnt[sz]+=c; }
      }
    }
    for(int i=0;i<n;i++) sc[i]=__builtin_popcount(outm[i]); autcount=0; used=0; autrec(0);
    u64 w=fact/autcount; wtot+=w;
    // check c3 against the score formula
    { i64 c3=(i64)n*(n-1)*(n-2)/6; for(int i=0;i<n;i++) c3-=(i64)sc[i]*(sc[i]-1)/2; if((u64)c3!=cnt[3]) chk++; }
    int par=0; { int b=0; for(int L=3;L<=n;L+=2){ if(cnt[L]&1) par|=1<<b; b++; } }
    for(int K=0;K<(1<<nl);K++){ int s=__builtin_popcount(par&K)&1; sgn[K]+= s?-(i64)w:(i64)w; }
  }
  for(int K=0;K<(1<<nl);K++) printf("%d %lld %llu\n", K, sgn[K], wtot);
  fprintf(stderr,"order %d: classes %llu, c3 check failures %llu\n", n, nt, chk);
  return 0; }
