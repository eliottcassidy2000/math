/* Recursive-tree tournaments: vertices 0..n-1 in order, transitive except child->parent arcs (parent < child).
   For every tree record H(T) and the leaf offsets delta_q = H(T + leaf at q) - 2 H(T) (new vertex beats only q).
   Output: for each odd value M, which of the offsets +1, -1, +3, -3 occur among trees with H = M. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
#define MAXN 14
static int n; static int par[MAXN];
static ll *dp; static uint32_t out[MAXN+1];
static ll H(int m){ int full=(1<<m)-1; memset(dp,0,sizeof(ll)*(size_t)(1<<m)*m);
  for(int v=0;v<m;v++) dp[(size_t)(1<<v)*m+v]=1;
  for(int S=1;S<=full;S++) for(int v=0;v<m;v++){ ll c=dp[(size_t)S*m+v]; if(!c) continue; uint32_t nb=out[v]&~S;
     while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|1u<<u)*m+u]+=c;} }
  ll h=0; for(int v=0;v<m;v++) h+=dp[(size_t)full*m+v]; return h; }
static void build(int m, const int *pp){ for(int i=0;i<m;i++) out[i]=0;
  for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(pp[j]==i) out[j]|=1u<<i; else out[i]|=1u<<j; } }
static int LIM; static unsigned char *flags; /* bit0:+1 bit1:-1 bit2:+3 bit3:-3 bit4: value occurs */
static ll trees=0;
static void rec(int i){
  if(i==n){ trees++; build(n,par); ll h=H(n); if(h>=LIM) return; flags[h]|=16;
    for(int q=0;q<n;q++){ int pp[MAXN+1]; memcpy(pp,par,sizeof(int)*n); pp[n]=q; build(n+1,pp); ll h2=H(n+1); ll d=h2-2*h;
      if(d==1) flags[h]|=1; else if(d==-1) flags[h]|=2; else if(d==3) flags[h]|=4; else if(d==-3) flags[h]|=8; }
    return; }
  for(int p=0;p<i;p++){ par[i]=p; rec(i+1);} }
int main(int argc,char**argv){ n=atoi(argv[1]); LIM=atoi(argv[2]);
  dp=malloc(sizeof(ll)*(size_t)(1<<(n+1))*(n+1)); flags=calloc(LIM,1); par[0]=-1; rec(1);
  int vals=0,both=0,plus=0,minus=0,none=0; int firstmiss=-1;
  for(int h=1;h<LIM;h+=2){ if(!(flags[h]&16)) continue; vals++; int p=flags[h]&1, m=(flags[h]>>1)&1;
    if(p&&m) both++; else if(p) plus++; else if(m) minus++; else { none++; if(firstmiss<0||h<firstmiss) {} }
  }
  printf("n=%d trees=%lld values<%d: %d | both +1,-1: %d | only +1: %d | only -1: %d | neither: %d\n",n,trees,LIM,vals,both,plus,minus,none);
  printf("values with neither unit offset (first 40):");
  int c=0; for(int h=1;h<LIM&&c<40;h+=2) if((flags[h]&16) && !(flags[h]&3)){ printf(" %d",h); c++; } printf("\n");
  printf("values with only +1 (first 30):"); c=0; for(int h=1;h<LIM&&c<30;h+=2) if((flags[h]&16) && (flags[h]&3)==1){ printf(" %d",h); c++; } printf("\n");
  printf("values with only -1 (first 30):"); c=0; for(int h=1;h<LIM&&c<30;h+=2) if((flags[h]&16) && (flags[h]&3)==2){ printf(" %d",h); c++; } printf("\n");
  return 0; }
