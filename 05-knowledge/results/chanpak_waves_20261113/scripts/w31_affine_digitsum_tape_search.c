/* Search for "linear tapes": one-apex hosts [g1][tape of gadgets B_bit][core][g2] + apex whose
   Hamiltonian-path count is exactly affine in X (tape = binary digits of X, MSB first, fixed length L).
   Gadgets: k=1,2,3 spine vertices with internal back arcs (span<=2) and apex marks (mark=1: apex beats it). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
#define MAXN 24
static int n; static uint32_t out[MAXN];
static ll *dp;
static ll Hcount(void){
  int full=(1<<n)-1; size_t sz=(size_t)(1<<n)*n; memset(dp,0,sz*sizeof(ll));
  for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
  for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ ll c=dp[(size_t)S*n+v]; if(!c) continue;
      uint32_t nb=out[v]&~S; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|1<<u)*n+u]+=c;} }
  ll H=0; for(int v=0;v<n;v++) H+=dp[(size_t)full*n+v]; return H; }
typedef struct { int k; int back[3][3]; int mark[3]; } gad;
static gad G[200]; static int ng=0;
static void addg(int k,int b01,int b02,int b12,int mk){ gad g; memset(&g,0,sizeof g); g.k=k;
  g.back[0][1]=b01; g.back[0][2]=b02; g.back[1][2]=b12; for(int i=0;i<k;i++) g.mark[i]=(mk>>i)&1; G[ng++]=g; }
/* spine: list of (vertex, mark) with back arcs given as pairs (later->earlier) */
static int sm[MAXN]; static int bk[MAXN][MAXN]; static int ns;
static void reset(){ ns=0; memset(bk,0,sizeof bk); }
static int addv(int mark){ sm[ns]=mark; return ns++; }
static void addgad(const gad*g){ int base=ns; for(int i=0;i<g->k;i++) addv(g->mark[i]);
  for(int i=0;i<g->k;i++) for(int j=i+1;j<g->k;j++) if(g->back[i][j]) bk[base+j][base+i]=1; }
static void build(){ n=ns+1; for(int i=0;i<n;i++) out[i]=0;
  for(int i=0;i<ns;i++) for(int j=i+1;j<ns;j++){ if(bk[j][i]) out[j]|=1u<<i; else out[i]|=1u<<j; }
  int x=ns; for(int i=0;i<ns;i++){ if(sm[i]) out[x]|=1u<<i; else out[i]|=1u<<x; } }
int main(int argc,char**argv){
  int L=atoi(argv[1]); int kmax=atoi(argv[2]);
  dp=malloc(sizeof(ll)*(size_t)(1<<MAXN)/8*MAXN/1); /* sized for n<=21 */
  for(int mk=0;mk<2;mk++) addg(1,0,0,0,mk);
  if(kmax>=2) for(int b=0;b<2;b++) for(int mk=0;mk<4;mk++) addg(2,b,0,0,mk);
  if(kmax>=3) for(int b01=0;b01<2;b01++) for(int b02=0;b02<2;b02++) for(int b12=0;b12<2;b12++) for(int mk=0;mk<8;mk++) addg(3,b01,b02,b12,mk);
  int hits=0; long tested=0;
  for(int i0=0;i0<ng;i0++) for(int i1=0;i1<ng;i1++){ if(i0==i1) continue; if(G[i0].k!=G[i1].k) continue;
    for(int cm=0;cm<8;cm++) for(int g1m=0;g1m<2;g1m++) for(int g2m=0;g2m<2;g2m++){
      ll Hs[64]; int ok=1;
      for(int X=0;X<(1<<L);X++){
        reset(); addv(g1m);
        for(int b=L-1;b>=0;b--) addgad(&G[((X>>b)&1)?i1:i0]);
        int c0=addv((cm>>0)&1), c1=addv((cm>>1)&1), c2=addv((cm>>2)&1); bk[c2][c0]=1; /* 3-cycle core */
        (void)c1; addv(g2m);
        build(); Hs[X]=Hcount(); tested++;
        if(X>=3){ ll B=Hs[2]-Hs[0]; B/=2; ll C=Hs[1]-Hs[0]-B; if(2*B!=Hs[2]-Hs[0]){ok=0;break;} ll pred=Hs[0]+B*X+C*__builtin_popcount(X); if(Hs[X]!=pred){ ok=0; break; } }
      }
      if(ok && Hs[1]!=Hs[0]){ hits++; if(hits<=20) printf("AFFINE+DIGITSUM: k=%d gadgets %d,%d core %d guards %d%d : H = %lld + %lld X + %lld s2(X)\n",G[i0].k,i0,i1,cm,g1m,g2m,Hs[0],(Hs[2]-Hs[0])/2,Hs[1]-Hs[0]-(Hs[2]-Hs[0])/2); }
    } }
  printf("L=%d kmax=%d gadgets=%d tested=%ld affine hits=%d\n",L,kmax,ng,tested,hits);
  return 0; }
