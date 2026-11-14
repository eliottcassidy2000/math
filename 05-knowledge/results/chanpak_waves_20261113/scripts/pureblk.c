// Wave 24 exploratory: free one-apex source, then blocks from a restricted set B (on (H,E,lam)):
//   kind 0: J and J O^k Z only (one Z per block: zero pair energy), k>=1;
//   kind 1: J and all J O v Z with exactly two Z's.
// Block action computed letter by letter. Usage: pureblk K kind
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static uint8_t *cov; static int kind;
static void stepO(i64*s){ s[1]=0; s[2]*=2; }
static void stepZ(i64*s){ i64 H=s[0],E=s[1],L=s[2],h=s[3]; s[0]=H+L-E; s[1]=h; s[2]=2*L-h; }
static void stepJ(i64*s){ i64 H=s[0],E=s[1],L=s[2]; s[0]=H+L-E; s[1]=H-E; s[2]=H+L; s[3]=H; }
static void blocks(i64*s);
static void visit(i64*s){ if(s[0]<=X) cov[s[0]>>1]=1; else return; blocks(s); }
static void runs(i64*t,int zleft,int len){ // after J O..., place remaining letters; end with Z
  if(t[0]>X) return;
  i64 u[4]; 
  // option: end here with Z (if exactly zleft==1)
  if(zleft==1){ for(int i=0;i<4;i++)u[i]=t[i]; stepZ(u); if(u[0]<=X) visit(u); }
  // option: O
  for(int i=0;i<4;i++)u[i]=t[i]; stepO(u); if(u[0]+u[2]<=X+u[2] && 2*u[2] <= 2*(X+1)) { if(u[2]<=X) runs(u,zleft,len+1); }
  // option: Z in the middle (kind 1 only, needs zleft==2)
  if(zleft==2){ for(int i=0;i<4;i++)u[i]=t[i]; stepZ(u); if(u[0]<=X) runs(u,1,len+1); }
}
static void blocks(i64*s){
  i64 t[4]; for(int i=0;i<4;i++)t[i]=s[i]; stepJ(t); if(t[0]>X) return; visit(t);   // J alone
  for(int i=0;i<4;i++)t[i]=s[i]; stepJ(t); stepO(t); runs(t,kind==0?1:2,2);
}
static void source(i64*s,int depth){ visit(s); i64 u[4]; if(depth>60) return;
  for(int i=0;i<4;i++)u[i]=s[i]; stepO(u); if(u[2]<=X) source(u,depth+1);
  for(int i=0;i<4;i++)u[i]=s[i]; i64 D=s[2]-s[1]; if(!(D==0&&s[2]==s[3]&&s[1]==s[3])){ stepZ(u); if(u[0]<=X) source(u,depth+1);} }
int main(int argc,char**argv){ int K=atoi(argv[1]); kind=atoi(argv[2]); X=(1LL<<K)-1; cov=calloc(X/2+2,1);
  i64 s[4]={1,1,1,1}; source(s,0);
  u64 m=0; i64 last=0; for(i64 n=1;n<=X;n+=2) if(!cov[n>>1]){m++;last=n;}
  printf("kind=%d K=%d misses=%llu largest=%lld\n",kind,K,m,last);
  for(int k=12;k<K;k++){ u64 mm=0,nn=0; for(i64 n=(1LL<<k)+1;n<(1LL<<(k+1));n+=2){nn++; if(!cov[n>>1]) mm++;} printf("  2^%d: %.4f missed\n",k,(double)mm/nn); }
  return 0; }
