// Representation counts r(N) for strict-block root-change words: initial run (any O/Z word, length <= L0), then blocks
// J w with w in a terminal set (J only after a complete block; "e" = empty block). A representation is a word ending at
// a block end (or within the initial run). Reports r over odd N in [X/2, X]: min/mean overall, by class mod 3, 5, 7,
// and for primes vs composites; also the N attaining the minimum.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static unsigned *cnt; static int L0, NB, EMPTYOK; static unsigned PB[64]; static int PLEN[64];
static int isprefix(unsigned code,int len){ for(int b=0;b<NB;b++) if(PLEN[b]>=len && (PB[b]&((1u<<len)-1))==code) return 1; return 0; }
static int isterm(unsigned code,int len){ if(len==0) return EMPTYOK; for(int b=0;b<NB;b++) if(PLEN[b]==len && PB[b]==code) return 1; return 0; }
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,unsigned code,int init){
  int term = init || isterm(code,run);
  if(term && cnt[H>>1]<4000000000u) cnt[H>>1]++;
  i64 D=L-E; if(D>X-H) return;
  if(init ? run<L0 : isprefix(code,run+1)) { if(2*L<=X-H) dfs(H,0,2*L,h,run+1,code,init); }
  { unsigned c2=code|(1u<<run);
    if((init ? run<L0 : isprefix(c2,run+1)) && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,run+1,c2,init); }
  if(term) dfs(H+D,H-E,H+L,H,0,0,0);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); L0=atoi(argv[2]); NB=0; EMPTYOK=0;
  for(int i=3;i<argc;i++){ const char*w=argv[i]; if(!strcmp(w,"e")){ EMPTYOK=1; continue; } unsigned c=0; int l=(int)strlen(w); for(int t=0;t<l;t++) if(w[t]=='Z') c|=1u<<t; PB[NB]=c; PLEN[NB]=l; NB++; }
  cnt=calloc(X/2+2,sizeof(unsigned)); unsigned char*comp=calloc(X/2+2,1);
  for(i64 p=3;p*p<=X;p+=2) if(!comp[p>>1]) for(i64 q=p*p;q<=X;q+=2*p) comp[q>>1]=1;
  dfs(1,1,1,1,0,0,1);
  double s=0; unsigned mn=~0u; i64 argmn=0, n0=0;
  double sc[3][7]={{0}}; double nc[3][7]={{0}}; unsigned mc[3][7]; for(int a=0;a<3;a++) for(int b=0;b<7;b++) mc[a][b]=~0u;
  int mods[3]={3,5,7}; double sp=0,np=0,scm=0,ncm=0; unsigned mp=~0u, mcomp=~0u; i64 argmp=0;
  for(i64 n=(X/2)|1;n<=X;n+=2){ unsigned r=cnt[n>>1]; s+=r; n0++; if(r<mn){mn=r;argmn=n;}
    for(int a=0;a<3;a++){ int c=n%mods[a]; sc[a][c]+=r; nc[a][c]++; if(r<mc[a][c]) mc[a][c]=r; }
    if(!comp[n>>1]){ sp+=r; np++; if(r<mp){mp=r;argmp=n;} } else { scm+=r; ncm++; if(r<mcomp) mcomp=r; } }
  printf("X=%lld window [X/2,X]: mean r %.1f, min r %u at N=%lld\n",X,s/n0,mn,argmn);
  printf("  primes: mean %.1f min %u (at %lld) | composites: mean %.1f min %u\n",sp/np,mp,argmp,scm/ncm,mcomp);
  for(int a=0;a<3;a++){ printf("  mod %d:",mods[a]); for(int c=0;c<mods[a];c++) printf("  [%d] mean %.1f min %u",c,sc[a][c]/nc[a][c],mc[a][c]); printf("\n"); }
  return 0; }
