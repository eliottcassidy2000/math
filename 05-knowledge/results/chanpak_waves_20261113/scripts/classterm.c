// Strict block words: initial run (any O/Z word, length <= L0), then blocks J w with w in a TERMINAL set (J only after a
// complete block; "e" = empty block). Values counted at block ends. Misses per window classified by primality.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int L0, NB, EMPTYOK; static unsigned PB[512]; static int PLEN[512]; static unsigned char *comp;
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static int isprefix(unsigned code,int len){ for(int b=0;b<NB;b++) if(PLEN[b]>=len && (PB[b]&((1u<<len)-1))==code) return 1; return 0; }
static int isterm(unsigned code,int len){ if(len==0) return EMPTYOK; for(int b=0;b<NB;b++) if(PLEN[b]==len && PB[b]==code) return 1; return 0; }
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,unsigned code,int init){
  int term = init || isterm(code,run);
  if(term) SETB(H);
  i64 D=L-E; if(D>X-H) return;
  if(init ? run<L0 : isprefix(code,run+1)) { if(2*L<=X-H) dfs(H,0,2*L,h,run+1,code,init); }
  { unsigned c2=code|(1u<<run);
    if((init ? run<L0 : isprefix(c2,run+1)) && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,run+1,c2,init); }
  if(term) dfs(H+D,H-E,H+L,H,0,0,0);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); L0=atoi(argv[2]); NB=0; EMPTYOK=0;
  for(int i=3;i<argc;i++){ const char*w=argv[i]; if(!strcmp(w,"e")){ EMPTYOK=1; continue; } unsigned c=0; int l=(int)strlen(w); for(int t=0;t<l;t++) if(w[t]=='Z') c|=1u<<t; PB[NB]=c; PLEN[NB]=l; NB++; }
  bits=calloc((X>>7)+2,sizeof(u64)); comp=calloc(X/2+2,1);
  for(i64 p=3;p*p<=X;p+=2) if(!comp[p>>1]) for(i64 q=p*p;q<=X;q+=2*p) comp[q>>1]=1;
  dfs(1,1,1,1,0,0,1);
  printf("terminal blocks:"); for(int i=3;i<argc;i++) printf(" %s",argv[i]); printf("  (L0=%d)\n",L0);
  for(int k=12;(1LL<<k)<X;k++){ i64 lo=1LL<<k,hi=1LL<<(k+1); if(hi>X) hi=X; i64 m=0,mp=0,odd=0,pr=0,m3=0,m35=0;
    for(i64 n=lo|1;n<hi;n+=2){ odd++; int isp=!comp[n>>1]; if(isp) pr++;
      if(!GETB(n)){ m++; if(isp) mp++; if(n%3==0) m3++; if(n%3==0||n%5==0) m35++; } }
    printf("  k=%2d missed share %.3f | primes: %.1f%% of misses vs %.1f%% of odds | multiples of 3: %.1f%% of misses (33.3%% of odds) | of 3 or 5: %.1f%%\n",
      k,(double)m/odd,m?100.0*mp/m:0.0,100.0*pr/odd,m?100.0*m3/m:0.0,m?100.0*m35/m:0.0); }
  return 0; }
