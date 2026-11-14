// Root-change words: initial run (any O/Z word, length <= L0), then blocks J w with w a prefix of an allowed word
// (prefix-closed allowed set; pass "ALL<c>" to allow every word of length <= c, "NOZJ<c>" for words of length <= c
// not starting with Z). Misses per dyadic window classified as prime / 7*prime / composite-other.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int L0, NB, MODE, CAP; static unsigned PB[512]; static int PLEN[512]; static unsigned char *comp;
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static int ok(unsigned code,int len){
  if(MODE==1) return len<=CAP;
  if(MODE==2) return len<=CAP && !(len>=1 && (code&1));
  for(int b=0;b<NB;b++) if(PLEN[b]>=len && (PB[b]&((1u<<len)-1))==code) return 1; return 0; }
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,unsigned code,int init){
  SETB(H); i64 D=L-E; if(D>X-H) return;
  if(init ? run<L0 : ok(code,run+1)) { if(2*L<=X-H) dfs(H,0,2*L,h,run+1,code,init); }
  { unsigned c2=code|(1u<<run);
    if((init ? run<L0 : ok(c2,run+1)) && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,run+1,c2,init); }
  dfs(H+D,H-E,H+L,H,0,0,0);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); L0=atoi(argv[2]); NB=0; MODE=0;
  for(int i=3;i<argc;i++){ const char*w=argv[i];
    if(!strncmp(w,"ALL",3)){ MODE=1; CAP=atoi(w+3); continue; } if(!strncmp(w,"NOZJ",4)){ MODE=2; CAP=atoi(w+4); continue; }
    unsigned c=0; int l=(int)strlen(w); for(int t=0;t<l;t++) if(w[t]=='Z') c|=1u<<t; PB[NB]=c; PLEN[NB]=l; NB++; }
  bits=calloc((X>>7)+2,sizeof(u64)); comp=calloc(X/2+2,1);
  for(i64 p=3;p*p<=X;p+=2) if(!comp[p>>1]) for(i64 q=p*p;q<=X;q+=2*p) comp[q>>1]=1;
  dfs(1,1,1,1,0,0,1);
  printf("family:"); for(int i=3;i<argc;i++) printf(" %s",argv[i]); printf("  (L0=%d)\n",L0);
  for(int k=12;(1LL<<k)<X;k++){ i64 lo=1LL<<k,hi=1LL<<(k+1); if(hi>X) hi=X; i64 m=0,mp=0,m7=0,lastp=0,odd=0,pr=0;
    for(i64 n=lo|1;n<hi;n+=2){ odd++; int isp=!comp[n>>1]; if(isp) pr++;
      if(!GETB(n)){ m++; if(isp){mp++;lastp=n;} else if(n%7==0 && !comp[(n/7)>>1]) m7++; } }
    printf("  k=%2d misses %8lld | prime %7lld (%.1f%% of misses; primes are %.1f%% of odds) | 7p %5lld | last prime miss %lld\n",
      k,m,mp,m?100.0*mp/m:0.0,100.0*pr/odd,m7,lastp); }
  return 0; }
