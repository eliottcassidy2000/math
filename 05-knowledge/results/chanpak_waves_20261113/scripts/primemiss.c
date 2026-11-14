// List misses of the strict 11-block orbit up to X with factorizations: which are primes, 7p, or products of two
// realized numbers (factor a in [3, sqrt] with a and N/a both attained by the orbit or small known counts).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int L0, NB, EMPTYOK; static unsigned PB[64]; static int PLEN[64];
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static int isprefix(unsigned code,int len){ for(int b=0;b<NB;b++) if(PLEN[b]>=len && (PB[b]&((1u<<len)-1))==code) return 1; return 0; }
static int isterm(unsigned code,int len){ if(len==0) return EMPTYOK; for(int b=0;b<NB;b++) if(PLEN[b]==len && PB[b]==code) return 1; return 0; }
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,unsigned code,int init){
  int term = init || isterm(code,run); if(term) SETB(H);
  i64 D=L-E; if(D>X-H) return;
  if(init ? run<L0 : isprefix(code,run+1)) { if(2*L<=X-H) dfs(H,0,2*L,h,run+1,code,init); }
  { unsigned c2=code|(1u<<run); if((init ? run<L0 : isprefix(c2,run+1)) && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,run+1,c2,init); }
  if(term) dfs(H+D,H-E,H+L,H,0,0,0);
}
static int realized(i64 n){ return !(n==7||n==21); } /* every odd n < 2^32 except 7, 21 is a count (THM-4875 (D)) */
int main(int argc,char**argv){
  X=atoll(argv[1]); L0=atoi(argv[2]); NB=0; EMPTYOK=0;
  for(int i=3;i<argc;i++){ const char*w=argv[i]; if(!strcmp(w,"e")){ EMPTYOK=1; continue; } unsigned c=0; int l=(int)strlen(w); for(int t=0;t<l;t++) if(w[t]=='Z') c|=1u<<t; PB[NB]=c; PLEN[NB]=l; NB++; }
  bits=calloc((X>>7)+2,sizeof(u64)); dfs(1,1,1,1,0,0,1);
  i64 tot=0, np=0, nprod=0, other=0, lastp=0, lastother=0;
  for(i64 n=3;n<=X;n+=2) if(!GETB(n)){ tot++;
      int isp=1; i64 fac=0; for(i64 d=3; d*d<=n; d+=2) if(n%d==0){ isp=0; if(!fac && realized(d) && realized(n/d)) fac=d; }
      if(isp){ np++; lastp=n; } else if(fac) nprod++; else { other++; lastother=n; } }
  printf("strict 11-block orbit up to %lld: %lld misses; %lld prime (largest %lld); %lld are products of two counts; %lld other (largest %lld)\n",X,tot,np,lastp,nprod,other,lastother);
  return 0; }
