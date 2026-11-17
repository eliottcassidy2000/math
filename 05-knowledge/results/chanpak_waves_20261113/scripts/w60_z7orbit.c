/* Single moved Z7-orbit: Q[B at X] with B = I7 (interval circulant {1,2,3}) or P7 (Paley/Fano {1,2,4}).
   H(Q[B]) = sum_r a_r(Q,X) c_r(B): a_r = walks in Q visiting X exactly r times and every other vertex once;
   c_r(B) = ordered r-covers of B by nonempty paths.  Z7 acts freely, so 7 | c_r and the orbit count is H/7. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 Hsub[128];
static void allH(int n,const int*o){ static u64 dp[128][7]; int N=1<<n; for(int S=0;S<N;S++) for(int v=0;v<n;v++) dp[S][v]=0; for(int v=0;v<n;v++) dp[1<<v][v]=1;
  for(int S=1;S<N;S++) for(int v=0;v<n;v++) if(dp[S][v]) for(int w=0;w<n;w++) if(!(S>>w&1)&&(o[v]>>w&1)) dp[S|1<<w][w]+=dp[S][v];
  Hsub[0]=1; for(int S=1;S<N;S++){ u64 h=0; for(int v=0;v<n;v++) h+=dp[S][v]; Hsub[S]=h; } }
static void covers(const int*o,u64*c){ /* c[r] for r=1..7: sum over surjective maps V->[r] of prod H(fibres) */
  allH(7,o); for(int r=0;r<=7;r++) c[r]=0; int f[7];
  for(int code=0;code<823543;code++){ int x=code, used=0, mask[7]={0}; for(int i=0;i<7;i++){ f[i]=x%7; x/=7; mask[f[i]]|=1<<i; used|=1<<f[i]; }
    int r=__builtin_popcount(used); if(used!=(1<<r)-1) continue; /* labels must be exactly 0..r-1 */
    u64 pr=1; for(int k=0;k<r;k++) pr*=Hsub[mask[k]]; c[r]+=pr; } }
static int isprime(u64 q){ if(q<2) return 0; for(u64 d=2;d*d<=q;d++) if(q%d==0) return 0; return 1; }
static int strongQ(int n,const int*o){ int full=(1<<n)-1,r=1,fr=1; while(fr){ int nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=o[v]; nf&=~r; r|=nf; fr=nf; } if(r!=full) return 0;
  r=1; fr=1; while(fr){ int nf=0; for(int v=0;v<n;v++) if(fr>>v&1) for(int u=0;u<n;u++) if(o[u]>>v&1) nf|=1<<u; nf&=~r; r|=nf; fr=nf; } return r==full; }
int main(int argc,char**argv){ u64 LIM=strtoull(argv[1],0,10);
  int I7[7]={0},P7[7]={0}; for(int i=0;i<7;i++) for(int d=1;d<7;d++){ int j=(i+d)%7; if(d==1||d==2||d==3) I7[i]|=1<<j; if(d==1||d==2||d==4) P7[i]|=1<<j; }
  u64 cI[8],cP[8]; covers(I7,cI); covers(P7,cP);
  printf("c_r(I7):"); for(int r=1;r<=7;r++) printf(" %llu",cI[r]); printf("\nc_r(P7):"); for(int r=1;r<=7;r++) printf(" %llu",cP[r]); printf("\n");
  printf("c_r/7 (I7):"); for(int r=1;r<=7;r++) printf(" %llu",cI[r]/7); printf("  all divisible by 7: %d\n", (int)(cI[1]%7==0&&cI[2]%7==0&&cI[3]%7==0));
  unsigned char *pr=calloc(LIM+1,1), *pr2=calloc(LIM+1,1); FILE*fi=fopen("w59_reps8.txt","r"); char line[256]; long cfg=0, strongcfg=0; u64 smallest=~0ULL; int sm_n=0,sm_o[8],sm_x=0,sm_b=0;
  static u64 dp[8][256][9];
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); int o[8]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10); if(n<2) continue; int st=strongQ(n,o);
    for(int X=0;X<n;X++){ /* walks: visit X r times (r<=7), others once */
      int others=((1<<n)-1)&~(1<<X); memset(dp,0,sizeof dp); u64 a[9]={0};
      for(int v=0;v<n;v++){ if(v==X) dp[v][0][1]=1; else dp[v][1<<v][0]=1; }
      for(int len=1; len<=n-1+7; len++) for(int v=0;v<n;v++) for(int m=0;m<(1<<n);m++){ if(m>>X&1) continue; for(int r=0;r<=7;r++){ u64 c=dp[v][m][r]; if(!c) continue; if(__builtin_popcount(m)+r!=len) continue;
            if(m==others && r>=1) a[r]+=c;
            for(int w=0;w<n;w++){ if(w==v||!(o[v]>>w&1)) continue; if(w==X){ if(r<7) dp[w][m][r+1]+=c; } else if(!(m>>w&1)) dp[w][m|1<<w][r]+=c; } } }
      int third=0; for(int r=3;r<=7;r++) if(a[r]) third=1;
      for(int b=0;b<2;b++){ u64*c=b?cP:cI; u64 H=0; for(int r=1;r<=7;r++) H+=a[r]*c[r]; u64 orb=H/7; cfg++; if(st) strongcfg++;
        if(st && orb<=LIM && isprime(orb)){ pr[orb]=1; if(!third) pr2[orb]=1; if(orb<smallest){ smallest=orb; sm_n=n; memcpy(sm_o,o,sizeof o); sm_x=X; sm_b=b; } } } } }
  long np=0,np2=0; for(u64 q=0;q<=LIM;q++){ if(pr[q]) np++; if(pr2[q]) np2++; }
  long primes=0; for(u64 q=2;q<=LIM;q++) if(isprime(q)) primes++;
  printf("configurations %ld (strong quotient %ld); distinct prime orbit counts <= %llu: %ld of %ld primes; with no third run at X: %ld\n",cfg,strongcfg,LIM,np,primes,np2);
  printf("smallest prime orbit count %llu: |Q|=%d X=%d block %s\n",smallest,sm_n,sm_x,sm_b?"P7":"I7");
  printf("first prime orbit counts:"); int k=0; for(u64 q=0;q<=LIM&&k<30;q++) if(pr[q]){ printf(" %llu",q); k++; } printf("\n");
  return 0; }
