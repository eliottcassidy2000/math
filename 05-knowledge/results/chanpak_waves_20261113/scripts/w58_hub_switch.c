/* Wave 58: the hub-arc switch of two-hub tournaments.  Flipping the hub arc changes H by 2(K(x,y) - K(y,x)) (every odd
   cycle through the hub arc has a transitive complement).  pairs: unit pairs ((2^p+1)^2, +2); density: how many odd N are
   lower members of a unit pair; residues: lower values mod 9. */
/* Two-hub tournaments (spine TT_m, hubs u1 (word X), u2 (word Y), arc e = 1 iff u1->u2).
   Hub-arc flip gap: H(X,Y,1) - H(X,Y,0) = 2 * (N(u1->P->u2... ) difference), see notes.
   Find, for targets T_p = (2^p+1)^2, configurations with {H, H'} = {T_p, T_p+2} or {T_p-2, T_p} (flip of e). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n;
static u64 count_order(const int*A,const int*B,int ab){ u64 dp[27]={0}, nd[27]; dp[0]=1;
  for(int i=0;i<n;i++){ memset(nd,0,sizeof nd);
    for(int st=0;st<27;st++){ u64 c=dp[st]; if(!c) continue; int s1=st%3, s2=st/3%3, s3=st/9;
      { int t1=(A[i]==0)?1:2; nd[t1+3*s2+9*s3]+=c; }
      if(s2!=0 || A[i]==1){ int t2=(B[i]==0)?1:2; nd[s1+3*t2+9*s3]+=c; }
      if(s3!=0 || B[i]==1){ nd[s1+3*s2+9*1]+=c; } }
    memcpy(dp,nd,sizeof dp); }
  u64 tot=0; for(int st=0;st<27;st++){ int s1=st%3,s2=st/3%3,s3=st/9; if(s1==2||s2==2||s3==2) continue; if(s2==0&&!ab) continue; tot+=dp[st]; } return tot; }
static int main_pairs(int argc,char**argv){ int L=atoi(argv[1]); int PMAX=atoi(argv[2]); u64 T[24]; for(int p=1;p<=PMAX;p++){ u64 q=(1ULL<<p)+1; T[p]=q*q; }
  int fA[24]={0}, fB[24]={0}; int X[40],Y[40];
  for(n=1;n<=L;n++){ u64 W=1ULL<<n;
    for(u64 wx=0;wx<W;wx++) for(u64 wy=0;wy<W;wy++){
      for(int i=0;i<n;i++){ X[i]=wx>>i&1; Y[i]=wy>>i&1; }
      u64 h1=count_order(X,Y,1)+count_order(Y,X,0); /* e=1: u1->u2 */
      u64 h0=count_order(X,Y,0)+count_order(Y,X,1);
      u64 lo=h0<h1?h0:h1, hi=h0<h1?h1:h0; if(hi-lo!=2) continue;
      for(int p=1;p<=PMAX;p++){
        if(lo==T[p] && fA[p]<2){ fA[p]++; printf("p=%d pair (T, T+2) = (%llu, %llu): m=%d x=",p,lo,hi,n); for(int i=0;i<n;i++) printf("%d",X[i]); printf(" y="); for(int i=0;i<n;i++) printf("%d",Y[i]); printf(" (lower has e=%d)\n", h0<h1?0:1); fflush(stdout); }
        if(hi==T[p] && fB[p]<2){ fB[p]++; printf("p=%d pair (T-2, T) = (%llu, %llu): m=%d x=",p,lo,hi,n); for(int i=0;i<n;i++) printf("%d",X[i]); printf(" y="); for(int i=0;i<n;i++) printf("%d",Y[i]); printf(" (lower has e=%d)\n", h0<h1?0:1); fflush(stdout); } } } }
  for(int p=1;p<=PMAX;p++) printf("p=%d T=%llu: (T,T+2) %d, (T-2,T) %d\n",p,T[p],fA[p],fB[p]); return 0; }
static int main_density(int argc,char**argv){ int L=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10); unsigned char*lo=calloc(LIM+3,1); int X[40],Y[40];
  for(n=1;n<=L;n++){ u64 W=1ULL<<n; for(u64 wx=0;wx<W;wx++) for(u64 wy=wx;wy<W;wy++){ for(int i=0;i<n;i++){X[i]=wx>>i&1;Y[i]=wy>>i&1;}
        u64 h1=count_order(X,Y,1)+count_order(Y,X,0), h0=count_order(X,Y,0)+count_order(Y,X,1); u64 a=h0<h1?h0:h1, b=h0<h1?h1:h0; if(b-a==2 && a<=LIM) lo[a]=1; } 
    u64 c=0, firstmiss=0; for(u64 N=1;N<=LIM;N+=2){ if(lo[N]) c++; else if(!firstmiss && N>1) firstmiss=N; }
    printf("spine<=%d: odd N<=%llu that are lower members of a hub-arc unit pair: %llu of %llu; first miss %llu\n",n,LIM,c,(LIM+1)/2,firstmiss); fflush(stdout); }
  int k=0; printf("misses below 400:"); for(u64 N=1;N<400;N+=2) if(!lo[N]){ printf(" %llu",N); if(++k>60) break; } printf("\n"); return 0; }
static int main_res(int argc,char**argv){ int L=atoi(argv[1]); long long r9[9]={0}, r27[27]={0}, all9[9]={0}; long long tot=0; int X[40],Y[40]; long long big9=0; u64 maxsmall9=0;
  for(n=1;n<=L;n++){ u64 W=1ULL<<n; for(u64 wx=0;wx<W;wx++) for(u64 wy=0;wy<W;wy++){ for(int i=0;i<n;i++){X[i]=wx>>i&1;Y[i]=wy>>i&1;}
        u64 h1=count_order(X,Y,1)+count_order(Y,X,0), h0=count_order(X,Y,0)+count_order(Y,X,1); all9[h0%9]++;
        u64 a=h0<h1?h0:h1, b=h0<h1?h1:h0; if(b-a!=2) continue; tot++; r9[a%9]++; r27[a%27]++; if(a%9==0){ if(a>100) big9++; if(a>maxsmall9) maxsmall9=a; } } }
  printf("unit pairs (spine<=%d): %lld; lower value mod 9:",L,tot); for(int r=0;r<9;r++) printf(" %d:%lld",r,r9[r]);
  printf("\nlower mod 27 (multiples of 9 only): 0:%lld 9:%lld 18:%lld\nlower values divisible by 9 above 100: %lld; largest lower value divisible by 9: %llu\n",r27[0],r27[9],r27[18],big9,maxsmall9);
  printf("all two-hub counts mod 9:"); for(int r=0;r<9;r++) printf(" %d:%lld",r,all9[r]); printf("\n"); return 0; }

/* Modes: pairs L PMAX | density L LIM | residues L */
int main(int argc,char**argv){ if(argc<2) return 1;
  if(!strcmp(argv[1],"pairs")) return main_pairs(argc-1,argv+1);
  if(!strcmp(argv[1],"density")) return main_density(argc-1,argv+1);
  if(!strcmp(argv[1],"residues")) return main_res(argc-1,argv+1); return 1; }
