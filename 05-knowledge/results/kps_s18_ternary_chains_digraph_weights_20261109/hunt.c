/* hill-climb search for a tournament with I(Omega, X) = target (alpha-sum m = (target-1)/X) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n, A[12][12];
static long long cyc[1<<12], F[1<<12], P[1<<12][12];
static long long Ix(long long x){
  int N=1<<n;
  for(int S=0;S<N;S++) cyc[S]=0;
  for(int s=0;s<n;s++){int base=1<<s;
    for(int S=0;S<N;S++) if((S&base)&&!(S&(base-1))) for(int v=0;v<n;v++) P[S][v]=0;
    P[base][s]=1;
    for(int S=0;S<N;S++){ if(!(S&base)||(S&(base-1))) continue;
      for(int v=0;v<n;v++){ long long c=P[S][v]; if(!c) continue;
        for(int w=s+1;w<n;w++) if(!(S>>w&1)&&A[v][w]) P[S|1<<w][w]+=c;
        if(A[v][s]&&__builtin_popcount(S)>=3) cyc[S]+=c; } } }
  F[0]=1;
  for(int S=1;S<N;S++){int m=__builtin_ctz(S); long long val=F[S&~(1<<m)]; int rest=S&~(1<<m);
    for(int T=rest;;T=(T-1)&rest){int U=T|(1<<m); int pc=__builtin_popcount(U);
      if(pc%2==1&&pc>=3&&cyc[U]) val+=x*cyc[U]*F[S&~U]; if(T==0)break;}
    F[S]=val;}
  return F[N-1];
}
static unsigned long long r=88172645463325252ULL; static unsigned long long xr(void){r^=r<<13;r^=r>>7;r^=r<<17;return r;}
int main(int argc,char**argv){ long long X=atoll(argv[1]), target=atoll(argv[2]); int nmin=atoi(argv[3]), nmax=atoi(argv[4]); long long iters=atoll(argv[5]); if(argc>6) r=strtoull(argv[6],0,10)|1;
  for(long long it=0; it<iters; it++){
    n=nmin+xr()%(nmax-nmin+1);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){int b=xr()&1; A[i][j]=b; A[j][i]=!b;}
    long long cur=Ix(X); long long d=llabs(cur-target);
    for(int step=0; step<400 && d; step++){
      int i=xr()%n, j=xr()%n; if(i==j) continue;
      int t=A[i][j]; A[i][j]=A[j][i]; A[j][i]=t;
      long long v=Ix(X); long long d2=llabs(v-target);
      if(d2<=d || (xr()%20==0)) { d=d2; cur=v; } else { t=A[i][j]; A[i][j]=A[j][i]; A[j][i]=t; }
    }
    if(d==0){ printf("FOUND n=%d rows ",n); for(int i=0;i<n;i++){for(int j=0;j<n;j++) printf("%d",A[i][j]); printf(" ");} printf("\n"); return 0; }
  }
  printf("not found\n"); return 0; }
