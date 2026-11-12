/* Auditor's literal recount. stdin lines: "c r j m O mask_0 ... mask_{m-1}" (mask i bit k: i -> k in B).
   Builds L = B + v (v -> y iff y in O, else y -> v) and L o 0^t (chain z1->..->zt, B beats all z, all z beat v),
   counts Hamiltonian paths by brute DP over (set, last vertex) on an adjacency matrix,
   and prints "c r j t H expected ok". */
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll dp[1<<15][15];
static int A[15][15];
static ll count(int N){
    for(int S=0;S<(1<<N);S++) for(int a=0;a<N;a++) dp[S][a]=0;
    for(int a=0;a<N;a++) dp[1<<a][a]=1;
    for(int S=1;S<(1<<N);S++) for(int a=0;a<N;a++){ ll x=dp[S][a]; if(!x) continue;
        for(int b=0;b<N;b++) if(!((S>>b)&1) && A[a][b]) dp[S|(1<<b)][b]+=x; }
    ll h=0; for(int a=0;a<N;a++) h+=dp[(1<<N)-1][a]; return h;
}
int main(void){
    ll c,r; int j,m,O; long bad=0,tot=0;
    while(scanf("%lld %lld %d %d %d",&c,&r,&j,&m,&O)==5){
        int mk[12]; for(int i=0;i<m;i++) scanf("%d",&mk[i]);
        for(int t=1;t<=3;t++){
            int N=m+1+t, v=m;
            for(int a=0;a<N;a++) for(int b=0;b<N;b++) A[a][b]=0;
            for(int a=0;a<m;a++) for(int b=0;b<m;b++) if(a!=b && ((mk[a]>>b)&1)) A[a][b]=1;
            /* sanity: tournament on B */
            for(int a=0;a<m;a++) for(int b=a+1;b<m;b++) if(A[a][b]+A[b][a]!=1){ printf("NOT A TOURNAMENT\n"); return 2; }
            for(int y=0;y<m;y++){ if((O>>y)&1) A[v][y]=1; else A[y][v]=1; }
            for(int z=0;z<t;z++){ int Z=m+1+z; for(int y=0;y<m;y++) A[y][Z]=1; A[Z][v]=1; for(int z2=z+1;z2<t;z2++) A[Z][m+1+z2]=1; }
            ll h=count(N), want=r+c*(1LL<<(j+t)); tot++;
            if(h!=want) bad++;
            printf("%lld %lld %d %d %lld %lld %s\n",c,r,j,t,h,want,h==want?"ok":"FAIL");
        }
    }
    fprintf(stderr,"checks %ld failures %ld\n",tot,bad);
    return 0;
}
