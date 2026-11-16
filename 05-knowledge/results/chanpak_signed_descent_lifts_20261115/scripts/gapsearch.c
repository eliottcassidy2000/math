/* Local search for small nontrivial ratios H(S)/H(S - x): is (3 H(S-x) + 1)/2 a lower bound for H(S) whenever
   x lies on a cycle?  Tournament on n vertices, x = vertex n-1.  Minimise f = 2 H(S) - 3 H(S-x) over tournaments
   with H(S) > H(S-x) by random arc flips (simulated annealing).  usage: gapsearch n iters seed */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long ll;
static int n; static int adj[16];
static ll dp[1<<15][15];
static ll ham(int nv){ int N=1<<nv; for(int S=1;S<N;S++) memset(dp[S],0,sizeof(ll)*nv); for(int v=0;v<nv;v++) dp[1<<v][v]=1;
    int msk=N-1;
    for(int S=1;S<N;S++) for(int v=0;v<nv;v++){ ll c=dp[S][v]; if(!c) continue; int t=adj[v]&msk&~S; while(t){int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b;} }
    ll h=0; for(int v=0;v<nv;v++) h+=dp[N-1][v]; return h; }
static double rnd(void){ return rand()/(RAND_MAX+1.0); }
int main(int argc,char**argv){ n=atoi(argv[1]); long iters=atol(argv[2]); srand(atoi(argv[3]));
    ll best=1LL<<60; 
    for(int restart=0; restart<40; restart++){
        memset(adj,0,sizeof(adj));
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(rnd()<0.5) adj[i]|=1<<j; else adj[j]|=1<<i; }
        ll hS=ham(n), hT=ham(n-1); double cur=(hS>hT)?(double)(2*hS-3*hT)/hT:1e9; 
        for(long it=0; it<iters/40; it++){
            int i=rand()%n, j=rand()%n; if(i==j) continue;
            adj[i]^=1<<j; adj[j]^=1<<i;
            ll a=ham(n), b=ham(n-1); double val=(a>b)?(double)(2*a-3*b)/b:1e9;
            double temp=0.05*(1.0-(double)it/(iters/40));
            if(val<=cur || rnd()<exp((cur-val)/ (temp+1e-9))){ cur=val; hS=a; hT=b;
                if(a>b && 2*a-3*b<best){ best=2*a-3*b; printf("n=%d: H(S)=%lld H(S-x)=%lld  2H(S)-3H(S-x)=%lld  ratio %.5f\n",n,a,b,best,(double)a/b); fflush(stdout); }
                if(a>b && 2*a-3*b<=0){ printf("COUNTEREXAMPLE adj:"); for(int v=0;v<n;v++) printf(" %d",adj[v]); printf("\n"); return 0; } }
            else { adj[i]^=1<<j; adj[j]^=1<<i; }
        }
    }
    printf("n=%d: least 2H(S)-3H(S-x) found: %lld\n",n,best); return 0; }
