/* One head over a transitive tape with one small strong block.
   Base B(i, C, j) = TT_i -> C -> TT_j (ordinal sum), C a fixed strong tournament; head x with out-set X.
   For every placement i (order n = i + |C| + j fixed) all 2^n extension counts H(B + x) are computed by
       H(B + x) = H(B) + sum_{y in X} (U(y) - E(y)) - sum_{pairs e,f in X} (W_ef + W_fe);
   reported: per placement the longest run of consecutive odd counts; for the union over placements the longest
   run and the maximal fully covered intervals.
   usage: tape_block n blocktype      blocktype: 0 = none (pure tape), 3 = C3, 4 = strong 4-tournament, 5 = 1010-fan (H=9) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long ll;
#define MAXN 17
static int m; static int out[MAXN], inn[MAXN];
static u64 *E, *St, *Hs; static u64 Hval, U[MAXN], W[MAXN][MAXN];
#define EE(S,v) E[(size_t)(S)*MAXN+(v)]
#define SS(S,v) St[(size_t)(S)*MAXN+(v)]
static void compute(void){
    int N=1<<m, full=N-1;
    memset(E,0,sizeof(u64)*(size_t)N*MAXN); memset(St,0,sizeof(u64)*(size_t)N*MAXN);
    for(int v=0;v<m;v++){ EE(1<<v,v)=1; SS(1<<v,v)=1; inn[v]=0; }
    for(int v=0;v<m;v++) for(int w=0;w<m;w++) if((out[w]>>v)&1) inn[v]|=1<<w;
    for(int S=1;S<N;S++){ for(int v=0;v<m;v++){ u64 c=EE(S,v); if(c){ int t=out[v]&~S; while(t){int b=t&-t; EE(S|b,__builtin_ctz(b))+=c; t^=b;} }
                                             u64 d=SS(S,v); if(d){ int t=inn[v]&~S; while(t){int b=t&-t; SS(S|b,__builtin_ctz(b))+=d; t^=b;} } } }
    for(int S=0;S<N;S++){ u64 h=0; for(int v=0;v<m;v++) h+=EE(S,v); Hs[S]=h; } Hs[0]=1; Hval=Hs[full];
    memset(U,0,sizeof(U)); memset(W,0,sizeof(W));
    for(int S=0;S<N;S++){ int C=full^S; if(!C) continue; u64 hs=Hs[S];
        for(int y=0;y<m;y++) if((C>>y)&1) U[y]+=hs*SS(C,y);
        if(S) for(int e=0;e<m;e++) if((S>>e)&1){ u64 te=EE(S,e); if(!te) continue; for(int f=0;f<m;f++) if((C>>f)&1) W[e][f]+=te*SS(C,f); } }
}
static int cmpll(const void*a,const void*b){ ll x=*(const ll*)a,y=*(const ll*)b; return (x>y)-(x<y); }
int main(int argc,char**argv){
    int n=atoi(argv[1]), bt=atoi(argv[2]); m=n; int N=1<<n, full=N-1;
    E=malloc(sizeof(u64)*(size_t)N*MAXN); St=malloc(sizeof(u64)*(size_t)N*MAXN); Hs=malloc(sizeof(u64)*N);
    ll *val=malloc(sizeof(ll)*N), *srt=malloc(sizeof(ll)*N);
    int bs = bt==0?0:(bt==3?3:(bt==4?4:5));
    /* block arcs inside (local indices 0..bs-1), default forward i->j for i<j, with listed reversed arcs */
    int rev[4][2]; int nrev=0;
    if(bt==3){ rev[0][0]=0; rev[0][1]=2; nrev=1; }                 /* 0->1->2->0 */
    if(bt==4){ rev[0][0]=0; rev[0][1]=3; nrev=1; }                 /* TT4 with arc (0,3) reversed: strong, H = 5 */
    if(bt==5){ rev[0][0]=0; rev[0][1]=2; rev[1][0]=2; rev[1][1]=4; nrev=2; }  /* two triangles sharing a vertex, H = 9 */
    ll LIM = 1LL<<(n+3); unsigned char *cov=calloc((size_t)LIM,1);
    for(int i=0;i+bs<=n;i++){
        for(int a=0;a<n;a++){ out[a]=0; for(int b=a+1;b<n;b++) out[a]|=1<<b; }
        for(int r=0;r<nrev;r++){ int a=i+rev[r][0], b=i+rev[r][1]; out[a]&=~(1<<b); out[b]|=1<<a; }
        compute(); ll d=(ll)Hval; val[0]=d;
        for(int X=1;X<N;X++){ int low=__builtin_ctz(X); int rest=X&(X-1); ll v=val[rest]+(ll)U[low]-(ll)EE(full,low);
            int t=rest; while(t){ int b=t&-t; int f=__builtin_ctz(b); v-=(ll)(W[low][f]+W[f][low]); t^=b; } val[X]=v; if(v>0&&v<LIM) cov[v]=1; }
        memcpy(srt,val,sizeof(ll)*N); qsort(srt,N,sizeof(ll),cmpll);
        int best=1,cur=1; ll bsv=srt[0],cs=srt[0]; long distinct=1;
        for(int k=1;k<N;k++){ if(srt[k]==srt[k-1]) continue; distinct++; if(srt[k]==srt[k-1]+2) cur++; else { cur=1; cs=srt[k]; } if(cur>best){ best=cur; bsv=cs; } }
        printf("n=%d block=%d placement i=%d: H(base)=%lld ; distinct counts %ld of %d ; max count %lld ; longest run %d: %lld..%lld\n",n,bt,i,d,distinct,N,srt[N-1],best,bsv,bsv+2*(best-1));
        if(bs==0) break;
    }
    /* union over placements: maximal covered intervals of odd numbers */
    int best=0,cur=0; ll bstart=0,cstart=0; long tot=0;
    for(ll v=1;v<LIM;v+=2){ if(cov[v]){ tot++; if(!cur) cstart=v; cur++; if(cur>best){ best=cur; bstart=cstart; } } else cur=0; }
    printf("UNION n=%d block=%d: %ld odd counts below %lld ; longest fully covered interval: %d odd numbers %lld..%lld\n",n,bt,tot,LIM,best,bstart,bstart+2*(ll)(best-1));
    return 0;
}
