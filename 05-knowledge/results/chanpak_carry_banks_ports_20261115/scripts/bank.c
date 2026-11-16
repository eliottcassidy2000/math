/* One-vertex lifts of the three branches of the signed descent D (six_tutte_collatz_20261116/COUNTS.md):
     D(h) = (h+1)/2 (h = 1 mod 4), (3h+3)/4 (h = 3 mod 8), (3h-1)/4 (h = 7 mod 8).
   Inverse branches:  b1: d -> 2d-1 ;  b2: d = 3 mod 6 -> (4d-3)/3 ;  b3: d = 5 mod 6 -> (4d+1)/3 ;  and b0: d -> 2d+1
   (the second branch of the two-sign halving descent).
   For every base tournament B (H(B) = d) and every out-set X of a new vertex x:
       H(B + x) = H(B) + sum_{y in X} (U_B(y) - E_B(y)) - sum_{pairs e,f in X} W_ef .
   A branch is LIFTED at d by one vertex if some base with H = d has an extension with the branch target.
   usage: lift_census nb ext DMAX < classes   (ext = 0: the listed tournaments as bases; ext = 1: all their one-vertex
   extensions as bases).  Prints, for each branch, the d <= DMAX in its domain that occur as H of a base, and those
   among them without a lift; one witness per lifted d is kept. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef long long ll;
static int m; static int out[10], inn[10];
static u64 E[1<<10][10], St[1<<10][10], Hs[1<<10];
static u64 Hval, U[10], W[10][10];
static void compute(void){
    int N=1<<m, full=N-1;
    for(int S=0;S<N;S++){ memset(E[S],0,sizeof(u64)*m); memset(St[S],0,sizeof(u64)*m); }
    for(int v=0;v<m;v++){ E[1<<v][v]=1; St[1<<v][v]=1; inn[v]=0; }
    for(int v=0;v<m;v++) for(int w=0;w<m;w++) if((out[w]>>v)&1) inn[v]|=1<<w;
    for(int S=1;S<N;S++){ for(int v=0;v<m;v++){ u64 c=E[S][v]; if(c){ int t=out[v]&~S; while(t){int b=t&-t; E[S|b][__builtin_ctz(b)]+=c; t^=b;} }
                                             u64 d=St[S][v]; if(d){ int t=inn[v]&~S; while(t){int b=t&-t; St[S|b][__builtin_ctz(b)]+=d; t^=b;} } } }
    for(int S=0;S<N;S++){ u64 h=0; for(int v=0;v<m;v++) h+=E[S][v]; Hs[S]=h; } Hs[0]=1; Hval=Hs[full];
    memset(U,0,sizeof(U)); memset(W,0,sizeof(W));
    for(int S=0;S<N;S++){ int C=full^S; if(!C) continue; u64 hs=Hs[S];
        for(int y=0;y<m;y++) if((C>>y)&1) U[y]+=hs*St[C][y];
        if(S) for(int e=0;e<m;e++) if((S>>e)&1){ u64 te=E[S][e]; if(!te) continue; for(int f=0;f<m;f++) if((C>>f)&1) W[e][f]+=te*St[C][f]; } }
}

/* Complete banks.  For a base tournament B and a new vertex x with out-set X:  H(B+x) = H(B) + 2 cut_W(X).
   A COMPLETE k-BANK is a context X0 and a set S of k vertices (disjoint from X0) such that the 2^k out-sets
   X0 + Y, Y a subset of S, give 2^k consecutive odd counts.  For every base: the largest k, with one witness.
   usage: bank m TOPK < masks      (lines of m out-masks) */
static long nB=0; static long hist[12]; static int TOPK=6; static int shown[12];
static void scan(void){
    int N=1<<m, full=N-1; static ll val[1<<10];
    ll d=(ll)Hval; val[0]=d;
    for(int X=1;X<N;X++){ int low=__builtin_ctz(X); int rest=X&(X-1); ll v=val[rest]+(ll)U[low]-(ll)E[full][low];
        int t=rest; while(t){ int b=t&-t; int f=__builtin_ctz(b); v-=(ll)(W[low][f]+W[f][low]); t^=b; } val[X]=v; }
    int bestk=0, bS=0, bX=0; ll bmin=0;
    for(int S=1;S<N;S++){ int k=__builtin_popcount(S); if(k<2||k<=bestk) continue; ll need=2*((1LL<<k)-1);
        int comp=full^S;
        for(int X0=comp;;X0=(X0-1)&comp){
            ll mn=val[X0], mx=val[X0]; int ok=1;
            for(int Y=S;Y;Y=(Y-1)&S){ ll v=val[X0|Y]; if(v<mn) mn=v; if(v>mx) mx=v; if(mx-mn>need){ ok=0; break; } }
            if(ok && mx-mn==need){ /* distinctness */ unsigned char seen[1024]; memset(seen,0,(size_t)(1<<k)); int dist=1;
                ll v0=val[X0]; int idx=(int)((v0-mn)/2); seen[idx]=1;
                for(int Y=S;Y;Y=(Y-1)&S){ int id=(int)((val[X0|Y]-mn)/2); if(seen[id]){ dist=0; break; } seen[id]=1; }
                if(dist){ bestk=k; bS=S; bX=X0; bmin=mn; break; } }
            if(X0==0) break; }
    }
    hist[bestk]++;
    if(bestk>=2 && shown[bestk]<TOPK){ shown[bestk]++; printf("BANK k=%d counts %lld..%lld ; H(base) = %lld ; S=%d X0=%d ; base:",bestk,bmin,bmin+2*((1LL<<bestk)-1),d,bS,bX); for(int i=0;i<m;i++) printf(" %d",out[i]); printf("\n"); }
}
int main(int argc,char**argv){
    m=atoi(argv[1]); if(argc>2) TOPK=atoi(argv[2]);
    char line[256];
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<m;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} out[i]=(int)v; p=e; } if(!ok) continue; compute(); scan(); nB++; }
    printf("# bases of order %d: %ld ; largest complete bank k: number of bases:",m,nB); for(int k=0;k<12;k++) if(hist[k]) printf(" %d:%ld",k,hist[k]); printf("\n");
    return 0;
}
