/* Cone census: all rooted tournaments (L, v), L = B + v, with response M = U_B(O) a power of two.
   Offsets c0 = H(B) - W(O) (THM-5126 (1')); then c0 + 2^(log2 M + t) is a path count for every t >= 1.
   usage: cone_census nb ext rmax [c] < reps      (optional odd c: roots with response c 2^j instead of 2^j)      (nb = order of the listed tournaments; ext = 1: use all one-vertex
   extensions of the listed tournaments as B (order nb+1), ext = 0: the listed tournaments themselves)
   prints for every odd offset r with |r| <= rmax the least exponent j (M = 2^j) and one witness (out-masks of B, O). */
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
#define RMAX 20001
static int bestj[2*RMAX]; static int bestB[2*RMAX][10]; static int bestO[2*RMAX]; static int bestm[2*RMAX];
static long nroots=0; static u64 codd=1;
static void scan(int rmax){
    int N=1<<m; static ll wsum[1<<10]; static u64 usum[1<<10];
    wsum[0]=0; usum[0]=0;
    for(int O=1;O<N;O++){ int low=__builtin_ctz(O); int rest=O&(O-1); ll w=wsum[rest]; u64 u=usum[rest]+U[low];
        int t=rest; while(t){ int b=t&-t; int f=__builtin_ctz(b); w+=(ll)(W[low][f]+W[f][low]); t^=b; }
        wsum[O]=w; usum[O]=u;
        if(u && (u%codd)==0 && !((u/codd)&((u/codd)-1))){ ll c0=(ll)Hval-w; if(c0>-rmax&&c0<rmax){ int j=__builtin_ctzll(u); int idx=(int)(c0+RMAX); nroots++;
            if(j<bestj[idx]){ bestj[idx]=j; bestO[idx]=O; bestm[idx]=m; for(int i=0;i<m;i++) bestB[idx][i]=out[i]; } } } }
}
int main(int argc,char**argv){
    int nb=atoi(argv[1]), ext=atoi(argv[2]), rmax=atoi(argv[3]); if(rmax>=RMAX) rmax=RMAX-1;
    if(argc>4) codd=strtoull(argv[4],0,10);
    for(int i=0;i<2*RMAX;i++) bestj[i]=99;
    char line[256]; int base[10]; long nB=0;
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<nb;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} base[i]=(int)v; p=e; } if(!ok) continue;
        if(!ext){ m=nb; for(int i=0;i<m;i++) out[i]=base[i]; compute(); scan(rmax); nB++; }
        else { m=nb+1; for(int X=0;X<(1<<nb);X++){ for(int i=0;i<nb;i++){ out[i]=base[i]; if(!((X>>i)&1)) out[i]|=1<<nb; } out[nb]=X; compute(); scan(rmax); nB++; } }
    }
    printf("# bases of order %d: %ld ; roots with response %llu 2^j found: %ld\n",m,nB,codd,nroots);
    for(int r=-rmax+1;r<rmax;r++){ int idx=r+RMAX; if(bestj[idx]<99){ printf("r=%d j=%d m=%d O=%d c=%llu B:",r,bestj[idx],bestm[idx],bestO[idx],codd); for(int i=0;i<bestm[idx];i++) printf(" %d",bestB[idx][i]); printf("\n"); } }
    return 0;
}
