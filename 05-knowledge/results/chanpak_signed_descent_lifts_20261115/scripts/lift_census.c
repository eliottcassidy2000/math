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
static int DMAX; static unsigned char *occ; static unsigned char *hit[4]; static int (*wit[4])[11];
static long nB=0;
static void scan(void){
    int N=1<<m, full=N-1; static ll val[1<<10];
    ll d=(ll)Hval; if(d>DMAX) return; occ[d]=1;
    ll tgt[4]; tgt[0]=2*d+1; tgt[1]=2*d-1; tgt[2]=(d%6==3)?(4*d-3)/3:-1; tgt[3]=(d%6==5)?(4*d+1)/3:-1;
    val[0]=d;
    for(int X=1;X<N;X++){ int low=__builtin_ctz(X); int rest=X&(X-1); ll v=val[rest]+(ll)U[low]-(ll)E[full][low];
        int t=rest; while(t){ int b=t&-t; int f=__builtin_ctz(b); v-=(ll)(W[low][f]+W[f][low]); t^=b; }
        val[X]=v;
        for(int k=0;k<4;k++) if(v==tgt[k] && !hit[k][d]){ hit[k][d]=1; for(int i=0;i<m;i++) wit[k][d][i]=out[i]; wit[k][d][m]=-1; wit[k][d][10]=X; } }
}
int main(int argc,char**argv){
    int nb=atoi(argv[1]), ext=atoi(argv[2]); DMAX=atoi(argv[3]);
    occ=calloc(DMAX+1,1); for(int k=0;k<4;k++){ hit[k]=calloc(DMAX+1,1); wit[k]=calloc(DMAX+1,sizeof(int[11])); }
    char line[256]; int base[10];
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<nb;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} base[i]=(int)v; p=e; } if(!ok) continue;
        if(!ext){ m=nb; for(int i=0;i<m;i++) out[i]=base[i]; compute(); scan(); nB++; }
        else { m=nb+1; for(int X=0;X<(1<<nb);X++){ for(int i=0;i<nb;i++){ out[i]=base[i]; if(!((X>>i)&1)) out[i]|=1<<nb; } out[nb]=X; compute(); scan(); nB++; } }
    }
    const char*nm[4]={"b0: d -> 2d+1","b1: d -> 2d-1","b2: d = 3 (mod 6) -> (4d-3)/3","b3: d = 5 (mod 6) -> (4d+1)/3"};
    printf("# bases of order %d: %ld ; DMAX %d\n",m,nB,DMAX);
    for(int k=0;k<4;k++){ long dom=0,lift=0; printf("## %s\n",nm[k]);
        for(int d=1;d<=DMAX;d+=2){ if(!occ[d]) continue; if(k==2&&d%6!=3) continue; if(k==3&&d%6!=5) continue; dom++; if(hit[k][d]) lift++; }
        printf("domain (d occurring as H of a base): %ld ; lifted by one vertex: %ld ; not lifted:",dom,lift);
        int shown=0; for(int d=1;d<=DMAX;d+=2){ if(!occ[d]) continue; if(k==2&&d%6!=3) continue; if(k==3&&d%6!=5) continue; if(!hit[k][d]){ if(shown<60) printf(" %d",d); shown++; } }
        printf("  (%d in all)\n",shown);
        for(int d=1;d<=DMAX;d+=2) if(hit[k][d]){ printf("W %d %d :",k,d); for(int i=0;i<10&&wit[k][d][i]>=0;i++) printf(" %d",wit[k][d][i]); printf(" X=%d\n",wit[k][d][10]); }
    }
    return 0;
}
