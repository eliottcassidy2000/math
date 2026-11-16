/* Synchronised banks of one head over  (upper tape) -> C -> (lower tape).
   For a block C and a head-contact pattern X on C the count is
     H = K(tapes) + 2 [ cut_C(X) + a_o(X)(alpha+u1) + a_e(X) alpha + b_o(X)(beta+z1) + b_e(X) beta ],
   a_par(X) = sum_{s not in X} lam_par[s],  b_par(X) = sum_{s in X} kap_par[s]  (see master_formula_check.py).
   Patterns with the same signature differ by the constant 2 (cut_C(X) - cut_C(X')) for EVERY pair of tape words.
   A synchronised bank of length L = L patterns of equal signature whose cut values are L consecutive integers.
   Variants: V1 key a_o+a_e (upper tape only, lowest tape vertex not a contact); V2 key (a_o,a_e) (upper tape, any word);
             V3 key (a_o+a_e, b_o+b_e) (two tapes, u1 = z1 = 0);  V4 key (a_o,a_e,b_o,b_e) (two tapes, any words).
   usage: sync_bank m < masks */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned long long u64;
#define MM 10
static int m, out[MM];
static ll h[1<<MM]; static ll Ee[1<<MM][MM]; static ll Pth[1<<MM][MM];
static ll Wp[MM][MM][2], lam[MM][2], kap[MM][2], G[2];
typedef struct { ll k1,k2,k3,k4,cut; int X; } rec;
static rec R[1<<MM];
static int cmp(const void*a,const void*b){ const rec*x=a,*y=b;
    if(x->k1!=y->k1) return x->k1<y->k1?-1:1; if(x->k2!=y->k2) return x->k2<y->k2?-1:1;
    if(x->k3!=y->k3) return x->k3<y->k3?-1:1; if(x->k4!=y->k4) return x->k4<y->k4?-1:1;
    if(x->cut!=y->cut) return x->cut<y->cut?-1:1; return 0; }
static int bestL[5]; static int bestOut[5][MM]; static ll bestKey[5][4], bestStart[5]; static int bestPat[5][64]; static long hist[5][64]; static ll bestH[5];
static void block(void){
    int N=1<<m, full=N-1;
    for(int S=0;S<N;S++) memset(Ee[S],0,sizeof(ll)*m);
    for(int v=0;v<m;v++) Ee[1<<v][v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<m;v++){ ll c=Ee[S][v]; if(c){ int t=out[v]&~S; while(t){ int b=t&-t; Ee[S|b][__builtin_ctz(b)]+=c; t^=b; } } }
    for(int S=0;S<N;S++){ ll s=0; for(int v=0;v<m;v++) s+=Ee[S][v]; h[S]=s; } h[0]=1;
    memset(Wp,0,sizeof(Wp)); G[0]=h[full]; G[1]=0;
    for(int s=0;s<m;s++){
        for(int S=0;S<N;S++) memset(Pth[S],0,sizeof(ll)*m);
        Pth[1<<s][s]=1;
        for(int S=1;S<N;S++) if((S>>s)&1) for(int v=0;v<m;v++){ ll c=Pth[S][v]; if(!c) continue;
            int par=__builtin_popcount(S)&1; ll w=c*h[full^S]; Wp[s][v][par]+=w; G[par]+=w;
            int t=out[v]&~S; while(t){ int b=t&-t; Pth[S|b][__builtin_ctz(b)]+=c; t^=b; } }
    }
    for(int s=0;s<m;s++) for(int p=0;p<2;p++){ ll a=0,b=0; for(int t=0;t<m;t++){ a+=Wp[t][s][p]; b+=Wp[s][t][p]; } lam[s][p]=a; kap[s][p]=b; }
    static ll cut[1<<MM], ao[1<<MM], ae[1<<MM], bo[1<<MM], be[1<<MM];
    for(int X=0;X<N;X++){ ll c=0,xo=0,xe=0,yo=0,ye=0;
        for(int a=0;a<m;a++){ if((X>>a)&1){ yo+=kap[a][1]; ye+=kap[a][0]; for(int b=0;b<m;b++) if(!((X>>b)&1)) c+=Wp[a][b][0]; } else { xo+=lam[a][1]; xe+=lam[a][0]; } }
        cut[X]=c; ao[X]=xo; ae[X]=xe; bo[X]=yo; be[X]=ye; }
    for(int V=1;V<=4;V++){
        for(int X=0;X<N;X++){ R[X].X=X; R[X].cut=cut[X];
            if(V==1){ R[X].k1=ao[X]+ae[X]; R[X].k2=R[X].k3=R[X].k4=0; }
            if(V==2){ R[X].k1=ao[X]; R[X].k2=ae[X]; R[X].k3=R[X].k4=0; }
            if(V==3){ R[X].k1=ao[X]+ae[X]; R[X].k2=bo[X]+be[X]; R[X].k3=R[X].k4=0; }
            if(V==4){ R[X].k1=ao[X]; R[X].k2=ae[X]; R[X].k3=bo[X]; R[X].k4=be[X]; } }
        qsort(R,N,sizeof(rec),cmp);
        int best=1,cur=1,bi=0,ci=0;
        for(int i=1;i<N;i++){ int same=(R[i].k1==R[i-1].k1&&R[i].k2==R[i-1].k2&&R[i].k3==R[i-1].k3&&R[i].k4==R[i-1].k4);
            if(same&&R[i].cut==R[i-1].cut) continue;
            if(same&&R[i].cut==R[i-1].cut+1) cur++; else { cur=1; ci=i; }
            if(cur>best){ best=cur; bi=ci; } }
        hist[V][best<63?best:63]++;
        if(best>bestL[V]){ bestL[V]=best; for(int i=0;i<m;i++) bestOut[V][i]=out[i]; bestKey[V][0]=R[bi].k1; bestKey[V][1]=R[bi].k2; bestKey[V][2]=R[bi].k3; bestKey[V][3]=R[bi].k4; bestStart[V]=R[bi].cut; bestH[V]=h[full];
            int k=0; ll last=R[bi].cut-1; for(int i=bi;i<N&&k<best&&k<64;i++){ if(R[i].cut==last) continue; if(R[i].cut!=last+1) break; bestPat[V][k++]=R[i].X; last=R[i].cut; } }
    }
}
int main(int argc,char**argv){
    m=atoi(argv[1]); char line[256]; long nb=0;
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<m;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} out[i]=(int)v; p=e; } if(!ok) continue; block(); nb++; }
    printf("# blocks of order %d: %ld\n",m,nb);
    const char*nm[5]={"","V1 upper tape, u1=0","V2 upper tape, any word","V3 two tapes, u1=z1=0","V4 two tapes, any words"};
    for(int V=1;V<=4;V++){ printf("%s: longest synchronised bank %d (cut values %lld..%lld; H(C) = %lld; key %lld %lld %lld %lld); block:",nm[V],bestL[V],bestStart[V],bestStart[V]+bestL[V]-1,bestH[V],bestKey[V][0],bestKey[V][1],bestKey[V][2],bestKey[V][3]);
        for(int i=0;i<m;i++) printf(" %d",bestOut[V][i]); printf(" ; patterns:"); for(int k=0;k<bestL[V]&&k<64;k++) printf(" %d",bestPat[V][k]); printf("\n");
        printf("   histogram (bank length: blocks):"); for(int i=1;i<64;i++) if(hist[V][i]) printf(" %d:%ld",i,hist[V][i]); printf("\n"); }
    return 0;
}
