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

/* Runs of consecutive odd counts among the one-vertex extensions of a single base tournament B:
   S(B) = { H(B + x) : out-set X of x }.  For every base the longest run v, v+2, ..., v+2(L-1) inside S(B).
   usage: ext_runs nb ext TOPK < classes   (ext = 0: listed tournaments as bases; 1: their one-vertex extensions) */
typedef struct { int L; ll start; ll d; int masks[10]; int m; } rec;
static rec top[64]; static rec topr[64]; static double topq[64]; static int TOPK=20; static int bestL[4096]; static ll bestS[4096]; static long nB=0; static long hist[512];
static int cmpll(const void*a,const void*b){ ll x=*(const ll*)a,y=*(const ll*)b; return (x>y)-(x<y); }
static void scan(void){
    int N=1<<m, full=N-1; static ll val[1<<10]; static ll srt[1<<10];
    ll d=(ll)Hval; val[0]=d;
    for(int X=1;X<N;X++){ int low=__builtin_ctz(X); int rest=X&(X-1); ll v=val[rest]+(ll)U[low]-(ll)E[full][low];
        int t=rest; while(t){ int b=t&-t; int f=__builtin_ctz(b); v-=(ll)(W[low][f]+W[f][low]); t^=b; } val[X]=v; }
    memcpy(srt,val,sizeof(ll)*N); qsort(srt,N,sizeof(ll),cmpll);
    int best=1,cur=1; ll bs=srt[0],cs=srt[0];
    for(int i=1;i<N;i++){ if(srt[i]==srt[i-1]) continue; if(srt[i]==srt[i-1]+2){ cur++; } else { cur=1; cs=srt[i]; } if(cur>best){ best=cur; bs=cs; } }
    hist[best<511?best:511]++;
    if(d<4096 && best>bestL[d]){ bestL[d]=best; bestS[d]=bs; }
    { /* best ratio run */ int c2=1; ll s2=srt[0]; double bq=0; ll bqs=0; int bql=0;
      for(int i=1;i<=N;i++){ int brk=(i==N)||(srt[i]!=srt[i-1]&&srt[i]!=srt[i-1]+2); if(i<N&&srt[i]==srt[i-1]) continue; if(!brk){ c2++; continue; }
          if(s2>=9){ double q=(double)(s2+2*c2)/(double)s2; if(q>bq){ bq=q; bqs=s2; bql=c2; } } if(i<N){ c2=1; s2=srt[i]; } }
      int pos2=-1; for(int k=0;k<TOPK;k++) if(bq>topq[k]){ pos2=k; break; }
      if(pos2>=0){ for(int k=TOPK-1;k>pos2;k--){ topr[k]=topr[k-1]; topq[k]=topq[k-1]; } topq[pos2]=bq; topr[pos2].L=bql; topr[pos2].start=bqs; topr[pos2].d=d; topr[pos2].m=m; for(int i=0;i<m;i++) topr[pos2].masks[i]=out[i]; } }
    /* insert into top list */
    int pos=-1; for(int k=0;k<TOPK;k++) if(best>top[k].L){ pos=k; break; }
    if(pos>=0){ for(int k=TOPK-1;k>pos;k--) top[k]=top[k-1]; top[pos].L=best; top[pos].start=bs; top[pos].d=d; top[pos].m=m; for(int i=0;i<m;i++) top[pos].masks[i]=out[i]; }
}
int main(int argc,char**argv){
    int nb=atoi(argv[1]), ext=atoi(argv[2]); if(argc>3) TOPK=atoi(argv[3]);
    char line[256]; int base[10];
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<nb;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} base[i]=(int)v; p=e; } if(!ok) continue;
        if(!ext){ m=nb; for(int i=0;i<m;i++) out[i]=base[i]; compute(); scan(); nB++; }
        else { m=nb+1; for(int X=0;X<(1<<nb);X++){ for(int i=0;i<nb;i++){ out[i]=base[i]; if(!((X>>i)&1)) out[i]|=1<<nb; } out[nb]=X; compute(); scan(); nB++; } }
    }
    printf("# bases of order %d: %ld\n",m,nB);
    printf("# histogram of the longest run (run length: number of bases):"); for(int i=1;i<512;i++) if(hist[i]) printf(" %d:%ld",i,hist[i]); printf("\n");
    for(int k=0;k<TOPK;k++) if(top[k].L){ printf("TOP run %d from %lld to %lld ; H(base) = %lld ; base:",top[k].L,top[k].start,top[k].start+2*(top[k].L-1),top[k].d); for(int i=0;i<top[k].m;i++) printf(" %d",top[k].masks[i]); printf("\n"); }
    for(int k=0;k<TOPK;k++) if(topr[k].L){ printf("RATIO %.4f run %d from %lld to %lld ; H(base) = %lld ; base:",topq[k],topr[k].L,topr[k].start,topr[k].start+2*(topr[k].L-1),topr[k].d); for(int i=0;i<topr[k].m;i++) printf(" %d",topr[k].masks[i]); printf("\n"); }
    printf("PERD"); for(int d=1;d<4096;d+=2) if(bestL[d]) printf(" %d:%d@%lld",d,bestL[d],bestS[d]); printf("\n");
    return 0;
}
