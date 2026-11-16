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

/* Cross-seed synchronised banks.  A seed is (block C, head pattern X on C).  For tapes above and below,
     H = h (1 + 2 N_up + 2 N_low) + 2 [ cut + a_o (alpha+u1) + a_e alpha + b_o (beta+z1) + b_e beta
                                         + G_e ((alpha+u1)(beta+z1) + alpha beta) + G_o ((alpha+u1) beta + alpha (beta+z1)) ].
   Seeds with equal invariants give counts differing by 2 (cut - cut') for every pair of tape words:
     W1 (upper tape only, u1 = 0): key (h, a_o+a_e);          W2 (upper tape only): key (h, a_o, a_e);
     W3 (two tapes, u1 = z1 = 0): key (h, G_e+G_o, a_o+a_e, b_o+b_e);   W4 (two tapes): key (h, G_e, G_o, a_o, a_e, b_o, b_e).
   usage: cross_bank m1 file1 [m2 file2 ...] */
typedef struct { ll k[7]; ll cut; int blk; int X; unsigned char m; } rec;
static rec *R; static long nR=0, capR=0; static int VAR=1;
static int nk[5]={0,2,3,4,7};
static void key(const rec*r,ll*o){ const ll*k=r->k; /* stored: h, Ge, Go, ao, ae, bo, be */
    if(VAR==1){ o[0]=k[0]; o[1]=k[3]+k[4]; } else if(VAR==2){ o[0]=k[0]; o[1]=k[3]; o[2]=k[4]; }
    else if(VAR==3){ o[0]=k[0]; o[1]=k[1]+k[2]; o[2]=k[3]+k[4]; o[3]=k[5]+k[6]; } else { for(int i=0;i<7;i++) o[i]=k[i]; } }
static int cmp(const void*a,const void*b){ ll x[7],y[7]; key(a,x); key(b,y); for(int i=0;i<nk[VAR];i++) if(x[i]!=y[i]) return x[i]<y[i]?-1:1;
    const rec*p=a,*q=b; if(p->cut!=q->cut) return p->cut<q->cut?-1:1; return 0; }
static int (*blocks)[MM]; static long nblk=0, capB=0; static ll HMAX=(1LL<<60); static long skipped=0;
static void block(void){
    int N=1<<m, full=N-1;
    for(int S=0;S<N;S++) memset(Ee[S],0,sizeof(ll)*m);
    for(int v=0;v<m;v++) Ee[1<<v][v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<m;v++){ ll c=Ee[S][v]; if(c){ int t=out[v]&~S; while(t){ int b=t&-t; Ee[S|b][__builtin_ctz(b)]+=c; t^=b; } } }
    for(int S=0;S<N;S++){ ll s=0; for(int v=0;v<m;v++) s+=Ee[S][v]; h[S]=s; } h[0]=1;
    if(h[full]>HMAX){ skipped++; return; }
    memset(Wp,0,sizeof(Wp)); G[0]=h[full]; G[1]=0;
    for(int s=0;s<m;s++){
        for(int S=0;S<N;S++) memset(Pth[S],0,sizeof(ll)*m);
        Pth[1<<s][s]=1;
        for(int S=1;S<N;S++) if((S>>s)&1) for(int v=0;v<m;v++){ ll c=Pth[S][v]; if(!c) continue;
            int par=__builtin_popcount(S)&1; ll w=c*h[full^S]; Wp[s][v][par]+=w; G[par]+=w;
            int t=out[v]&~S; while(t){ int b=t&-t; Pth[S|b][__builtin_ctz(b)]+=c; t^=b; } }
    }
    for(int s=0;s<m;s++) for(int p=0;p<2;p++){ ll a=0,b=0; for(int t=0;t<m;t++){ a+=Wp[t][s][p]; b+=Wp[s][t][p]; } lam[s][p]=a; kap[s][p]=b; }
    if(nblk==capB){ capB=capB?2*capB:1024; blocks=realloc(blocks,sizeof(int[MM])*capB); }
    for(int i=0;i<MM;i++) blocks[nblk][i]=i<m?out[i]:-1;
    for(int X=0;X<N;X++){ ll c=0,xo=0,xe=0,yo=0,ye=0;
        for(int a=0;a<m;a++){ if((X>>a)&1){ yo+=kap[a][1]; ye+=kap[a][0]; for(int b=0;b<m;b++) if(!((X>>b)&1)) c+=Wp[a][b][0]; } else { xo+=lam[a][1]; xe+=lam[a][0]; } }
        if(nR==capR){ capR=capR?2*capR:(1<<16); R=realloc(R,sizeof(rec)*capR); }
        rec*r=&R[nR++]; r->k[0]=h[full]; r->k[1]=G[0]; r->k[2]=G[1]; r->k[3]=xo; r->k[4]=xe; r->k[5]=yo; r->k[6]=ye; r->cut=c; r->blk=(int)nblk; r->X=X; r->m=(unsigned char)m; }
    nblk++;
}
int main(int argc,char**argv){
    if(getenv("HMAX")) HMAX=atoll(getenv("HMAX"));
    for(int a=1;a+1<argc;a+=2){ m=atoi(argv[a]); FILE*f=fopen(argv[a+1],"r"); char line[256];
        while(fgets(line,sizeof(line),f)){ char*p=line; int ok=1; for(int i=0;i<m;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} out[i]=(int)v; p=e; } if(!ok) continue; block(); } fclose(f); }
    printf("# seeds (block, pattern): %ld from %ld blocks (block count at most %lld; %ld blocks skipped)\n",nR,nblk,HMAX,skipped);
    if(getenv("DUMP")){ FILE*d=fopen(getenv("DUMP"),"w"); for(long i=0;i<nR;i++){ rec*r=&R[i]; fprintf(d,"%lld %lld %lld %lld %lld %lld %lld %lld %d %d %d",r->k[0],r->k[1],r->k[2],r->k[3],r->k[4],r->k[5],r->k[6],r->cut,r->blk,r->X,(int)r->m); for(int t=0;t<r->m;t++) fprintf(d," %d",blocks[r->blk][t]); fprintf(d,"\n"); } fclose(d); return 0; }
    const char*nm[5]={"","W1 upper tape, u1=0: key (h, a)","W2 upper tape: key (h, a_o, a_e)","W3 two tapes, u1=z1=0: key (h, G, a, b)","W4 two tapes: key (h, G_e, G_o, a_o, a_e, b_o, b_e)"};
    for(VAR=1;VAR<=4;VAR++){ qsort(R,nR,sizeof(rec),cmp);
        long best=1,cur=1,bi=0,ci=0; long hist[64]; memset(hist,0,sizeof(hist)); long runs=0; int perh[256]; ll perhcut[256]; memset(perh,0,sizeof(perh));
        for(long i=1;i<=nR;i++){ int same=0; if(i<nR){ ll x[7],y[7]; key(&R[i],x); key(&R[i-1],y); same=1; for(int t=0;t<nk[VAR];t++) if(x[t]!=y[t]) same=0; }
            if(i<nR&&same&&R[i].cut==R[i-1].cut) continue;
            if(i<nR&&same&&R[i].cut==R[i-1].cut+1){ cur++; continue; }
            hist[cur<63?cur:63]++; runs++; { ll hh=R[ci].k[0]; if(hh<256 && cur>perh[hh]){ perh[hh]=(int)cur; perhcut[hh]=R[ci].cut; } } if(cur>best){ best=cur; bi=ci; } cur=1; ci=i; }
        printf("   best bank by block count h (h:length@first cut):"); for(int hh=1;hh<256;hh+=2) if(perh[hh]>=2) printf(" %d:%d@%lld",hh,perh[hh],perhcut[hh]); printf("\n");
        ll kk[7]; key(&R[bi],kk);
        printf("%s\n   longest bank %ld: cut values %lld..%lld ; key",nm[VAR],best,R[bi].cut,R[bi].cut+best-1); for(int t=0;t<nk[VAR];t++) printf(" %lld",kk[t]); printf("\n   maximal runs by length:"); for(int i=2;i<64;i++) if(hist[i]) printf(" %d:%ld",i,hist[i]); printf("\n");
        ll last=R[bi].cut-1; long shown=0; for(long i=bi;i<nR&&shown<best;i++){ if(R[i].cut==last) continue; if(R[i].cut!=last+1) break; last=R[i].cut; shown++;
            printf("   seed cut=%lld order %d block:",R[i].cut,R[i].m); for(int t=0;t<R[i].m;t++) printf(" %d",blocks[R[i].blk][t]); printf(" X=%d\n",R[i].X); if(shown>=12&&best>14) { printf("   ...\n"); break; } }
    }
    return 0;
}
