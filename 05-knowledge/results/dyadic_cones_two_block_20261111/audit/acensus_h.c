/* Auditor's own census of rooted tournaments (L, v), B = L - v, O = N+(v).
   Offset and response computed from end/start path tables of B:
     Q(e,f) = #ordered pairs (P1,P2) of disjoint nonempty paths covering B, P1 ends at e, P2 starts at f
     c0 = H(B) - sum_{e != f in O} Q(e,f) ;  M = sum_{y in O} sum_{P1 possibly empty} ...
   (THM-5126 (1')); validated separately by literal counting (averify.c).
   Modes:  acensus lab n          all 2^(n(n-1)/2) labelled tournaments of order n as bases (n <= 7)
           acensus cls nb         bases read from stdin (nb out-masks per line)
           acensus ext nb         all 2^nb one-vertex extensions of each listed base
   Output: lines "c r j" = least exponent j over roots with response c 2^j (c odd <= 255), |offset| <= 2047.  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define MX 10
#define R 2048
static int n; static int adj[MX][MX];
static ll endp[1<<MX][MX], stp[1<<MX][MX], hs[1<<MX];
static ll U[MX], Q[MX][MX], Hb;
static unsigned char best[128][2*R+1];
static long long nroots=0;
static unsigned char hval[1<<20];
static int curmask[MX]; static int wm[16][601][MX+2]; static int haveW[16][601];
static void tables(void){
    int F=(1<<n);
    for(int S=0;S<F;S++) for(int a=0;a<n;a++){ endp[S][a]=0; stp[S][a]=0; }
    for(int a=0;a<n;a++){ endp[1<<a][a]=1; stp[1<<a][a]=1; }
    /* increasing popcount order is automatic with increasing S for adding bits */
    for(int S=1;S<F;S++) for(int a=0;a<n;a++){
        if(!((S>>a)&1)) continue;
        if(S==(1<<a)) continue;
        int P=S&~(1<<a); ll e=0, s=0;
        for(int b=0;b<n;b++) if((P>>b)&1){ if(adj[b][a]) e+=endp[P][b]; if(adj[a][b]) s+=stp[P][b]; }
        endp[S][a]=e; stp[S][a]=s; }
    hs[0]=1; for(int S=1;S<F;S++){ ll h=0; for(int a=0;a<n;a++) h+=endp[S][a]; hs[S]=h; }
    Hb=hs[F-1];
    for(int a=0;a<n;a++){ U[a]=0; for(int b=0;b<n;b++) Q[a][b]=0; }
    for(int C=1;C<F;C++){ int S=(F-1)^C;
        for(int y=0;y<n;y++) if((C>>y)&1) U[y]+=hs[S]*stp[C][y];
        if(S) for(int e=0;e<n;e++) if((S>>e)&1) for(int f=0;f<n;f++) if((C>>f)&1) Q[e][f]+=endp[S][e]*stp[C][f]; }
}
static void roots(void){
    int F=1<<n;
    for(int O=1;O<F;O++){
        ll M=0, W=0;
        for(int y=0;y<n;y++) if((O>>y)&1){ M+=U[y]; for(int x=0;x<n;x++) if(x!=y && ((O>>x)&1)) W+=Q[x][y]; }
        ll c0=Hb-W;
        if(M>0){ ll hv=c0+2*M; if(hv>0 && hv<(1<<20)) hval[hv]=1; }
        if(M<=0) continue;
        int j=0; ll c=M; while(!(c&1)){ c>>=1; j++; }
        if(c>255 || c0>2047 || c0< -2047) continue;
        nroots++;
        unsigned char *p=&best[c>>1][c0+R]; if(j<*p){ *p=(unsigned char)j; if(c<=31 && c0>=-300 && c0<=300){ int*w=wm[c>>1][c0+300]; haveW[c>>1][c0+300]=1; w[0]=n; w[1]=O; for(int a=0;a<n;a++){ int mk=0; for(int b=0;b<n;b++) if(adj[a][b]) mk|=1<<b; w[2+a]=mk; } } }
    }
}
static void setfrommasks(const int *mask){ for(int a=0;a<n;a++) for(int b=0;b<n;b++) adj[a][b]=(a!=b)&&((mask[a]>>b)&1); }
int main(int argc,char**argv){
    memset(best,255,sizeof best);
    long nb_bases=0;
    if(!strcmp(argv[1],"lab")){ n=atoi(argv[2]); int E=n*(n-1)/2;
        for(long X=0;X<(1L<<E);X++){ int k=0; for(int a=0;a<n;a++) for(int b=a+1;b<n;b++){ int bit=(X>>k)&1; k++; adj[a][b]=bit; adj[b][a]=!bit; } for(int a=0;a<n;a++) adj[a][a]=0;
            tables(); roots(); nb_bases++; }
    } else { int nb=atoi(argv[2]); int ext=!strcmp(argv[1],"ext"); char line[512]; int base[MX];
        while(fgets(line,sizeof line,stdin)){ char *p=line; int ok=1; for(int i=0;i<nb;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} base[i]=(int)v; p=e; } if(!ok) continue;
            if(!ext){ n=nb; setfrommasks(base); tables(); roots(); nb_bases++; }
            else { n=nb+1; for(int X=0;X<(1<<nb);X++){ int mk[MX]; for(int i=0;i<nb;i++){ mk[i]=base[i]; if(!((X>>i)&1)) mk[i]|=1<<nb; } mk[nb]=X; setfrommasks(mk); tables(); roots(); nb_bases++; } }
        } }
    if(getenv("HVALS")){ for(int i=1;i<(1<<20);i+=2) if(hval[i]) printf("%d\n",i); return 0; }
    printf("# bases %ld order %d roots-in-window %lld\n",nb_bases,n,nroots);
    for(int ci=0;ci<128;ci++) for(int r=-2047;r<=2047;r++) if(best[ci][r+R]!=255) printf("%d %d %d\n",2*ci+1,r,best[ci][r+R]);
    for(int ci=0;ci<16;ci++) for(int r=-300;r<=300;r++) if(haveW[ci][r+300]){ int*w=wm[ci][r+300]; printf("W %d %d %d %d %d",2*ci+1,r,best[ci][r+R],w[0],w[1]); for(int a=0;a<w[0];a++) printf(" %d",w[2+a]); printf("\n"); }
    return 0;
}
