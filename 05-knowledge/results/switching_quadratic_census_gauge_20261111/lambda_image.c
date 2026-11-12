/* Image and fibres of the endpoint graph: switching classes of labelled tournaments on n <= 7 points -> even graphs.
   usage: lambda_image n
   For each class (representative: the member in which vertex n-1 is a source) computes Lambda (THM-5139) and
   prints the number of distinct labelled graphs obtained, the distribution of fibre sizes, and the number of
   classes by sorted degree sequence of Lambda. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, npairs, pidx[8][8]; static uint32_t vmask[8];
static uint8_t evec(uint32_t T){
    int out[8]; int full=(1<<n)-1; for(int i=0;i<n;i++) out[i]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if((T>>pidx[i][j])&1) out[i]|=1<<j; else out[j]|=1<<i; }
    static uint32_t dp[128][7]; uint32_t tend[8], sst[8];
    memset(dp,0,sizeof(dp)); for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ uint32_t c=dp[S][v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b;} }
    for(int v=0;v<n;v++) tend[v]=dp[full][v];
    memset(dp,0,sizeof(dp)); for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ uint32_t c=dp[S][v]; if(!c) continue; int t=(~out[v])&full&~(1<<v)&~S; while(t){int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b;} }
    for(int v=0;v<n;v++) sst[v]=dp[full][v];
    uint8_t e=0; for(int v=0;v<n;v++) if((tend[v]+sst[v])&1) e|=1<<v; return e;
}
static int cmp64(const void*a,const void*b){ uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b; return x<y?-1:(x>y); }
int main(int argc,char**argv){
    n=atoi(argv[1]); npairs=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) pidx[i][j]=pidx[j][i]=npairs++;
    for(int v=0;v<n;v++){ vmask[v]=0; for(int u=0;u<n;u++) if(u!=v) vmask[v]|=1u<<pidx[u][v]; }
    uint32_t NT=1u<<npairs; size_t ncl=0; uint64_t *keys=malloc(sizeof(uint64_t)*(NT>>(n-1)));
    long degcount[4096]; memset(degcount,0,sizeof(degcount)); static char degname[4096][32]; int ndeg=0;
    for(uint32_t T=0;T<NT;T++){
        int src=1; for(int v=0;v<n-1;v++) if((T>>pidx[v][n-1])&1) src=0; if(!src) continue;
        uint8_t e0=evec(T); uint64_t key=0; int deg[8];
        for(int w=0;w<n;w++){ uint8_t col=e0^evec(T^vmask[w]); key|=((uint64_t)col)<<(8*w); deg[w]=__builtin_popcount(col); }
        keys[ncl++]=key;
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(deg[j]<deg[i]){int t=deg[i];deg[i]=deg[j];deg[j]=t;}
        char nm[32]; for(int i=0;i<n;i++) nm[i]='0'+deg[i]; nm[n]=0;
        int k; for(k=0;k<ndeg;k++) if(!strcmp(degname[k],nm)) break; if(k==ndeg){ strcpy(degname[ndeg++],nm); } degcount[k]++;
    }
    qsort(keys,ncl,sizeof(uint64_t),cmp64);
    size_t distinct=0; long fib[40000]; memset(fib,0,sizeof(fib)); size_t i=0; long maxf=0;
    while(i<ncl){ size_t j=i; while(j<ncl&&keys[j]==keys[i]) j++; long f=(long)(j-i); if(f<40000) fib[f]++; if(f>maxf) maxf=f; distinct++; i=j; }
    printf("n=%d classes=%zu distinct endpoint graphs=%zu ; fibre size : number of graphs ->",n,ncl,distinct);
    for(long f=1;f<=maxf&&f<40000;f++) if(fib[f]) printf(" %ld:%ld",f,fib[f]); printf("\n");
    printf("   classes by degree sequence of Lambda:"); for(int k=0;k<ndeg;k++) printf(" %s:%ld",degname[k],degcount[k]); printf("\n");
    return 0;
}
