/* 2-adic structure of H and hc on a switching cube.
   usage: two_adic_cube n seed [mode]
   mode 0 (default): random tournament from seed (xorshift64*); mode 1: transitive; mode 2: circulant from seed bits (n odd)
   Vertex n-1 is never switched (a switching set and its complement give the same tournament), so the cube has
   dimension n-1. For F in {H, hc}: c_U = sum_{S subset U} (-1)^{|U|-|S|} F(T^S)  (monomial coefficients in the
   0/1 switching variables). Prints min v_2(c_U) for each |U| and how many U attain it, and the F_2-degree of the
   binary digits 1..4 of H. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n;
static uint64_t rng;
static uint64_t nxt(void){ rng^=rng>>12; rng^=rng<<25; rng^=rng>>27; return rng*2685821657736338717ULL; }
static int64_t *dp;
static void counts(const int *out, int64_t *H, int64_t *hc){
    int full=(1<<n)-1; size_t sz=(size_t)(full+1)*n;
    memset(dp,0,sz*sizeof(int64_t));
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    for(int S=1;S<=full;S++){ int64_t *row=dp+(size_t)S*n; for(int v=0;v<n;v++){ int64_t c=row[v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; int w=__builtin_ctz(b); dp[(size_t)(S|b)*n+w]+=c; t^=b;} } }
    int64_t h=0; for(int v=0;v<n;v++) h+=dp[(size_t)full*n+v]; *H=h;
    memset(dp,0,sz*sizeof(int64_t)); dp[(size_t)1*n+0]=1;
    for(int S=1;S<=full;S+=2){ int64_t *row=dp+(size_t)S*n; for(int v=0;v<n;v++){ int64_t c=row[v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; int w=__builtin_ctz(b); dp[(size_t)(S|b)*n+w]+=c; t^=b;} } }
    int64_t c=0; for(int v=1;v<n;v++) if(out[v]&1) c+=dp[(size_t)full*n+v]; *hc=c;
}
static int v2(int64_t x){ if(x==0) return 99; int v=0; while(!(x&1)){x>>=1;v++;} return v; }
int main(int argc,char**argv){
    n=atoi(argv[1]); rng=strtoull(argv[2],0,10)*0x9E3779B97F4A7C15ULL+88172645463325252ULL; int mode=argc>3?atoi(argv[3]):0;
    for(int i=0;i<20;i++) nxt();
    int out0[20]; memset(out0,0,sizeof(out0));
    if(mode==2){ int half=(n-1)/2; int conn=0; for(int d=1;d<=half;d++){ if(nxt()>>63) conn|=1<<d; else conn|=1<<(n-d);} for(int i=0;i<n;i++) for(int d=1;d<n;d++) if((conn>>d)&1) out0[i]|=1<<((i+d)%n); }
    else for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ int b = mode==1 ? 1 : (int)(nxt()>>63); if(b) out0[i]|=1<<j; else out0[j]|=1<<i; }
    int k=n-1; size_t NC=(size_t)1<<k; int full=(1<<n)-1;
    dp=malloc(((size_t)1<<n)*n*sizeof(int64_t));
    int64_t *TH=malloc(NC*sizeof(int64_t)), *TC=malloc(NC*sizeof(int64_t)), *raw=malloc(NC*sizeof(int64_t));
    for(size_t X=0;X<NC;X++){
        int out[20];
        for(int i=0;i<n;i++){ int xi=(X>>i)&1; int o=out0[i]; int m = xi ? ((o&(int)X) | (~o & ~(int)X & full)) : ((o&~(int)X) | (~o & (int)X & full)); out[i]=m&~(1<<i)&full; }
        counts(out,&TH[X],&TC[X]); raw[X]=TH[X];
    }
    printf("n=%d seed=%s mode=%d  H(T)=%lld hc(T)=%lld\n",n,argv[2],mode,(long long)TH[0],(long long)TC[0]);
    for(int which=0;which<2;which++){
        int64_t *c = which? TC:TH;
        for(int i=0;i<k;i++) for(size_t m=0;m<NC;m++) if((m>>i)&1) c[m]-=c[m^((size_t)1<<i)];
        int minv[32]; long cntmin[32]; long viol=0; for(int j=0;j<32;j++){minv[j]=99;cntmin[j]=0;}
        for(size_t m=1;m<NC;m++){ int j=__builtin_popcountll(m); int v=v2(c[m]); int need=(j+1)/2-(which?1:0); if(v<need) viol++; if(v<minv[j]){minv[j]=v;cntmin[j]=1;} else if(v==minv[j]) cntmin[j]++; }
        printf("%s: violations of v_2(c_U) >= ceil(|U|/2)%s: %ld ;  min v_2 by |U| (count attaining):", which?"hc":"H ", which?"-1":"", viol);
        for(int j=1;j<=k;j++){ if(minv[j]==99) printf(" %d:inf",j); else printf(" %d:%d(%ld)",j,minv[j],cntmin[j]); }
        printf("\n");
    }
    unsigned char *bit=malloc(NC);
    printf("F_2-degree of binary digit k of H on the cube:");
    for(int kk=1;kk<=4;kk++){
        for(size_t m=0;m<NC;m++) bit[m]=(raw[m]>>kk)&1;
        for(int i=0;i<k;i++) for(size_t m=0;m<NC;m++) if((m>>i)&1) bit[m]^=bit[m^((size_t)1<<i)];
        int deg=-1; for(size_t m=0;m<NC;m++) if(bit[m]){ int j=__builtin_popcountll(m); if(j>deg) deg=j; }
        printf(" k=%d:%d",kk,deg);
    }
    printf("\n");
    return 0;
}
