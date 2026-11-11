/* Statistics of the quadratic function of H mod 4 on random switching classes.
   usage: rank_stats n trials seed
   For a random tournament T: lambda_a = (H(T^a)-H(T))/2 mod 2, M_ab = (H(T)+H(T^a)+H(T^b)+H(T^ab))/2 mod 2.
   By THM-5139, (H(T^X)-H(T))/2 = sum_{a in X} lambda_a + sum_{a<b in X} M_ab (mod 2) for every X.
   Prints: histogram of corank of M on F_2^n/<1> (dimension n-1; corank = n-1-rank), and of the class imbalance
   G = #{X: H(T^X)=H(T) mod 4} - #{X: H(T^X)=H(T)+2 mod 4} over the 2^(n-1) members, computed from the form
   (classes: G=0, or |G| = 2^(n-1-rank/2)), with the sign relative to H = 1 mod 4. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint64_t rng;
static uint64_t nxt(void){ rng^=rng>>12; rng^=rng<<25; rng^=rng>>27; return rng*2685821657736338717ULL; }
static uint64_t *dp;
static uint64_t Hcount(const int *out){
    int full=(1<<n)-1; memset(dp,0,(size_t)(full+1)*n*sizeof(uint64_t));
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    for(int S=1;S<=full;S++){ uint64_t *row=dp+(size_t)S*n; for(int v=0;v<n;v++){ uint64_t c=row[v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; int w=__builtin_ctz(b); dp[(size_t)(S|b)*n+w]+=c; t^=b;} } }
    uint64_t h=0; for(int v=0;v<n;v++) h+=dp[(size_t)full*n+v]; return h;
}
static void sw(const int *o0,int X,int *out){ int full=(1<<n)-1; for(int i=0;i<n;i++){ int o=o0[i]; int m=((X>>i)&1)? ((o&X)|(~o&~X&full)) : ((o&~X)|(~o&X&full)); out[i]=m&~(1<<i)&full; } }
int main(int argc,char**argv){
    n=atoi(argv[1]); long trials=atol(argv[2]); rng=strtoull(argv[3],0,10)*0x9E3779B97F4A7C15ULL+1234567ULL; for(int i=0;i<20;i++) nxt();
    dp=malloc(((size_t)1<<n)*n*sizeof(uint64_t));
    long corank[32]; memset(corank,0,sizeof(corank)); long zeroG=0,posG=0,negG=0; double sumG=0, sumabs=0; long h1=0;
    for(long t=0;t<trials;t++){
        int o0[20]; memset(o0,0,sizeof(o0)); int out[20];
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(nxt()>>63) o0[i]|=1<<j; else o0[j]|=1<<i; }
        uint64_t H0=Hcount(o0); uint64_t Ha[20]; int lam[20]; int M[20]; memset(M,0,sizeof(M));
        for(int a=0;a<n;a++){ sw(o0,1<<a,out); Ha[a]=Hcount(out); lam[a]=(int)(((Ha[a]-H0)>>1)&1); }
        for(int a=0;a<n;a++) for(int b=a+1;b<n;b++){ sw(o0,(1<<a)|(1<<b),out); uint64_t Hab=Hcount(out); int m=(int)(((H0+Ha[a]+Ha[b]+Hab)>>1)&1); if(m){M[a]|=1<<b; M[b]|=1<<a;} }
        /* rank of M over F_2 */
        int rows[20]; memcpy(rows,M,sizeof(int)*n); int rank=0;
        for(int col=0;col<n;col++){ int piv=-1; for(int r=rank;r<n;r++) if((rows[r]>>col)&1){piv=r;break;} if(piv<0) continue; int tmp=rows[piv]; rows[piv]=rows[rank]; rows[rank]=tmp; for(int r=0;r<n;r++) if(r!=rank && ((rows[r]>>col)&1)) rows[r]^=rows[rank]; rank++; }
        corank[n-1-rank]++;
        /* imbalance over X subset of {0..n-2} by the quadratic formula */
        long G=0; int k=n-1;
        for(int X=0;X<(1<<k);X++){ int q=0; for(int a=0;a<k;a++) if((X>>a)&1){ q^=lam[a]; q^=__builtin_parity(M[a]&X&((1<<a)-1)); } G+= q? -1:1; }
        int sgn = ((H0&3)==1)? 1:-1;   /* G counts members with the residue of T; convert to (#H=1) - (#H=3) */
        G*=sgn; if(G==0) zeroG++; else if(G>0) posG++; else negG++;
        sumG+=(double)G/(double)(1<<k); sumabs+= (G<0?-(double)G:(double)G)/(double)(1<<k); if((H0&3)==1) h1++;
    }
    printf("n=%d trials=%ld seed=%s\n",n,trials,argv[3]);
    printf("corank of M on F_2^n/<1> (dim %d):",n-1); for(int c=0;c<n;c++) if(corank[c]) printf(" %d:%ld",c,corank[c]); printf("\n");
    printf("class imbalance (#H=1 - #H=3 mod 4)/2^(n-1): zero %ld, positive %ld, negative %ld ; mean %.6f ; mean |.| %.6f ; fraction of sampled T with H=1 mod 4: %.4f\n",zeroG,posG,negG,sumG/trials,sumabs/trials,(double)h1/trials);
    return 0;
}
