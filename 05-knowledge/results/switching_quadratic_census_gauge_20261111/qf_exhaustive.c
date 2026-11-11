/* Exhaustive check of the endpoint-graph theorem (THM-5139) on ALL labelled tournaments of order n <= 7.
   usage: qf_exhaustive n
   For every tournament T (bit p(i,j)=1 iff i->j, i<j):
     e[T]  = parity vector of s(v)+t(v)  (Hamiltonian paths starting / ending at v)
     H mod 4, hc mod 2 (Hamiltonian cycles)
   Checks, for every T:
     (a) Lambda(T) := (e[T] xor e[T^w])_w is symmetric, zero diagonal, even row sums
     (b) Lambda(T^w) = Lambda(T) for every w   (class invariance => e affine on the class)
     (c) hc[T] + hc[T^a] = 1 + e_T(a)                              (mod 2)
     (d) hc[T] + hc[T^a] + hc[T^b] + hc[T^ab] = Lambda_ab          (mod 2)
     (e) sum over the 8 switchings at subsets of {a,b,c} of hc = 0 (mod 2), of H = 0 (mod 4)
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, npairs, pidx[8][8];
static uint32_t vmask[8];
static void compute(uint32_t T, uint8_t *e, uint8_t *h4, uint8_t *hcp){
    int out[8]; int full=(1<<n)-1;
    for(int i=0;i<n;i++) out[i]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if((T>>pidx[i][j])&1) out[i]|=1<<j; else out[j]|=1<<i; }
    static uint32_t dp[128][7];
    uint32_t tend[8], sstart[8];
    /* ends */
    memset(dp,0,sizeof(dp));
    for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ uint32_t c=dp[S][v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; int w=__builtin_ctz(b); dp[S|b][w]+=c; t^=b;} }
    for(int v=0;v<n;v++) tend[v]=dp[full][v];
    /* starts = ends of the converse */
    memset(dp,0,sizeof(dp));
    for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ uint32_t c=dp[S][v]; if(!c) continue; int t=(~out[v])&full&~(1<<v)&~S; while(t){int b=t&-t; int w=__builtin_ctz(b); dp[S|b][w]+=c; t^=b;} }
    for(int v=0;v<n;v++) sstart[v]=dp[full][v];
    uint32_t H=0; uint8_t ev=0;
    for(int v=0;v<n;v++){ H+=tend[v]; if((tend[v]+sstart[v])&1) ev|=1<<v; }
    /* Hamiltonian cycles through vertex 0 */
    memset(dp,0,sizeof(dp)); dp[1][0]=1;
    for(int S=1;S<=full;S+=2) for(int v=0;v<n;v++){ uint32_t c=dp[S][v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; int w=__builtin_ctz(b); dp[S|b][w]+=c; t^=b;} }
    uint32_t hc=0; for(int v=1;v<n;v++) if(out[v]&1) hc+=dp[full][v];
    *e=ev; *h4=H&3; *hcp=hc&1;
}
int main(int argc,char**argv){
    n=atoi(argv[1]); if(n<3||n>7){fprintf(stderr,"3<=n<=7\n");return 1;}
    npairs=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ pidx[i][j]=pidx[j][i]=npairs++; }
    for(int v=0;v<n;v++){ vmask[v]=0; for(int u=0;u<n;u++) if(u!=v) vmask[v]|=1u<<pidx[u][v]; }
    uint32_t NT=1u<<npairs;
    uint8_t *E=malloc(NT), *H4=malloc(NT), *HC=malloc(NT);
    for(uint32_t T=0;T<NT;T++) compute(T,&E[T],&H4[T],&HC[T]);
    unsigned long long bad_a=0,bad_b=0,bad_c=0,bad_d=0,bad_e_hc=0,bad_e_H=0, cnt=0;
    unsigned long long h1=0,h3=0;
    for(uint32_t T=0;T<NT;T++){
        cnt++; if(H4[T]==1) h1++; else h3++;
        uint8_t L[8];
        for(int w=0;w<n;w++) L[w]=E[T]^E[T^vmask[w]];   /* column w */
        for(int w=0;w<n;w++){
            if((L[w]>>w)&1) bad_a++;
            for(int u=0;u<n;u++) if(((L[w]>>u)&1)!=((L[u]>>w)&1)) bad_a++;
        }
        uint8_t rowsum=0; for(int w=0;w<n;w++) rowsum^=L[w]; if(rowsum) bad_a++;
        for(int w=0;w<n;w++){ uint32_t T2=T^vmask[w]; for(int x=0;x<n;x++){ if((E[T2]^E[T2^vmask[x]])!=L[x]) bad_b++; } }
        for(int a=0;a<n;a++) if((HC[T]^HC[T^vmask[a]])!=(1^((E[T]>>a)&1))) bad_c++;
        for(int a=0;a<n;a++) for(int b=a+1;b<n;b++){
            uint8_t v=HC[T]^HC[T^vmask[a]]^HC[T^vmask[b]]^HC[T^vmask[a]^vmask[b]];
            if(v!=((L[b]>>a)&1)) bad_d++;
            for(int c=b+1;c<n;c++){
                uint8_t x=0; unsigned s=0;
                for(int m=0;m<8;m++){ uint32_t T2=T; if(m&1)T2^=vmask[a]; if(m&2)T2^=vmask[b]; if(m&4)T2^=vmask[c]; x^=HC[T2]; s+=H4[T2]; }
                if(x) bad_e_hc++; if(s&3) bad_e_H++;
            }
        }
    }
    /* histogram, over switching classes, of the number of members with H = 1 (mod 4) */
    { uint16_t *c1=calloc(NT,sizeof(uint16_t)); unsigned long long hist[130]; memset(hist,0,sizeof(hist));
      for(uint32_t T=0;T<NT;T++){ uint32_t R=T; for(int v=0;v<n-1;v++) if(!((T>>pidx[v][n-1])&1)) {;} else R^=vmask[v]; /* v -> n-1: switch v */
          if(H4[T]==1) c1[R]++; }
      unsigned long long ncl=0; for(uint32_t T=0;T<NT;T++){ int src=1; for(int v=0;v<n-1;v++) if((T>>pidx[v][n-1])&1) src=0; if(src){ ncl++; hist[c1[T]]++; } }
      printf("switching classes: %llu, each with %d members; number of members with H = 1 (mod 4) : number of classes ->",ncl,1<<(n-1));
      for(int i=0;i<130;i++) if(hist[i]) printf(" %d:%llu",i,hist[i]); printf("\n"); free(c1); }
    printf("n=%d tournaments=%llu  H=1 mod 4: %llu  H=3 mod 4: %llu\n",n,cnt,h1,h3);
    printf("violations: (a) Lambda symmetric/zero diagonal/even rows: %llu  (b) class invariance: %llu  (c) D_a hc = 1+e(a): %llu  (d) D_a D_b hc = Lambda_ab: %llu  (e) triple sums: hc %llu, H mod 4 %llu\n",bad_a,bad_b,bad_c,bad_d,bad_e_hc,bad_e_H);
    return 0;
}
