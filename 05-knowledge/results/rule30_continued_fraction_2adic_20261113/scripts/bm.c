// Bit-packed Berlekamp-Massey over GF(2): linear complexity profile of a 0/1 byte file.
// Reports max |L(n) - n/2| and the histogram of jump sizes (= degrees of the continued-fraction partial quotients).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
static long N; static u64 *R; // R bit j = s_{N-1-j}
static inline int getR(long j){ return (R[j>>6]>>(j&63))&1; }
static inline u64 wordR(long j){ // 64 bits of R starting at bit j
  long w=j>>6; int sh=j&63; u64 lo=R[w]>>sh; if(sh==0) return lo; return lo | (R[w+1]<<(64-sh)); }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"rb"); fseek(f,0,SEEK_END); N=ftell(f); fseek(f,0,SEEK_SET); unsigned char*s=malloc(N); fread(s,1,N,f); fclose(f);
  if(argc>2){ long n2=atol(argv[2]); if(n2<N) N=n2; }
  long nw=N/64+3; R=calloc(nw+2,8); for(long j=0;j<N;j++) if(s[N-1-j]) R[j>>6]|=1ULL<<(j&63);
  u64 *C=calloc(nw+2,8), *B=calloc(nw+2,8), *Tm=calloc(nw+2,8);
  C[0]=1; B[0]=1; long L=0, m=1; long hist[200]={0}; long maxdev=0, maxdevn=0; long lastL=0;
  long Clen=1, Blen=1; // in bits (degree+1) upper bounds
  for(long n=0;n<N;n++){
    // discrepancy d = sum_{i=0}^{L} C_i s_{n-i};  s_{n-i} = R bit (N-1-n+i)
    long base=N-1-n; u64 acc=0;
    for(long i=0;i<=L;i+=64){ u64 cw=C[i>>6]; u64 rw=wordR(base+i); long rem=L+1-i; if(rem<64) cw&=((1ULL<<rem)-1); acc^=cw&rw; }
    int d=__builtin_popcountll(acc)&1;
    if(d){
      if(2*L<=n){ long newClen = (Blen+m > Clen? Blen+m : Clen); long words=(Clen+63)/64+1; memcpy(Tm,C,words*8);
        // C = C + x^m B
        long ws=m>>6; int bs=m&63; long bw=(Blen+63)/64+1;
        for(long w=0;w<bw;w++){ u64 v=B[w]; if(bs==0) C[w+ws]^=v; else { C[w+ws]^=v<<bs; C[w+ws+1]^=v>>(64-bs);} }
        long newL=n+1-L; long jump=newL-L; if(jump<200) hist[jump]++;
        long bwords=(Blen+63)/64+1; if(bwords>words) memset(B,0,bwords*8);
        memcpy(B,Tm,words*8); Blen=Clen; Clen=newClen;
        L=newL; m=1;
      } else {
        long ws=m>>6; int bs=m&63; long bw=(Blen+63)/64+1;
        for(long w=0;w<bw;w++){ u64 v=B[w]; if(bs==0) C[w+ws]^=v; else { C[w+ws]^=v<<bs; C[w+ws+1]^=v>>(64-bs);} }
        if(Blen+m>Clen) Clen=Blen+m; m++;
      }
    } else m++;
    long dev = 2*L-(n+1); if(dev<0) dev=-dev; if(dev>maxdev){maxdev=dev; maxdevn=n+1;}
  }
  printf("N=%ld final L=%ld (N/2=%ld); max |2L(n)-n| = %ld at n=%ld\n",N,L,N/2,maxdev,maxdevn);
  printf("jump-size (partial quotient degree) histogram:"); for(int k=0;k<200;k++) if(hist[k]) printf(" %d:%ld",k,hist[k]); printf("\n");
  return 0;
}
