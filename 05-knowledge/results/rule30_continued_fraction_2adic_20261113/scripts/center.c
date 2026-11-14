// Rule 30 single-seed center column c_0..c_{N-1}, bit-packed evolution restricted to the shrinking dependency cone.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
int main(int argc,char**argv){
  long N=atol(argv[1]);
  long W=2*N+2;                     // cells -N..N mapped to bit (i+N+1); word array
  long nw=W/64+2;
  u64 *a=calloc(nw,8), *b=calloc(nw,8);
  long off=N+1; a[off>>6]|=1ULL<<(off&63);
  unsigned char *c=malloc(N);
  for(long t=0;t<N;t++){
    c[t]=(a[off>>6]>>(off&63))&1;
    // needed cells at time t for center at time N-1: |i| <= N-1-t ; also within light cone |i|<=t. update range words
    long rad = t+1; if(N-1-t < rad) rad = N-1-t+1; 
    long lo=(off-rad-1)>>6, hi=(off+rad+1)>>6; if(lo<0) lo=0; if(hi>=nw-1) hi=nw-2;
    for(long w=lo; w<=hi; w++){
      u64 x=a[w];
      u64 left = (x<<1) | (w>0? a[w-1]>>63 : 0);      // bit j gets cell j-1 (left neighbour)
      u64 right= (x>>1) | (a[w+1]<<63);               // bit j gets cell j+1
      b[w]= left ^ (x | right);
    }
    for(long w=lo; w<=hi; w++) a[w]=b[w];
  }
  FILE*f=fopen(argv[2],"wb"); fwrite(c,1,N,f); fclose(f);
  long ones=0; for(long t=0;t<N;t++) ones+=c[t];
  printf("N=%ld ones=%ld first 64:",N,ones); for(int t=0;t<64&&t<N;t++) printf("%d",c[t]); printf("\n");
  return 0;
}
