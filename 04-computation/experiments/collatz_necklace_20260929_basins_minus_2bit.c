/* 3x-1 sheet basins to N=2^32 with 2 bits per integer (basin 1,2,3; 0 = unset). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static uint8_t *A; 
static inline int getb(uint64_t n){ return (A[n>>2]>>((n&3)*2))&3; }
static inline void setb(uint64_t n,int v){ A[n>>2] |= (uint8_t)(v<<((n&3)*2)); }
static int basin_of_odd(uint64_t o){ if(o==1) return 1; if(o==5||o==7) return 2; if(o==17||o==25||o==37||o==55||o==41||o==61||o==91) return 3; return 0; }
int main(int argc,char**argv){
  uint64_t N = argc>1 ? strtoull(argv[1],0,10) : (1ULL<<32);
  A = calloc(N/4+2,1); if(!A){fprintf(stderr,"alloc fail\n");return 1;}
  for(uint64_t n=1;n<=N;n++){
    u128 v=n;
    for(;;){
      u128 o=v; while(!(o&1)) o>>=1;
      if(o<=91){ int b=basin_of_odd((uint64_t)o); if(b){ setb(n,b); break; } }
      if(v<n){ setb(n,getb((uint64_t)v)); break; }
      if(v&1) v=(3*v-1)/2; else v>>=1;
    }
  }
  int S=0; while((1ULL<<(S+1))<=N) S++;
  uint64_t tot[4]={0,0,0,0};
  printf("N=%llu\nrange dens1 dens2 dens3\n",(unsigned long long)N);
  for(int s=20;s<=S;s++){
    uint64_t lo=1ULL<<s, hi=(1ULL<<(s+1))-1; if(hi>N) hi=N; uint64_t b[4]={0,0,0,0};
    for(uint64_t n=lo;n<=hi;n++) b[getb(n)]++;
    uint64_t w=hi-lo+1;
    printf("[2^%d,2^%d) %.7f %.7f %.7f unset=%llu\n",s,s+1,(double)b[1]/w,(double)b[2]/w,(double)b[3]/w,(unsigned long long)b[0]);
  }
  for(uint64_t n=1;n<=N;n++) tot[getb(n)]++;
  printf("cumulative [1,N]: %.8f %.8f %.8f unset=%llu\n",(double)tot[1]/N,(double)tot[2]/N,(double)tot[3]/N,(unsigned long long)tot[0]);
  return 0;
}
