/* 3x+1 sheet: for every n<=N, the index j of the first power of two 4^j hit by the orbit
   (first trunk element; always an even power for n not itself a power of two). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
int main(int argc,char**argv){
  uint64_t N = argc>1 ? strtoull(argv[1],0,10) : (1ULL<<30);
  uint8_t *j = malloc(N+1); if(!j){fprintf(stderr,"alloc fail\n");return 1;}
  memset(j,0xff,N+1);
  for(uint64_t k=0;(1ULL<<k)<=N;k++) j[1ULL<<k] = 0xfe; /* powers of two themselves: excluded (density zero) */
  uint64_t bad=0;
  for(uint64_t n=1;n<=N;n++){
    if(j[n]!=0xff) continue;
    u128 v=n;
    for(;;){
      if((v&(v-1))==0){ /* power of two: first trunk hit */
        int k=0; u128 t=v; while(t>1){t>>=1;k++;}
        if(k%2) j[n]=(uint8_t)((k+1)/2); else { bad++; j[n]=0xfd; } /* T-form: first trunk hit is 2^(2j-1), entered from (4^j-1)/3 */
        break; }
      if(v<n){ j[n]=j[(uint64_t)v]; if(j[n]==0xfe||j[n]==0xfd) bad++; break; }
      if(v&1) v=(3*v+1)/2; else v>>=1;
    }
  }
  printf("N=%llu bad(odd-power first hits)=%llu\n",(unsigned long long)N,(unsigned long long)bad);
  int S=0; while((1ULL<<(S+1))<=N) S++;
  static uint64_t cnt[256];
  printf("per dyadic range: density of first trunk hit at 4^j, j=2..12\n");
  for(int s=8;s<=S;s++){
    uint64_t lo=1ULL<<s, hi=(1ULL<<(s+1))-1; if(hi>N) hi=N; uint64_t w=hi-lo+1;
    memset(cnt,0,sizeof cnt);
    for(uint64_t n=lo;n<=hi;n++) cnt[j[n]]++;
    printf("[2^%d,2^%d):",s,s+1);
    for(int q=2;q<=12;q++) printf(" %.5f",(double)cnt[q]/w);
    printf("  (tail>12: %.6f)\n",(double)( w - cnt[0]-cnt[1]-cnt[2]-cnt[3]-cnt[4]-cnt[5]-cnt[6]-cnt[7]-cnt[8]-cnt[9]-cnt[10]-cnt[11]-cnt[12]-cnt[0xfe]-cnt[0xfd]-cnt[0xff])/w);
  }
  memset(cnt,0,sizeof cnt);
  for(uint64_t n=1;n<=N;n++) cnt[j[n]]++;
  printf("cumulative [1,N] e_j = density of first trunk hit at 4^j:\n");
  double cum=0;
  for(int q=1;q<=40;q++){ if(cnt[q]==0) continue; cum+= (double)cnt[q]/N; printf("  j=%d 4^j=%llu e_j=%.7f cum=%.7f  D(4^j)=1-cum_{<j}=%.7f\n",q,(unsigned long long)1<<(2*q),(double)cnt[q]/N,cum,1.0-(cum-(double)cnt[q]/N)); }
  free(j); return 0;
}
