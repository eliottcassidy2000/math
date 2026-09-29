/* 3x-1 sheet on positive integers, T(n)=n/2 (even), (3n-1)/2 (odd).
   For every n<=N record the first ray-hit: the first orbit element v with
   odd part in the odd cycle set {1; 5,7; 17,25,37,55,41,61,91}, coded as
   (cycle-odd-index, v_2(v)). Basin = cycle of that odd part.
   Memo: if the orbit drops below n, inherit the code of the smaller value. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static const uint64_t oddcyc[10]={1,5,7,17,25,37,55,41,61,91};
static const int cycof[10]={0,1,1,2,2,2,2,2,2,2};
static int oddindex(uint64_t o){ for(int i=0;i<10;i++) if(oddcyc[i]==o) return i; return -1; }
int main(int argc,char**argv){
  uint64_t N = argc>1 ? strtoull(argv[1],0,10) : (1ULL<<29);
  uint16_t *code = malloc((N+1)*sizeof(uint16_t));
  if(!code){fprintf(stderr,"alloc fail\n");return 1;}
  memset(code,0xff,(N+1)*sizeof(uint16_t));
  /* code = idx*64 + k  where v = oddcyc[idx]*2^k is the first ray element hit */
  uint64_t maxsteps=0; u128 maxv=0;
  for(uint64_t n=1;n<=N;n++){
    u128 v=n; uint64_t steps=0;
    for(;;){
      /* check ray membership of v */
      uint64_t k=0; u128 o=v; while(!(o&1)){o>>=1;k++;}
      if(o<=91){ int id=oddindex((uint64_t)o); if(id>=0){ code[n]=(uint16_t)(id*64+(k>63?63:k)); break; } }
      if(v<n){ code[n]=code[(uint64_t)v]; break; }
      if(v&1) v=(3*v-1)/2; else v>>=1;
      steps++; if(v>maxv) maxv=v;
    }
    if(steps>maxsteps) maxsteps=steps;
  }
  /* dyadic-range statistics */
  printf("N=%llu maxsteps=%llu maxv~2^%d\n",(unsigned long long)N,(unsigned long long)maxsteps,(int)(127-__builtin_clzll((uint64_t)(maxv>>64))) );
  int S=0; while((1ULL<<(S+1))<=N) S++;
  printf("range basin1 basin2 basin3 | ray entries (idx,k):count top 12\n");
  uint64_t tot[3]={0,0,0};
  for(int s=0;s<=S;s++){
    uint64_t lo=1ULL<<s, hi=(1ULL<<(s+1))-1; if(hi>N) hi=N;
    uint64_t b[3]={0,0,0};
    static uint64_t cnt[640]; memset(cnt,0,sizeof cnt);
    for(uint64_t n=lo;n<=hi;n++){ uint16_t c=code[n]; if(c==0xffff){continue;} int id=c/64; b[cycof[id]]++; cnt[c]++; }
    uint64_t w=hi-lo+1;
    printf("[2^%d,2^%d) %llu %llu %llu  dens %.6f %.6f %.6f\n",s,s+1,(unsigned long long)b[0],(unsigned long long)b[1],(unsigned long long)b[2],(double)b[0]/w,(double)b[1]/w,(double)b[2]/w);
    for(int i=0;i<3;i++) tot[i]+=b[i];
    if(s==S){ /* print entry spectrum for the top range and cumulative */
      printf("top-range ray-entry spectrum (odd,k,density):\n");
      for(int c=0;c<640;c++) if(cnt[c]>0 && (double)cnt[c]/w>0.002) printf("  odd=%llu k=%d dens=%.5f\n",(unsigned long long)oddcyc[c/64],c%64,(double)cnt[c]/w);
    }
  }
  printf("cumulative densities [1,N]: %.7f %.7f %.7f\n",(double)tot[0]/N,(double)tot[1]/N,(double)tot[2]/N);
  /* cumulative entry spectrum */
  static uint64_t cnt[640]; memset(cnt,0,sizeof cnt);
  for(uint64_t n=1;n<=N;n++){ uint16_t c=code[n]; if(c!=0xffff) cnt[c]++; }
  printf("cumulative ray-entry spectrum (odd,k,density) for dens>1e-4:\n");
  for(int c=0;c<640;c++) if(cnt[c]>0 && (double)cnt[c]/N>1e-4) printf("  odd=%llu k=%d dens=%.6f\n",(unsigned long long)oddcyc[c/64],c%64,(double)cnt[c]/N);
  free(code); return 0;
}
