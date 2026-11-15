// For random tournaments of order n, the change of S = 3H + P2 under every single arc flip and every cyclic-triangle
// reversal: report how often |dS| is small (a switching method at unit scale needs small fixed steps).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static u64 s[2];
static u64 rnd(void){ u64 s1=s[0]; const u64 s0=s[1]; s[0]=s0; s1^=s1<<23; s[1]=s1^s0^(s1>>17)^(s0>>26); return s[1]+s0; }
static int n; static u64 *dp,*H;
static i64 S(unsigned*out){ int N=1<<n; memset(dp,0,sizeof(u64)*N*n);
  for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1; H[0]=1;
  for(int S_=1;S_<N;S_++){ u64 t=0; u64*row=dp+(size_t)S_*n; unsigned m=S_; while(m){int v=__builtin_ctz(m); m&=m-1; u64 x=row[v]; if(!x) continue; t+=x; unsigned nb=out[v]&~S_; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S_|(1<<u))*n+u]+=x;}} H[S_]=t; }
  u64 P2=0; for(int S_=1;S_<N-1;S_++) P2+=H[S_]*H[(N-1)^S_]; return 3*(i64)H[N-1]+(i64)P2; }
int main(int argc,char**argv){ n=atoi(argv[1]); int T=atoi(argv[2]); s[0]=12345; s[1]=67890;
  dp=malloc(sizeof(u64)*(1<<n)*n); H=malloc(sizeof(u64)*(1<<n));
  long tot_f=0, small_f[5]={0}, tot_t=0, small_t[5]={0}; long hf[26]={0}, ht[26]={0}; i64 bounds[5]={2,10,100,1000,10000}; double meanabs=0; long cnt=0;
  for(int it=0;it<T;it++){ unsigned out[32]; memset(out,0,sizeof out);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(rnd()&1) out[i]|=1u<<j; else out[j]|=1u<<i; }
    i64 s0=S(out);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ unsigned o2[32]; memcpy(o2,out,sizeof o2); o2[i]^=1u<<j; o2[j]^=1u<<i; // flip i-j
        i64 d=S(o2)-s0; if(d<0)d=-d; if(d<=25) hf[d]++; tot_f++; meanabs+=(double)d/s0; cnt++; for(int b=0;b<5;b++) if(d<=bounds[b]) small_f[b]++; }
    for(int a=0;a<n;a++) for(int b=a+1;b<n;b++) for(int c=b+1;c<n;c++){ int ab=(out[a]>>b)&1, bc=(out[b]>>c)&1, ca=(out[c]>>a)&1;
        if(!((ab&&bc&&ca)||(!ab&&!bc&&!ca))) continue; unsigned o2[32]; memcpy(o2,out,sizeof o2);
        o2[a]^=1u<<b; o2[b]^=1u<<a; o2[b]^=1u<<c; o2[c]^=1u<<b; o2[c]^=1u<<a; o2[a]^=1u<<c;   // reverse the triangle
        i64 d=S(o2)-s0; if(d<0)d=-d; if(d<=25) ht[d]++; tot_t++; for(int q=0;q<5;q++) if(d<=bounds[q]) small_t[q]++; } }
  printf("n=%d, %d random tournaments: mean |dS|/S over single flips %.4f\n",n,T,meanabs/cnt);
  printf("  single arc flips   (%ld): |dS|<=2: %ld  <=10: %ld  <=100: %ld  <=1000: %ld  <=10000: %ld\n",tot_f,small_f[0],small_f[1],small_f[2],small_f[3],small_f[4]);
  printf("  triangle reversals (%ld): |dS|<=2: %ld  <=10: %ld  <=100: %ld  <=1000: %ld  <=10000: %ld\n",tot_t,small_t[0],small_t[1],small_t[2],small_t[3],small_t[4]);
  printf("  |dS| histogram 0..24 flips: "); for(int q=0;q<25;q++) printf("%ld ",hf[q]); printf("\n  |dS| histogram 0..24 triangles: "); for(int q=0;q<25;q++) printf("%ld ",ht[q]); printf("\n");
  return 0; }
