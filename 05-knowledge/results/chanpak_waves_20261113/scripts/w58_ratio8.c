/* order-8 tournaments (extensions of order-7 reps): counts h with a vertex t of B(T,t) = 3h; report h mod 3 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 st[256][8], Hs[256]; static int nn, outm[8];
static void run(void){ int N=1<<nn; memset(st,0,sizeof st); for(int v=0;v<nn;v++) st[1<<v][v]=1;
  /* st[S][v] = paths covering S starting at v: v -> w ... */
  for(int S=1;S<N;S++){ if(__builtin_popcount(S)<2) continue; for(int v=0;v<nn;v++) if(S>>v&1){ int R=S&~(1<<v); u64 t=0; int o=outm[v]&R; while(o){ int w=__builtin_ctz(o); o&=o-1; t+=st[R][w]; } st[S][v]=t; } }
  Hs[0]=1; for(int S=1;S<N;S++){ u64 a=0; for(int v=0;v<nn;v++) a+=st[S][v]; Hs[S]=a; } }
int main(void){ FILE*fi=fopen("w57_reps.txt","r"); char line[512]; static unsigned char has[1<<16]; int reps7[456][7],k7=0;
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); int o[8]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10);
    nn=n; memcpy(outm,o,sizeof(int)*n); run(); int full=(1<<n)-1; u64 h=Hs[full];
    for(int v=0;v<n;v++){ u64 B=0; for(int S=1;S<=full;S++) if((S>>v&1)&&st[S][v]) B+=st[S][v]*Hs[full&~S]; if(B==3*h && h<(1<<16)) has[h]=1; }
    if(n==7) memcpy(reps7[k7++],o,sizeof(int)*7); }
  for(int r=0;r<k7;r++) for(int S=0;S<128;S++){ nn=8; for(int i=0;i<7;i++){ outm[i]=reps7[r][i]; if(!(S>>i&1)) outm[i]|=1<<7; } outm[7]=S; run();
      u64 h=Hs[255]; for(int v=0;v<8;v++){ u64 B=0; for(int T=1;T<=255;T++) if((T>>v&1)&&st[T][v]) B+=st[T][v]*Hs[255&~T]; if(B==3*h && h<(1<<16)) has[h]=1; } }
  int c=0,c3=0; printf("counts h with a ratio-3 root (B = 3h), orders <= 8:"); for(int h=1;h<(1<<16);h++) if(has[h]){ printf(" %d",h); c++; if(h%3==0) c3++; } printf("\n%d counts, %d divisible by 3\n",c,c3); return 0; }
