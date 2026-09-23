// Undirected shortcut-Collatz game: from m, moves T(x) (forward) and predecessors 2x, (2x-1)/3 [x=2 mod 3] (backward).
// BFS for the minimal number of moves to reach any value < m (values capped at CAP*m). Compare with forward stopping time.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define HS (1<<22)
static uint64_t keys[HS]; static int stamp[HS]; static int cur=0;
static int seen(uint64_t x){ uint64_t h=(x*11400714819323198485ull)>>42; while(1){ if(stamp[h]!=cur){stamp[h]=cur;keys[h]=x;return 0;} if(keys[h]==x) return 1; h=(h+1)&(HS-1);} }
int main(int argc,char**argv){
  uint64_t N=atoll(argv[1]); uint64_t CAP=atoll(argv[2]); int DMAX=atoi(argv[3]);
  static uint64_t q1[1<<21], q2[1<<21];
  int hist[64]={0}; int histT[256]={0}; int worstU=0; uint64_t worstm=0; long long sumU=0,sumT=0,cnt=0; int unresolved=0;
  for(uint64_t m=3;m<=N;m++){
    // forward stopping time
    uint64_t x=m; int st=0; while(x>=m){ x=(x&1)?(3*x+1)/2:x/2; st++; if(st>250) break; }
    cur++; int n1=0; q1[n1++]=m; seen(m); int depth=0, found=0;
    while(n1>0 && depth<DMAX && !found){
      int n2=0; depth++;
      for(int i=0;i<n1 && !found;i++){
        uint64_t y=q1[i]; uint64_t nb[3]; int nn=0;
        nb[nn++]=(y&1)?(3*y+1)/2:y/2;
        nb[nn++]=2*y;
        if(y%3==2) nb[nn++]=(2*y-1)/3;
        for(int t=0;t<nn;t++){ uint64_t z=nb[t]; if(z<m && z>=1){found=1;break;} if(z>CAP*m) continue; if(seen(z)) continue; if(n2<(1<<21)) q2[n2++]=z; }
      }
      memcpy(q1,q2,sizeof(uint64_t)*n2); n1=n2;
    }
    if(!found){unresolved++; continue;}
    hist[depth<63?depth:63]++; if(depth>worstU){worstU=depth;worstm=m;}
    sumU+=depth; sumT+=st; cnt++;
  }
  printf("N=%llu CAP=%llu: undirected min-descent depth: mean %.3f, max %d (at m=%llu), unresolved %d;  forward stopping time mean %.3f\n",
    (unsigned long long)N,(unsigned long long)CAP,(double)sumU/cnt,worstU,(unsigned long long)worstm,unresolved,(double)sumT/cnt);
  printf("depth histogram:"); for(int d=1;d<40;d++) if(hist[d]) printf(" %d:%d",d,hist[d]); printf("\n");
  return 0;
}
