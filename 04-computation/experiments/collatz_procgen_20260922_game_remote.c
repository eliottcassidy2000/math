// Steinhaus remoteness (length of optimal play: winner wins as fast as possible, loser loses as slowly
// as possible) for Althofer's 3n+-1 game, on all odd non-multiples of 3 up to CAP, by LEVEL SWEEPS:
//   level r: an unresolved x becomes N with remoteness r if some child is P with remoteness r-1
//            (child 1 = terminal, remoteness 0); it becomes P with remoteness r if both children are
//            N and the larger remoteness is r-1.   (Children above CAP are unknown.)
// Capped remoteness is an upper bound for the true remoteness (exact when all optimal lines stay
// <= CAP); stability between two caps is reported by the driver.
// Then statistics on the fully resolved prefix: remoteness distribution vs log2 n, record holders,
// and which move optimal play uses (winner: fastest winning child; loser: slowest child).
// usage: game_remote CAP [CAP2-for-comparison-file] ; writes rem_<CAP>.bin (uint16 per index)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

static uint64_t CAP, MAXIDX;
static uint16_t *R;   // 0 unresolved, else remoteness+1
static inline uint64_t idx_of(uint64_t n){ return n/3; }
static inline uint64_t n_of(uint64_t i){ return 3*i + 1 + (i&1); }
static inline uint64_t oddpart(uint64_t x){ return x >> __builtin_ctzll(x); }
static inline uint64_t dchild(uint64_t n){ return oddpart((n&3)==1 ? 3*n+1 : 3*n-1); }
static inline uint64_t uchild(uint64_t n){ return ((n&3)==1 ? 3*n-1 : 3*n+1) >> 1; }
static inline int rem(uint64_t y){ if(y==1) return 0; if(y>CAP) return -1; int r=R[idx_of(y)]; return r? r-1 : -1; }
static int rem_start(uint64_t n){   // any odd start
  if(n==1) return 1;
  if(n%3) return rem(n);
  int ra=rem(dchild(n)), rb=rem(uchild(n)); int best=-1;
  if(ra>=0 && ra%2==0) best=ra+1; if(rb>=0 && rb%2==0 && (best<0 || rb+1<best)) best=rb+1;
  if(best>=0) return best; if(ra>=0 && rb>=0) return 1+(ra>rb?ra:rb); return -1;
}
int main(int argc,char**argv){
  CAP=strtoull(argv[1],0,10); if(CAP%2==0) CAP--;
  const char* cmpf = argc>2? argv[2]:NULL;
  MAXIDX=idx_of(CAP)+1; R=calloc(MAXIDX,2);
  time_t t0=time(0);
  R[0]=1; // terminal 1: remoteness 0
  uint64_t unresolved=MAXIDX-1; int r;
  for(r=1;r<65000;r++){
    uint64_t newc=0;
    for(uint64_t i=1;i<MAXIDX;i++){
      if(R[i]) continue; uint64_t x=n_of(i); if(x>CAP) break;
      int ra=rem(dchild(x)), rb=rem(uchild(x));
      if((ra==r-1 && ra%2==0) || (rb==r-1 && rb%2==0)){ R[i]=(uint16_t)(r+1); newc++; continue; }
      if(ra>=0 && rb>=0 && ra%2==1 && rb%2==1 && (ra>rb?ra:rb)==r-1){ R[i]=(uint16_t)(r+1); newc++; }
    }
    // note: nodes set during this sweep have remoteness r; they are not used as r-1 in the same sweep
    // because the tests compare with r-1 only.
    unresolved-=newc;
    if(!newc && r>2) break;
  }
  printf("# game_remote CAP=%llu: levels=%d unresolved=%llu (%lds)\n",(unsigned long long)CAP,r,(unsigned long long)unresolved,(long)(time(0)-t0));
  uint64_t B=3; while(B<=CAP && rem_start(B)>=0) B+=2;
  printf("fully resolved prefix: every odd start < %llu\n",(unsigned long long)B);
  // distribution by octave
  printf("remoteness by octave (odd n in [2^k,2^{k+1}), n<B): mean, max (argmax)\n");
  for(int k=1;k<64;k++){ uint64_t lo=1ULL<<k, hi=lo<<1; if(lo>=B) break; double s=0; uint64_t c=0; int mx=-1; uint64_t am=0;
    for(uint64_t n=lo+1;n<hi && n<B;n+=2){ int q=rem_start(n); s+=q; c++; if(q>mx){mx=q; am=n;} }
    printf(" k=%2d: mean %.2f  max %d at n=%llu%s\n",k,s/c,mx,(unsigned long long)am, hi>B? " (partial octave)":""); }
  // record holders (running max)
  printf("remoteness records (odd n < B):"); { int mx=-1; for(uint64_t n=1;n<B;n+=2){ int q=rem_start(n); if(q>mx){ mx=q; printf(" %llu:%d",(unsigned long long)n,q);} } } printf("\n");
  // optimal move usage on odd n<B with both children resolved
  { uint64_t Nn=0,Pn=0, winUp=0, winUpForced=0, winTie=0, loseUp=0, loseTie=0;
    for(uint64_t n=3;n<B;n+=2){ uint64_t a=dchild(n), b=uchild(n); int ra=rem(a), rb=rem(b); int q=rem_start(n); if(q<0) continue;
      if(q%2==1){ Nn++; int pa=(ra>=0&&ra%2==0), pb=(rb>=0&&rb%2==0);
        if(pb && (!pa || rb<ra)) { winUp++; if(!pa) winUpForced++; } else if(pa&&pb&&ra==rb) winTie++; }
      else { Pn++; if(ra<0||rb<0) continue; if(rb>ra) loseUp++; else if(ra==rb) loseTie++; } }
    printf("optimal play (odd 3<=n<B): N-positions %llu: fastest win is ascending in %.4f (forced: only ascent wins %.4f; ties %.4f)\n",(unsigned long long)Nn,(double)winUp/Nn,(double)winUpForced/Nn,(double)winTie/Nn);
    printf("                           P-positions %llu: slowest loss is ascending in %.4f (ties %.4f)\n",(unsigned long long)Pn,(double)loseUp/Pn,(double)loseTie/Pn); }
  // principal variation statistics for some starts
  { uint64_t starts[]={7,13,15,23,27,41417,817837,834437,0}; for(int s=0;starts[s];s++){ uint64_t n=starts[s]; if(n>=B) continue; int q=rem_start(n); printf("PV from %llu (remoteness %d):",(unsigned long long)n,q);
      int ups[2]={0,0}, ply=0; uint64_t x=n; while(x!=1){ uint64_t a=dchild(x), b=uchild(x); int ra=rem(a), rb=rem(b); int rx=rem_start(x); uint64_t nx; int up;
        if(rx%2==1){ int pa=(a==1)||(ra>=0&&ra%2==0), pb=(rb>=0&&rb%2==0); if(a==1){nx=a;up=0;} else if(pb && (!pa || rb<ra)){nx=b;up=1;} else {nx=a;up=0;} }
        else { if(rb>ra){nx=b;up=1;} else {nx=a;up=0;} }
        if(ply<24) printf(" %s%llu",up?"^":"v",(unsigned long long)nx); ups[ply&1]+=up; ply++; x=nx; if(ply>1000) break; }
      printf("%s  [plies %d; ascents by first player %d, by second %d]\n", ply>24?" ...":"", ply, ups[0], ups[1]); } }
  if(cmpf){ // compare with a smaller-cap run on its resolved prefix
    FILE*F=fopen(cmpf,"rb"); uint64_t c2; fread(&c2,8,1,F); uint64_t m2=idx_of(c2)+1; uint16_t *R2=malloc(m2*2); fread(R2,2,m2,F); fclose(F);
    uint64_t same=0,diff=0,res2=0; int maxd=0; for(uint64_t i=1;i<m2;i++){ if(!R2[i]) continue; res2++; if(R2[i]==R[i]) same++; else { diff++; int d=(int)R2[i]-(int)R[i]; if(d>maxd) maxd=d; if(R2[i]<R[i]) printf("  !! smaller cap gives SMALLER remoteness at n=%llu\n",(unsigned long long)n_of(i)); } }
    printf("comparison with cap %llu: of %llu positions resolved there, remoteness equal for %llu, different for %llu (max excess %d)\n",(unsigned long long)c2,(unsigned long long)res2,(unsigned long long)same,(unsigned long long)diff,maxd); }
  char fn[256]; snprintf(fn,sizeof fn,"rem_%llu.bin",(unsigned long long)CAP); FILE*F=fopen(fn,"wb"); fwrite(&CAP,8,1,F); fwrite(R,2,MAXIDX,F); fclose(F);
  return 0;
}
