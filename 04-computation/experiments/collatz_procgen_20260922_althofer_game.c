// Althofer 3n+-1 game: position odd n; mover picks 3n+1 or 3n-1, halves to odd; reaching 1 wins for the mover.
// Retrograde fixpoint on odd n <= CAP: W=1 (N: mover wins), W=2 (P: mover loses), 0 undetermined (children beyond CAP unknown).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static inline uint64_t op(uint64_t x){ while(!(x&1)) x>>=1; return x; }
int main(int argc,char**argv){
  uint64_t CAP=atoll(argv[1]); uint64_t M=(CAP+1)/2+1; // index i -> n=2i+1
  uint8_t *W=calloc(M,1);
  int changed=1, it=0;
  while(changed){ changed=0; it++;
    for(uint64_t i=1;i<M;i++){ uint64_t n=2*i+1; if(n>CAP) break; if(W[i]) continue;
      uint64_t a=op(3*n+1), b=op(3*n-1);
      if(a==1||b==1){W[i]=1;changed=1;continue;}
      int wa = (a<=CAP)? W[(a-1)/2] : 0, wb=(b<=CAP)? W[(b-1)/2] : 0;
      if(wa==2||wb==2){W[i]=1;changed=1;continue;}
      if(wa==1&&wb==1){W[i]=2;changed=1;continue;}
    }
  }
  uint64_t cN=0,cP=0,c0=0, firstund=0; for(uint64_t i=1;i<M;i++){ uint64_t n=2*i+1; if(n>CAP) break; if(W[i]==1)cN++; else if(W[i]==2)cP++; else {c0++; if(!firstund) firstund=n;} }
  printf("CAP=%llu iterations=%d: N=%llu P=%llu undetermined=%llu (first undetermined n=%llu)\n",(unsigned long long)CAP,it,(unsigned long long)cN,(unsigned long long)cP,(unsigned long long)c0,(unsigned long long)firstund);
  // P-positions list (first 40) and residues
  printf("first P-positions:"); int k=0; for(uint64_t i=1;i<M && k<40;i++){ if(W[i]==2){ printf(" %llu",(unsigned long long)(2*i+1)); k++; } } printf("\n");
  // P density by residue mod 8, 16
  for(int mod=8; mod<=32; mod*=2){ printf("P fraction by n mod %d:",mod); for(int r=1;r<mod;r+=2){ uint64_t t=0,p=0; for(uint64_t n=r;n<=CAP;n+=mod){ if(n<3) continue; uint64_t i=(n-1)/2; t++; if(W[i]==2)p++; } printf(" %d:%.3f",r,t?(double)p/t:0); } printf("\n"); }
  return 0;
}
