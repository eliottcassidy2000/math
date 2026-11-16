/* Wave 43 (mac-mini-2026-11-16-odometer). Includes the exact two-apex count formula from the 7-21 repository (code/twoapex/apexlib.h, absolute path below). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int m;
#include "/Users/e/Documents/GitHub/7-21/code/twoapex/apexlib.h"
static ll H(unsigned head, int h){ m=h+7; unsigned u=head|(45u<<h); unsigned w=((1u<<h)-1)|((~68u&0x7Fu)<<h); return formula(u,w,0); }
static int cmp(const void*a,const void*b){ ll x=*(ll*)a,y=*(ll*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){ int h=atoi(argv[1]); unsigned N=1u<<h; ll *v=malloc(sizeof(ll)*N); unsigned *id=malloc(sizeof(unsigned)*N);
  ll (*pr)[2]=malloc(sizeof(ll)*2*N); for(unsigned i=0;i<N;i++){ pr[i][0]=H(i,h); pr[i][1]=i; }
  qsort(pr,N,sizeof(ll)*2,cmp);
  /* binary search helper */
  long q=0; int shown=0;
  for(unsigned i=0;i<N;i++){ ll b=pr[i][0]; if(i&&pr[i-1][0]==b) continue; int ok=1; ll ids[4]; ids[0]=pr[i][1];
    for(int k=1;k<4;k++){ ll t=b+6*k; long lo=0,hi=N-1,f=-1; while(lo<=hi){ long mid=(lo+hi)/2; if(pr[mid][0]<t) lo=mid+1; else if(pr[mid][0]>t) hi=mid-1; else { f=mid; break; } } if(f<0){ok=0;break;} ids[k]=pr[f][1]; }
    if(ok){ q++; if(shown<5){ shown++; printf("heads %lld %lld %lld %lld: counts %lld %lld %lld %lld -> with the window: every odd count %lld..%lld\n",ids[0],ids[1],ids[2],ids[3],b,b+6,b+12,b+18,b,b+24); } } }
  printf("head length %d: %ld base counts B with B+6, B+12, B+18 also base counts (2-level odometers)\n",h,q); return 0; }
