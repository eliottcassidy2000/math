// Distinct counts realized by the proved free-address four-window 908u + 9g(u) + {752,754,756,758}, u >= 1
// (other session's LINEAR_WINDOWS.md), counted up to X; plus collision statistics of F(u) = 908u + 9g(u):
// exact collisions and near-collisions |F(u) - F(u')| <= 6 (window overlaps).
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 g(i64 U){ if(U==0) return 1; int k=64-__builtin_clzll((u64)U); i64 s=1;
  for(int t=0;t<k;t++){ int bt=(U>>t)&1; if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
static int cmp(const void*a,const void*b){ i64 x=*(const i64*)a,y=*(const i64*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  i64 X=atoll(argv[1]); i64 U=X/908+2;
  unsigned char *hit=calloc(X/2+2,1);
  i64 *v=malloc(sizeof(i64)*U); i64 nv=0;
  for(i64 u=1;u<U;u++){ i64 F=908*u+9*g(u); v[nv++]=F; for(int d=752; d<=758; d+=2){ i64 n=F+d; if(n<=X) hit[n>>1]=1; } }
  i64 cntodd=0; for(i64 n=1;n<=X;n+=2) if(hit[n>>1]) cntodd++;
  qsort(v,nv,sizeof(i64),cmp);
  i64 exact=0, near=0; for(i64 i=1;i<nv;i++){ i64 d=v[i]-v[i-1]; if(d==0) exact++; if(d<=6) near++; }
  printf("X=%lld: the four-window family realizes %lld odd counts <= X, i.e. %.5f of the odd numbers (4/917 = %.5f, 4/908 = %.5f)\n",X,cntodd,(double)cntodd/(X/2),4.0/917,4.0/908);
  printf("  F(u), u < %lld: exact repeats %lld, adjacent gaps <= 6 (incl. repeats) %lld (fraction %.5f)\n",U,exact,near,(double)near/nv);
  return 0; }
