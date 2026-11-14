// Multiplicity of F(u) = A*u + B*g(u) over 1 <= u <= U: sort the values and count repeats.
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 g(i64 U){ if(U==0) return 1; int k=64-__builtin_clzll((u64)U); i64 s=1;
  for(int t=0;t<k;t++){ int bt=(U>>t)&1; if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
static int cmp(const void*a,const void*b){ i64 x=*(const i64*)a,y=*(const i64*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  i64 U=atoll(argv[1]); i64 A=argc>2?atoll(argv[2]):908, B=argc>3?atoll(argv[3]):9;
  i64 *v=malloc(sizeof(i64)*(U+1)); for(i64 u=1;u<=U;u++) v[u-1]=A*u+B*g(u);
  qsort(v,U,sizeof(i64),cmp);
  i64 distinct=0,maxm=0,m=1,repeated=0;
  for(i64 i=0;i<U;i++){ if(i>0 && v[i]==v[i-1]) m++; else { if(i>0){ if(m>maxm) maxm=m; if(m>1) repeated+=m; } distinct++; m=1; } }
  if(m>maxm) maxm=m; if(m>1) repeated+=m;
  printf("F(u)=%lld u + %lld g(u), 1<=u<=%lld: distinct values %lld, max multiplicity %lld, arguments in repeated values %lld\n",A,B,U,distinct,maxm,repeated);
  return 0; }
