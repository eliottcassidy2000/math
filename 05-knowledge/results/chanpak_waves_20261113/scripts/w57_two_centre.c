/* Wave 57 recomputation of known families (credit THM-4876 (6) one centre; Two-Hub Conjecture, wave 5).
   Usage: ./w57_one_centre MAXLEN LIMIT ; ./w57_two_centre MAXLEN LIMIT */
/* Two contact vertices x,y on the toothpick TT_n.  X_i=1 iff x->i, Y_i=1 iff y->i, e=1 iff x->y.
   H = sum over the two orders (x first / y first) of #ordered 3-block set partitions (R1,R2,R3) of [n]
   with: R1 empty or last(R1) -> first apex; R2 empty ? (first apex -> second apex) : (first apex -> min R2, max R2 -> second apex);
   R3 empty or second apex -> min R3.   DP over i with block states {empty, ok, bad}. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n;
static u64 count_order(const int*A,const int*B,int ab){ /* A = first apex contact word, B = second; ab = 1 iff first apex -> second */
  u64 dp[27]={0}, nd[27]; dp[0]=1; /* state = s1 + 3 s2 + 9 s3, s in {0 empty,1 ok,2 bad} */
  for(int i=0;i<n;i++){ memset(nd,0,sizeof nd);
    for(int st=0;st<27;st++){ u64 c=dp[st]; if(!c) continue; int s1=st%3, s2=st/3%3, s3=st/9;
      /* to R1: max cond: i -> A, i.e. A[i]==0 */
      { int t1 = (A[i]==0)?1:2; nd[t1+3*s2+9*s3]+=c; }
      /* to R2: min cond A[i]==1 (A -> i) when empty; max cond: i -> B, B[i]==0 */
      if(s2!=0 || A[i]==1){ int t2=(B[i]==0)?1:2; nd[s1+3*t2+9*s3]+=c; }
      /* to R3: min cond B[i]==1 when empty; no max cond */
      if(s3!=0 || B[i]==1){ nd[s1+3*s2+9*1]+=c; } }
    memcpy(dp,nd,sizeof dp); }
  u64 tot=0; for(int st=0;st<27;st++){ int s1=st%3,s2=st/3%3,s3=st/9; if(s1==2||s2==2||s3==2) continue; if(s2==0 && !ab) continue; tot+=dp[st]; }
  return tot; }
int main(int argc,char**argv){ int L=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10); unsigned char*hit=calloc(LIM/2+2,1);
  int X[32],Y[32];
  for(n=1;n<=L;n++){ u64 W=1ULL<<n;
    for(u64 wx=0;wx<W;wx++) for(u64 wy=wx;wy<W;wy++) for(int e=0;e<2;e++){ /* x<->y symmetric: wy>=wx suffices */
      for(int i=0;i<n;i++){ X[i]=wx>>i&1; Y[i]=wy>>i&1; }
      u64 h=count_order(X,Y,e)+count_order(Y,X,!e); if(h<=LIM && (h&1)) hit[h/2]=1; }
    u64 miss=0; for(u64 N=1;N<=LIM;N+=2) if(!hit[N/2]) miss++;
    printf("n<=%d (order n+2): odd N <= %llu missed %llu:",n,LIM,miss); int c=0; for(u64 N=1;N<=LIM&&c<30;N+=2) if(!hit[N/2]){ printf(" %llu",N); c++; } printf("\n"); fflush(stdout); }
  return 0; }
