/* sieve_bounds.c: maximal 2-adic sieve (descent + backward-tree branch certificates, J=0) with the explicit threshold
   B of every certificate used: descent at t: n > c_t/(2^t-3^a);  branch m = (2^i x_s - cu)/3^b with
   x_s = (3^a n + c_s)/2^s:  m < n  iff  n > (2^i c_s - 2^s cu)/(2^s 3^b - 2^i 3^a).  Reports max B per depth. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned __int128 u128; typedef __int128 i128;
#define B3 39
static uint64_t P3[48]; static uint64_t INV2; static int KMAX;
static long double maxB[64]; static uint64_t alive[64];
static long double best_cert_B; static int found;
static inline uint64_t mulmod(uint64_t a,uint64_t b,uint64_t m){return (uint64_t)(((u128)a*b)%m);}
static int lt1(int e2,int e3){ /* 2^e2 3^e3 < 1 */
  if(e2<=0&&e3<=0) return !(e2==0&&e3==0); if(e2>=0&&e3>=0) return 0;
  if(e2<0){ if(-e2>=127) return 1; u128 p2=((u128)1)<<(-e2),p3=1; for(int i=0;i<e3;i++){p3*=3;if(p3>=p2)return 0;} return p3<p2; }
  if(e2>=127) return 0; u128 p2=((u128)1)<<e2,p3=1; for(int i=0;i<-e3;i++){p3*=3;if(p3>p2)return 1;} return p2<p3; }
static u128 CS; static int S_, A_;
/* node value v = (2^i x_s - cu)/3^b ; ratio 2^(i-s) 3^(a-b) */
static int bs(uint64_t z,int p,int i,int b,u128 cu){
  int e2=i-S_, e3=A_-b;
  if(lt1(e2,e3)){
    long double num=(long double)(((u128)1)<<i)*(long double)CS - (long double)(((u128)1)<<S_)*(long double)cu;
    long double den=ldexpl(powl(3.0L,b),S_) - ldexpl(powl(3.0L,A_),i);
    long double Bv = num>0? num/den : 0; best_cert_B=Bv; return 1; }
  if(!lt1(e2+p,e3-p)) return 0; if(p==0) return 0; int r=(int)(z%3); if(r==0) return 0;
  if(r==2){ uint64_t w=((2*z-1)/3)%P3[p-1]; if(bs(w,p-1,i+1,b+1,2*cu+(u128)P3[b])) return 1; }
  return bs((2*z)%P3[p],p,i+1,b,2*cu);
}
static void dfs(int t,int a,uint64_t Y,u128 c){
  if(t>=1){ u128 p2=((u128)1)<<t,p3=1; int d=1; for(int i=0;i<a;i++){p3*=3;if(p3>=p2){d=0;break;}}
    if(d&&p3<p2){ long double Bv=(long double)c/(long double)(p2-p3); if(Bv>maxB[t]) maxB[t]=Bv; return; } }
  int p=a<B3?a:B3; CS=c; S_=t; A_=a;
  if(bs(Y%P3[p],p,0,0,0)){ if(best_cert_B>maxB[t]) maxB[t]=best_cert_B; return; }
  alive[t]++; if(t==KMAX) return;
  dfs(t+1,a,mulmod(Y,INV2,P3[B3]),c);
  dfs(t+1,a+1,mulmod((3*Y+1)%P3[B3],INV2,P3[B3]),3*c+(((u128)1)<<t));
}
int main(int argc,char**argv){ KMAX=atoi(argv[1]); P3[0]=1; for(int i=1;i<=40;i++)P3[i]=P3[i-1]*3; INV2=(P3[B3]+1)/2;
  dfs(0,0,0,0); long double run=0;
  for(int t=1;t<=KMAX;t++){ if(maxB[t]>run) run=maxB[t]; printf("%2d alive %10llu  max certificate threshold at depth t: 2^%.2Lf   running max 2^%.2Lf\n",t,(unsigned long long)alive[t], maxB[t]>0?log2l(maxB[t]):0.0L, run>0?log2l(run):0.0L); }
}
