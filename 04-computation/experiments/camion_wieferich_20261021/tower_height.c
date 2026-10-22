/* tower_height.c -- Camion tower data via the higher Kummer congruences (THM: for m >= 2K-1 the sequence
 * A(j) = eps^-j D(p + j(p-1)) is, mod p^K, a polynomial of degree <= 2K-2 in j).  Computes s(m) mod p^K for
 * m <= p + (2K-2)(p-1), interpolates, and evaluates at j_a = (p^a - p)/(p-1) (giving D(p^a) mod p^K) and at the
 * p-adic limit j_inf = -p/(p-1) (giving Lambda_p = lim eps^a D(p^a) mod p^K, up to the factor eps^(j)).
 * Prints v_p of D(p^a) (a = 1..K+2) and of Lambda_p; h(p) = v_p(Lambda_p) when < K.  Requires p > 2K.
 * Optional direct check of D(p^2) mod p^K when p^2 is small.  mac-mini-2026-10-21-fanopsp.  Usage: ./tower_height p K [direct] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128; typedef __int128 i128;
static u64 Q;
static u64 mulm(u64 a,u64 b){ return (u64)((u128)a*b%Q); }
static u64 pw(u64 b,u64 e){ u64 r=1%Q; b%=Q; while(e){ if(e&1) r=mulm(r,b); b=mulm(b,b); e>>=1; } return r; }
static u64 inv(u64 a){ /* a unit mod Q = p^K: extended Euclid */ i128 t=0,nt=1,r=Q,nr=a%Q; while(nr){ i128 qq=r/nr, tmp; tmp=t-qq*nt; t=nt; nt=tmp; tmp=r-qq*nr; r=nr; nr=tmp; } if(t<0) t+=Q; return (u64)t; }
static int vp(u64 x,u64 p,int K){ if(x==0) return K; int c=0; while(x%p==0){ x/=p; c++; } return c; }
static u64 *s_upto(long M){ u64 *T=malloc((M+1)*8),*s=malloc((M+1)*8),*row=malloc((M+2)*8);
  for(long j=0;j<=M;j++) T[j]=pw(2,(u64)j*(j-1)/2); row[0]=1;
  for(long m=1;m<=M;m++){ for(long i=m-1;i>=1;i--){ row[i]+=row[i-1]; if(row[i]>=Q) row[i]-=Q; } row[m]=1;
    u128 acc=0; for(long i=1;i<m;i++){ acc+=(u128)mulm(row[i],s[i])*T[m-i]; if((i&15)==0) acc%=Q; } acc%=Q; s[m]=(T[m]+Q-(u64)acc)%Q; }
  free(T); free(row); return s; }
static u64 cfun(long m){ return (m&1)? pw(2,(m-1)/2) : 0; }
int main(int argc,char**argv){ u64 p=atoll(argv[1]); int K=atoi(argv[2]); int direct=(argc>3);
  Q=1; for(int i=0;i<K;i++) Q*=p; int eps=(p%8==1||p%8==7)?1:-1; int r=2*K-1; /* points j = 0..2K-2 */
  long M=(long)p+(long)(2*K-2)*(long)(p-1);
  u64 *s=s_upto(M);
  u64 A[64]; for(int j=0;j<r;j++){ long m=(long)p+(long)j*(long)(p-1); u64 D=(s[m]+Q-cfun(m))%Q; A[j]=(eps==1||j%2==0)?D:(Q-D)%Q; }
  /* forward differences */
  u64 dif[64]; u64 tmp[64]; memcpy(tmp,A,sizeof(u64)*r); for(int i=0;i<r;i++){ dif[i]=tmp[0]; for(int j=0;j<r-1-i;j++) tmp[j]=(tmp[j+1]+Q-tmp[j])%Q; }
  /* evaluate Newton series at z (given mod Q as an element; binomial C(z,i) = prod (z-l)/i!) */
  #define EVAL(zv,res) { u64 acc=0, b=1%Q; for(int i=0;i<r;i++){ acc=(acc+mulm(b,dif[i]))%Q; b=mulm(mulm(b,((zv)+Q-(u64)i%Q)%Q),inv((u64)(i+1))); } res=acc; }
  printf("p=%llu K=%d eps=%+d  (s computed to m=%ld)\n",p,K,eps,M);
  u64 pa=p; for(int a=1;a<=K+2;a++){ /* j_a = (p^a - p)/(p-1) mod Q, computed exactly via the integer for small a, else via p-adic formula */
    u64 ja; { /* (p^a - p)/(p-1) = p + p^2 + ... + p^(a-1)  (a>=2), 0 for a=1 */ u64 acc=0, t=p%Q; for(int i=1;i<a;i++){ acc=(acc+t)%Q; t=mulm(t,p); } ja=acc; }
    u64 val; EVAL(ja,val); /* A(j) = eps^-j D: D(p^a) = eps^(j_a) * A(j_a); parity of j_a = (a-1) mod 2 */
    long par=(a-1)&1; u64 Dv=(eps==1||par==0)?val:(Q-val)%Q;
    printf("  a=%d: v_p(D(p^%d)) = %s%d\n",a,a,(vp(Dv,p,K)==K)?">=":"",vp(Dv,p,K)); pa*=p; }
  /* p-adic limit point j_inf = -p/(p-1) */
  u64 jinf=(Q-mulm(p%Q,inv((p-1)%Q)))%Q; u64 lam; EVAL(jinf,lam);
  printf("  Lambda_p (limit, up to sign): v_p = %s%d   => h(p) %s %d\n",(vp(lam,p,K)==K)?">=":"",vp(lam,p,K),(vp(lam,p,K)<K)?"=":">=",vp(lam,p,K));
  if(direct){ long m2=(long)(p*p); u64 *s2=s_upto(m2); u64 Dd=(s2[m2]+Q-cfun(m2))%Q; u64 ja=p%Q; u64 val; EVAL(ja,val); u64 Dv=(eps==1)?val:(Q-val)%Q;
    printf("  direct check D(p^2) mod p^%d: direct=%llu interpolated=%llu %s\n",K,Dd,Dv,(Dd==Dv)?"MATCH":"MISMATCH"); }
  return 0; }
