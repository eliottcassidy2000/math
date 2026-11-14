// Coverage of the hub-cut family: N = H(x) + 2^j L(x,k)  (j>=1), plus N = H(x), plus hub-root rays H + 2^j (1+X).
// H(x)=1+2W(x) one-centre count; L(x,k) = port of spine vertex k (DP verified against brute force).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 LIM; static int MMAX; static unsigned char *cov;
static int x[64];
static void mark(u64 N){ if(N<=LIM && (N&1)) cov[N>>1]=1; }
static u64 Lroot(int m,int k){
  // case A: u in P; state index: pa(-1,0,1)->0,1,2 ; pb 0/1 ; q 0/1  => 12 states
  u64 A[3][2][2]; memset(A,0,sizeof A); A[0][0][0]=1;
  for(int p=0;p<m;p++){ if(p==k) continue;
    u64 N[3][2][2]; memset(N,0,sizeof N);
    for(int pa=0;pa<3;pa++)for(int pb=0;pb<2;pb++)for(int q=0;q<2;q++){ u64 c=A[pa][pb][q]; if(!c) continue;
      N[1+x[p]][pb][q]+=c;                         // to Pa
      if(pb||x[p]==1) N[pa][1][q]+=c;              // to Pb
      if(q||p>k) N[pa][pb][1]+=c;                  // to Q
    }
    memcpy(A,N,sizeof A);
  }
  u64 a=0; for(int pb=0;pb<2;pb++)for(int q=0;q<2;q++) a+=A[0][pb][q]+A[1][pb][q];
  u64 B[3][2]; memset(B,0,sizeof B); B[0][0]=1;
  for(int p=0;p<m;p++){ if(p==k) continue;
    u64 N[3][2]; memset(N,0,sizeof N);
    for(int qa=0;qa<3;qa++)for(int qb=0;qb<2;qb++){ u64 c=B[qa][qb]; if(!c) continue;
      N[qa][qb]+=c;                                 // to P
      if(qa!=0||p>k) N[1+x[p]][qb]+=c;              // to Qa
      if(qb||x[p]==1) N[qa][1]+=c;                  // to Qb
    }
    memcpy(B,N,sizeof B);
  }
  u64 b=B[1][0]+B[1][1]; if(x[k]==0) b+=B[0][0]+B[0][1];
  return a+b;
}
static long long nwords=0;
static void dfs(int m,u64 w,u64 Xm,u64 Xm1){
  u64 H=1+2*w;
  if(H>LIM) return;
  nwords++;
  mark(H);
  for(u64 j=1, L=1+Xm; ; j++){ if(L> (LIM>>j)) break; mark(H+(L<<j)); }
  if(m>0 && 3*H<=LIM){
    for(int k=0;k<m;k++){ u64 L=Lroot(m,k); for(int j=1;j<62;j++){ if(L>(LIM>>j)) break; mark(H+(L<<j)); } }
  }
  if(m>=MMAX) return;
  if(Xm==0){ x[m]=1; dfs(m+1,w,1,0); return; }
  x[m]=0; dfs(m+1,w+(Xm-Xm1),2*Xm,Xm);
  x[m]=1; dfs(m+1,w,2*Xm+1,Xm);
}
int main(int argc,char**argv){
  LIM=strtoull(argv[1],0,10); MMAX=atoi(argv[2]);
  cov=calloc(LIM/2+2,1);
  dfs(0,0,0,0);
  long long miss=0; u64 last=0; int shown=0;
  printf("words visited %lld\n",nwords);
  printf("misses:");
  for(u64 N=1;N<=LIM;N+=2) if(!cov[N>>1]){ miss++; last=N; if(shown<70){printf(" %llu",N);shown++;} }
  printf("\ntotal misses %lld of %llu odd <= %llu ; largest miss %llu\n",miss,(LIM+1)/2,LIM,last);
  return 0;
}
