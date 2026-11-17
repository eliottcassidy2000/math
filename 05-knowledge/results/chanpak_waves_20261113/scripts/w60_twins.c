/* Profile twins: classes of the same order with equal ordered-cover counts c_r for all r >= 2 but different c_1 = H.
   Substituting either at a vertex X of any Q changes H(Q[.]) by (c_1 - c_1')·H(Q) (THM-3121): a module switch. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 Hs[256];
static void allH(int n,const int*o){ static u64 dp[256][8]; int N=1<<n; for(int S=0;S<N;S++) for(int v=0;v<n;v++) dp[S][v]=0; for(int v=0;v<n;v++) dp[1<<v][v]=1;
  for(int S=1;S<N;S++) for(int v=0;v<n;v++) if(dp[S][v]) for(int w=0;w<n;w++) if(!(S>>w&1)&&(o[v]>>w&1)) dp[S|1<<w][w]+=dp[S][v];
  Hs[0]=1; for(int S=1;S<N;S++){ u64 h=0; for(int v=0;v<n;v++) h+=dp[S][v]; Hs[S]=h; } }
/* ordered r-covers via subset convolution: f_r(S) = sum over first block A (nonempty, contains min of S? no: ordered) */
static void covers(int n,u64*c){ static u64 f[9][256]; int N=1<<n; for(int r=0;r<=n;r++) for(int S=0;S<N;S++) f[r][S]=0; f[0][0]=1;
  for(int r=1;r<=n;r++) for(int S=1;S<N;S++){ u64 t=0; for(int A=S;A;A=(A-1)&S) t+=Hs[A]*f[r-1][S&~A]; f[r][S]=t; }
  for(int r=1;r<=n;r++) c[r]=f[r][N-1]; }
typedef struct { int n; u64 c[9]; int o[8]; } Rec;
static int cmp(const void*a,const void*b){ const Rec*x=a,*y=b; if(x->n!=y->n) return x->n-y->n; for(int r=2;r<=x->n;r++) if(x->c[r]!=y->c[r]) return x->c[r]<y->c[r]?-1:1; return x->c[1]<y->c[1]?-1:(x->c[1]>y->c[1]); }
int main(void){ FILE*fi=fopen("w59_reps8.txt","r"); char line[256]; Rec*R=malloc(sizeof(Rec)*7500); int k=0;
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); int o[8]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10); allH(n,o); R[k].n=n; covers(n,R[k].c); memcpy(R[k].o,o,sizeof o); k++; }
  qsort(R,k,sizeof(Rec),cmp); int groups=0, unit=0;
  for(int i=0;i<k;){ int j=i; while(j<k && R[j].n==R[i].n){ int same=1; for(int r=2;r<=R[i].n;r++) if(R[j].c[r]!=R[i].c[r]) same=0; if(!same) break; j++; }
    /* distinct c_1 within group */
    u64 vals[64]; int nv=0; for(int t=i;t<j;t++){ int seen=0; for(int s=0;s<nv;s++) if(vals[s]==R[t].c[1]) seen=1; if(!seen&&nv<64) vals[nv++]=R[t].c[1]; }
    if(nv>1){ groups++; u64 mn=vals[0],mx=vals[0]; for(int s=0;s<nv;s++){ if(vals[s]<mn) mn=vals[s]; if(vals[s]>mx) mx=vals[s]; }
      int hasunit=0; for(int a=0;a<nv;a++) for(int b=0;b<nv;b++) if(vals[a]+2==vals[b]) hasunit=1; if(hasunit) unit++;
      if(groups<=12||hasunit){ printf("order %d twins: H values",R[i].n); for(int s=0;s<nv;s++) printf(" %llu",vals[s]); printf("  (c_2 = %llu)%s\n",R[i].c[2],hasunit?"  UNIT":""); } }
    i=j; }
  printf("profile-twin groups (orders <= 8): %d; groups containing an H-difference of exactly 2: %d\n",groups,unit); return 0; }
