/* Skeleton blow-up switch search.  Skeleton Q (iso reps, order k <= 7 from reps.txt), two tape slots s1 != s2 (blown up to TT_r1, TT_r2),
   flip one arc between two non-slot vertices.  H(Q[...]) = sum over surjective walks of prod Surj(r_i, visits_i) (THM-3121).
   Condition: W(+) - W(-) = 2 at visit vector (1,1) and 0 elsewhere (gap exactly 2 for all r1,r2 >= 1).
   Report the lower member's exponential array c_ij (bases i^r1 j^r2) and test diagonals for (2^p+1)^2 + eps, eps in {0,-2}. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define VMAX 6
static int k, outm[8]; static int s1, s2;
static ll W[VMAX+1][VMAX+1];
static ll memo[8][1<<8][VMAX+1][VMAX+1]; static char seen[8][1<<8][VMAX+1][VMAX+1];
static int sing_mask_full;
/* count walks: state (cur, mask of visited singletons, v1, v2) -> number of completions; we just do forward DP */
static void walks(void){ memset(W,0,sizeof W);
  static ll dp[8][1<<8][VMAX+1][VMAX+1]; memset(dp,0,sizeof dp);
  int smask=0; for(int v=0;v<k;v++) if(v!=s1&&v!=s2) smask|=1<<v; sing_mask_full=smask;
  for(int v=0;v<k;v++){ if(v==s1) dp[v][0][1][0]=1; else if(v==s2) dp[v][0][0][1]=1; else dp[v][1<<v][0][0]=1; }
  /* iterate by total length; simple relaxation in order of (popcount(mask)+v1+v2) */
  for(int len=1; len<= (k-2)+2*VMAX; len++){
    for(int cur=0;cur<k;cur++) for(int m=0;m<(1<<k);m++){ if(m&~smask) continue; for(int a=0;a<=VMAX;a++) for(int b=0;b<=VMAX;b++){ ll c=dp[cur][m][a][b]; if(!c) continue; if(__builtin_popcount(m)+a+b!=len) continue;
          if(m==smask && a>=1 && b>=1) W[a][b]+=c;
          for(int nx=0;nx<k;nx++){ if(nx==cur || !(outm[cur]>>nx&1)) continue;
            if(nx==s1){ if(a<VMAX) dp[nx][m][a+1][b]+=c; }
            else if(nx==s2){ if(b<VMAX) dp[nx][m][a][b+1]+=c; }
            else if(!(m>>nx&1)) dp[nx][m|1<<nx][a][b]+=c; } } } } }
static ll binom(int n,int r){ ll x=1; for(int i=0;i<r;i++) x=x*(n-i)/(i+1); return x; }
int main(int argc,char**argv){ int KMIN=atoi(argv[1]), KMAX=atoi(argv[2]); FILE*fi=fopen("w57_reps.txt","r"); char line[256]; long cfg=0, unit=0, hits=0;
  ll Wm[VMAX+1][VMAX+1], Wp[VMAX+1][VMAX+1];
  while(fgets(line,sizeof line,fi)){ char*p=line; k=strtol(p,&p,10); int o[8]; for(int i=0;i<k;i++) o[i]=strtol(p,&p,10); if(k<KMIN||k>KMAX) continue;
    for(s1=0;s1<k;s1++) for(s2=0;s2<k;s2++){ if(s1==s2) continue;
      for(int u=0;u<k;u++) for(int w=0;w<k;w++){ if(u==s1||u==s2||w==s1||w==s2||u==w) continue; if(!(o[u]>>w&1)) continue; /* arc u->w, flip it */
        cfg++;
        memcpy(outm,o,sizeof(int)*k); walks(); memcpy(Wm,W,sizeof W);
        outm[u]&=~(1<<w); outm[w]|=1<<u; walks(); memcpy(Wp,W,sizeof W);
        int ok=1; ll d11=Wp[1][1]-Wm[1][1]; if(d11!=2 && d11!=-2) ok=0;
        for(int a=1;a<=VMAX&&ok;a++) for(int b=1;b<=VMAX&&ok;b++){ if(a==1&&b==1) continue; if(Wp[a][b]!=Wm[a][b]) ok=0; }
        if(!ok) continue; unit++;
        ll (*L)[VMAX+1] = (d11>0)?Wm:Wp; /* lower member */
        /* exponential array c[i][j], i,j in 1..VMAX: Surj(r,v)=sum_i (-1)^(v-i) C(v,i) i^r */
        static ll c[VMAX+1][VMAX+1]; memset(c,0,sizeof c);
        for(int a=1;a<=VMAX;a++) for(int b=1;b<=VMAX;b++){ if(!L[a][b]) continue; for(int i=1;i<=a;i++) for(int j=1;j<=b;j++){ ll s=((a-i)%2?-1:1)*binom(a,i)*(((b-j)%2)?-1:1)*binom(b,j); c[i][j]+=L[a][b]*s; } }
        { int maxb=0; for(int i=1;i<=VMAX;i++) for(int j=1;j<=VMAX;j++) if(c[i][j]){ if(i>maxb) maxb=i; if(j>maxb) maxb=j; }
          printf("UNIT k=%d slots=(%d,%d) maxbase=%d c22=%lld c12=%lld c21=%lld c11=%lld c33=%lld c13=%lld c31=%lld c23=%lld c32=%lld\n",k,s1,s2,maxb,c[2][2],c[1][2],c[2][1],c[1][1],c[3][3],c[1][3],c[3][1],c[2][3],c[3][2]); }
        /* diagonal tests r1=p+d1, r2=p+d2, d in {0,1,2}: kappa_B = sum_{ij=B} c_ij i^d1 j^d2 ; need kappa_4=1,kappa_2=2,kappa_1 in {1,-1}, others 0 */
        for(int d1=0;d1<=2;d1++) for(int d2=0;d2<=2;d2++){ ll kap[50]; memset(kap,0,sizeof kap);
            for(int i=1;i<=VMAX;i++) for(int j=1;j<=VMAX;j++) if(c[i][j]){ ll f=1; for(int t=0;t<d1;t++) f*=i; for(int t=0;t<d2;t++) f*=j; kap[i*j]+=c[i][j]*f; }
            int good = (kap[4]==1 && kap[2]==2 && (kap[1]==1||kap[1]==-1)); for(int B=3;B<50&&good;B++) if(B!=4 && kap[B]) good=0;
            if(good){ hits++; printf("HIT k=%d Q=",k); for(int i=0;i<k;i++) printf("%d,",o[i]); printf(" slots=(%d,%d) flip %d->%d  d=(%d,%d) eps=%lld\n",s1,s2,u,w,d1,d2,kap[1]-1); fflush(stdout); } }
      } } }
  printf("configs %ld, unit switches (gap exactly 2 for all tape lengths) %ld, diagonal hits %ld\n",cfg,unit,hits); return 0; }
