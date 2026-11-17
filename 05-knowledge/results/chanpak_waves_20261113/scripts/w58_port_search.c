/* Base skeleton: two triangles in series, T1 = {0,1,2} (0->1->2->0), T2 = {3,4,5} (3->4->5->3), T1 => T2,
   optionally with m0 singleton(s) between.  Tape slots: one vertex of each triangle (s1 in T1, s2 in T2).
   Add ports x, a with arbitrary out-sets to the base; port arc flipped.  Walk DP (THM-3121) on the blow-up skeleton.
   Report unit switches (W+ - W- = +-2 at (1,1) only) and their minus arrays; flag (2^p+1)^2 diagonals. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define VMAX 6
static int k, outm[12]; static int s1, s2;
static ll W[VMAX+1][VMAX+1];
static void walks(void){ memset(W,0,sizeof W);
  static ll dp[12][1<<12][VMAX+1][VMAX+1];
  int smask=0; for(int v=0;v<k;v++) if(v!=s1&&v!=s2) smask|=1<<v;
  for(int v=0;v<k;v++) for(int m=0;m<(1<<k);m++) if(!(m&~smask)) for(int a=0;a<=VMAX;a++) for(int b=0;b<=VMAX;b++) dp[v][m][a][b]=0;
  for(int v=0;v<k;v++){ if(v==s1) dp[v][0][1][0]=1; else if(v==s2) dp[v][0][0][1]=1; else dp[v][1<<v][0][0]=1; }
  for(int len=1; len<= (k-2)+2*VMAX; len++)
    for(int cur=0;cur<k;cur++) for(int m=0;m<(1<<k);m++){ if(m&~smask) continue; int pm=__builtin_popcount(m); if(pm>len) continue;
      for(int a=0;a<=VMAX;a++) for(int b=0;b<=VMAX;b++){ if(pm+a+b!=len) continue; ll c=dp[cur][m][a][b]; if(!c) continue;
          if(m==smask && a>=1 && b>=1) W[a][b]+=c;
          for(int nx=0;nx<k;nx++){ if(nx==cur || !(outm[cur]>>nx&1)) continue;
            if(nx==s1){ if(a<VMAX) dp[nx][m][a+1][b]+=c; }
            else if(nx==s2){ if(b<VMAX) dp[nx][m][a][b+1]+=c; }
            else if(!(m>>nx&1)) dp[nx][m|1<<nx][a][b]+=c; } } } }
static ll binom(int n,int r){ ll x=1; for(int i=0;i<r;i++) x=x*(n-i)/(i+1); return x; }
int main(int argc,char**argv){ int mid=atoi(argv[1]); /* number of singletons between the triangles */
  int k0=6+mid; k=k0+2; int X=k0, A=k0+1; long unit=0,hits=0,cfg=0; ll Wm[VMAX+1][VMAX+1],Wp[VMAX+1][VMAX+1];
  for(int t1=0;t1<3;t1++) for(int t2=3;t2<6;t2++){ s1=t1; s2=t2;
    for(int xo=0;xo<(1<<k0);xo++) for(int ao=0;ao<(1<<k0);ao++){
      /* build base: T1 0..2, mids 6..6+mid-1, T2 3..5 ; ordinal T1 -> mids -> T2 */
      int base[12]={0}; base[0]|=2; base[1]|=4; base[2]|=1; base[3]|=1<<4; base[4]|=1<<5; base[5]|=1<<3;
      for(int i=0;i<3;i++) for(int j=3;j<6;j++) base[i]|=1<<j;
      for(int q=0;q<mid;q++){ int v=6+q; for(int i=0;i<3;i++) base[i]|=1<<v; for(int j=3;j<6;j++) base[v]|=1<<j; for(int q2=q+1;q2<mid;q2++) base[v]|=1<<(6+q2); }
      for(int v=0;v<k0;v++) outm[v]=base[v];
      outm[X]=xo; outm[A]=ao; for(int v=0;v<k0;v++){ if(!(xo>>v&1)) outm[v]|=1<<X; if(!(ao>>v&1)) outm[v]|=1<<A; }
      cfg++;
      outm[X]|=1<<A; walks(); memcpy(Wm,W,sizeof W);
      outm[X]&=~(1<<A); outm[A]|=1<<X; walks(); memcpy(Wp,W,sizeof W);
      ll d=Wp[1][1]-Wm[1][1]; int ok=(d==2||d==-2); for(int a=1;a<=VMAX&&ok;a++) for(int b=1;b<=VMAX&&ok;b++) if(!(a==1&&b==1) && Wp[a][b]!=Wm[a][b]) ok=0;
      if(!ok) continue; unit++;
      ll (*L)[VMAX+1]=(d>0)?Wm:Wp; static ll c[VMAX+1][VMAX+1]; memset(c,0,sizeof c);
      for(int a=1;a<=VMAX;a++) for(int b=1;b<=VMAX;b++){ if(!L[a][b]) continue; for(int i=1;i<=a;i++) for(int j=1;j<=b;j++){ ll s=((a-i)%2?-1:1)*binom(a,i)*(((b-j)%2)?-1:1)*binom(b,j); c[i][j]+=L[a][b]*s; } }
      int maxb=0; for(int i=1;i<=VMAX;i++) for(int j=1;j<=VMAX;j++) if(c[i][j]){ if(i>maxb) maxb=i; if(j>maxb) maxb=j; }
      if(unit<=40 || maxb<=2) { printf("UNIT slots=(%d,%d) xo=%d ao=%d gap=%lld maxbase=%d c11=%lld c12=%lld c21=%lld c22=%lld c13=%lld c31=%lld c23=%lld c32=%lld c33=%lld\n",s1,s2,xo,ao,d,maxb,c[1][1],c[1][2],c[2][1],c[2][2],c[1][3],c[3][1],c[2][3],c[3][2],c[3][3]); }
      for(int d1=0;d1<=2;d1++) for(int d2=0;d2<=2;d2++){ ll kap[60]; memset(kap,0,sizeof kap);
        for(int i=1;i<=VMAX;i++) for(int j=1;j<=VMAX;j++) if(c[i][j]){ ll f=1; for(int t=0;t<d1;t++) f*=i; for(int t=0;t<d2;t++) f*=j; kap[i*j]+=c[i][j]*f; }
        int good=(kap[4]==1&&kap[2]==2&&(kap[1]==1||kap[1]==-1)); for(int B=3;B<60&&good;B++) if(B!=4&&kap[B]) good=0;
        if(good){ hits++; printf("HIT slots=(%d,%d) xo=%d ao=%d d=(%d,%d) eps=%lld\n",s1,s2,xo,ao,d1,d2,kap[1]-1); } } } }
  printf("mid=%d configs %ld unit %ld hits %ld\n",mid,cfg,unit,hits); return 0; }
