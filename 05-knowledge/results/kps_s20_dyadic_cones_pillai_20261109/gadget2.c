/* gadget.c: for every labelled gadget (G,u,z) on g vertices (g<=GMAX), compute the 3x3 one-arc join matrix
   M(G,u,z): (h,e,R)(N,z) = M (h,e,R)(T,v), where N = T => G with u -> v reversed.
   Computed by brute force on three bases with independent states and checked on a fourth.
   Output line: g Gcode u z  m00 m01 m02 m10 m11 m12 m20 m21 m22   (Gcode = arc bits i<j: bit set if i->j) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll Hs[1<<14], endd[(1<<14)*14], st[(1<<14)*14];
static int n; static unsigned out[16];
static void dp(void){
  unsigned N=1u<<n; Hs[0]=1;
  for(unsigned S=1;S<N;S++){ ll tot=0;
    for(int v=0;v<n;v++) if(S>>v&1){ unsigned P=S&~(1u<<v); ll a=0,b=0;
      if(!P){a=1;b=1;} else for(int w=0;w<n;w++) if(P>>w&1){ if(out[w]>>v&1) a+=endd[P*n+w]; if(out[v]>>w&1) b+=st[P*n+w]; }
      endd[S*n+v]=a; st[S*n+v]=b; tot+=a; }
    Hs[S]=tot; }
}
static void state(int z, ll *h, ll *e, ll *R){
  unsigned N=1u<<n, F=N-1; *h=Hs[F]; *e=endd[F*n+z]; ll r=0;
  for(unsigned S=0;S<N;S++) if(S>>z&1) r+=st[S*n+z]*Hs[F&~S];
  *R=r;
}
/* bases: (out-masks, size, root) */
static unsigned B0[]={0}; /* single vertex: (1,1,1) */
static unsigned B1[]={2,0}; /* 0->1, root 0: (1,0,2) */
static unsigned B2[]={2,4,1}; /* C3 root 0: (3,1,3) */
static unsigned B3[]={6,4,0,7}; /* TT3 + ... root 3: we'll compute */
int main(int argc,char**argv){
  int GMAX=atoi(argv[1]);
  unsigned *bases[4]={B0,B1,B2,B3}; int bsz[4]={1,2,3,4}; int broot[4]={0,0,0,3};
  ll bs[4][3];
  for(int b=0;b<4;b++){ n=bsz[b]; for(int i=0;i<n;i++) out[i]=bases[b][i]; dp(); state(broot[b],&bs[b][0],&bs[b][1],&bs[b][2]); }
  /* solve: X = [bs0;bs1;bs2] (rows), M^T = X^-1 Y. use integer Cramer (det X) */
  ll X[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++)X[i][j]=bs[i][j];
  ll det=X[0][0]*(X[1][1]*X[2][2]-X[1][2]*X[2][1])-X[0][1]*(X[1][0]*X[2][2]-X[1][2]*X[2][0])+X[0][2]*(X[1][0]*X[2][1]-X[1][1]*X[2][0]);
  ll adj[3][3]; /* adj = det * X^-1 */
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){ int r0=(j+1)%3,r1=(j+2)%3,c0=(i+1)%3,c1=(i+2)%3; adj[i][j]=X[r0][c0]*X[r1][c1]-X[r0][c1]*X[r1][c0]; }
  fprintf(stderr,"base states: (%lld %lld %lld) (%lld %lld %lld) (%lld %lld %lld) (%lld %lld %lld) det %lld\n",
    bs[0][0],bs[0][1],bs[0][2],bs[1][0],bs[1][1],bs[1][2],bs[2][0],bs[2][1],bs[2][2],bs[3][0],bs[3][1],bs[3][2],det);
  long bad=0,tot=0;
  for(int g=1;g<=GMAX;g++){
    char fn[512]; sprintf(fn,"%s/cls_%d.txt",argv[2],g); FILE*cf=fopen(fn,"r"); char line[4096];
    while(fgets(line,sizeof line,cf)){ unsigned long long code; if(sscanf(line,"%llu",&code)!=1) continue;
      unsigned gout[16]={0}; int k=0;
      for(int i=0;i<g;i++)for(int j=i+1;j<g;j++){ if(code>>k&1) gout[i]|=1u<<j; else gout[j]|=1u<<i; k++; }
      for(int u=0;u<g;u++) for(int z=0;z<g;z++){
        ll Y[4][3];
        for(int b=0;b<4;b++){
          int tb=bsz[b]; n=tb+g;
          for(int i=0;i<tb;i++){ out[i]=bases[b][i]; for(int j=0;j<g;j++) out[i]|=1u<<(tb+j); }
          for(int j=0;j<g;j++) out[tb+j]=gout[j]<<tb;
          out[broot[b]]&=~(1u<<(tb+u)); out[tb+u]|=1u<<broot[b];
          dp(); state(tb+z,&Y[b][0],&Y[b][1],&Y[b][2]);
        }
        ll M[3][3]; int okdiv=1;
        for(int r=0;r<3;r++)for(int c=0;c<3;c++){ ll s=0; for(int i=0;i<3;i++) s+=adj[c][i]*Y[i][r]; if(s%det) okdiv=0; M[r][c]=s/det; }
        int ok=okdiv; for(int r=0;r<3;r++){ ll s=0; for(int c=0;c<3;c++) s+=M[r][c]*bs[3][c]; if(s!=Y[3][r]) ok=0; }
        tot++; if(!ok) bad++;
        printf("%d %llu %d %d",g,code,u,z); for(int r=0;r<3;r++)for(int c=0;c<3;c++) printf(" %lld",M[r][c]); putchar(10);
      }
    }
    fclose(cf);
  }
  fprintf(stderr,"gadgets %ld, linearity failures %ld\n",tot,bad);
  return 0;
}
