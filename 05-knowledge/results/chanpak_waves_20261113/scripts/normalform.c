#include <stdio.h>
#include <string.h>
#include <stdlib.h>
typedef unsigned long long u64;
static u64 Hof(const int*x,const int*y,int len,int e){
  u64 c[2][3][3][2]; memset(c,0,sizeof c); c[0][0][0][0]=1; c[1][0][0][0]=1;
  for(int i=0;i<len;i++){ u64 n[2][3][3][2]; memset(n,0,sizeof n);
    for(int o=0;o<2;o++){ int h1=o?y[i]:x[i], h2=o?x[i]:y[i];
      for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ u64 v=c[o][a][b][d]; if(!v) continue;
        n[o][h1?2:1][b][d]+=v; if(b||h1) n[o][a][h2?2:1][d]+=v; if(d||h2) n[o][a][b][1]+=v; } }
    memcpy(c,n,sizeof c); }
  u64 H=0; for(int o=0;o<2;o++){ int h12=o?!e:e; for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++) if(!(a==2||b==2||(b==0&&!h12))) H+=c[o][a][b][d]; }
  return H; }
#define MAXV 200000
static int nfL[MAXV]; static u64 nfcode[MAXV]; static int nfe[MAXV];
int main(void){
  int X[16],Y[16]; for(int i=0;i<MAXV;i++) nfL[i]=-1;
  for(int L=0;L<=10;L++) for(u64 code=0;code<(1ULL<<(2*L));code++){
    for(int i=0;i<L;i++){X[i]=(code>>(2*(L-1-i)))&1; Y[i]=(code>>(2*(L-1-i)+1))&1;}   // first letter = most significant (shortlex)
    for(int e=0;e<2;e++){ u64 h=Hof(X,Y,L,e); if(h<MAXV && nfL[h]<0){ nfL[h]=L; nfcode[h]=code; nfe[h]=e; } } }
  int bound=29523; /* = (3^10-1)/2 */ int tot=0, here=0, here_any=0; int firstbad[12]; int nb=0;
  for(int N=3;N<=bound;N+=2){ if(nfL[N]<1) continue; tot++;
    int L=nfL[N]; u64 code=nfcode[N]>>2; /* drop last letter (least significant) */
    for(int i=0;i<L-1;i++){X[i]=(code>>(2*(L-2-i)))&1; Y[i]=(code>>(2*(L-2-i)+1))&1;}
    int ok=0, okany=0;
    for(int e=0;e<2;e++){ u64 h=Hof(X,Y,L-1,e); if(h<MAXV && nfL[h]==L-1 && nfcode[h]==code && nfe[h]==e) ok=1; if(h<MAXV && nfL[h]==L-1 && nfcode[h]==code) okany=1; }
    if(ok) here++; if(okany) here_any++; else if(nb<12) firstbad[nb++]=N; }
  printf("odd N in [3,%d] with normal forms: %d; truncation is a normal form (same e): %d (%.1f%%); (any e): %d (%.1f%%)\n",bound,tot,here,100.0*here/tot,here_any,100.0*here_any/tot);
  printf("first non-hereditary N:"); for(int i=0;i<nb;i++) printf(" %d",firstbad[i]); printf("\n"); return 0; }
