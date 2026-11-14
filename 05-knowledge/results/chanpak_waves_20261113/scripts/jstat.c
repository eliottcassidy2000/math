#include <stdio.h>
#include <string.h>
#include <stdlib.h>
typedef long long ll;
static ll Hof(int*x,int*y,int len,int e){
  ll c[2][3][3][2]; memset(c,0,sizeof c); c[0][0][0][0]=1; c[1][0][0][0]=1;
  for(int i=0;i<len;i++){ ll n[2][3][3][2]; memset(n,0,sizeof n);
    for(int o=0;o<2;o++){ int h1=o?y[i]:x[i], h2=o?x[i]:y[i];
      for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ ll v=c[o][a][b][d]; if(!v) continue;
        n[o][h1?2:1][b][d]+=v; if(b||h1) n[o][a][h2?2:1][d]+=v; if(d||h2) n[o][a][b][1]+=v; } }
    memcpy(c,n,sizeof c); }
  ll H=0; for(int o=0;o<2;o++){ int h12=o?!e:e; for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++) if(!(a==2||b==2||(b==0&&!h12))) H+=c[o][a][b][d]; }
  return H; }
int main(void){
  int X[20],Y[20]; ll tot=0, j2=0, j0=0, small=0; ll hist[9]={0};
  for(int L=1;L<=10;L++) for(unsigned long long code=0;code<(1ULL<<(2*L));code++){
    for(int i=0;i<L;i++){X[i]=(code>>(2*i))&1; Y[i]=(code>>(2*i+1))&1;}
    ll J=Hof(X,Y,L,1)-Hof(X,Y,L,0); ll a=J<0?-J:J; tot++;
    if(a==2) j2++; if(a==0) j0++; if(a<=16) small++; if(a<=16 && a%2==0) hist[a/2]++; }
  printf("word pairs (L<=10): %lld; |J|=0: %.1f%%, |J|=2: %.1f%%, |J|<=16: %.1f%%\n",tot,100.0*j0/tot,100.0*j2/tot,100.0*small/tot);
  printf("|J| histogram (0,2,...,16):"); for(int i=0;i<9;i++) printf(" %lld",hist[i]); printf("\n"); return 0; }
