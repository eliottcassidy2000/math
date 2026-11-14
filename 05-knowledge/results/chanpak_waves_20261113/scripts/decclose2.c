// P25: is every odd n (23..B, n!=21) reached by a one-vertex extension of a realizer of m with n = 2m + 3^k ?
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
typedef unsigned long long u64;
static int X[32],Y[32],L;
static u64 Hof(int*x,int*y,int len,int e){
  u64 c[2][3][3][2]; memset(c,0,sizeof c); c[0][0][0][0]=1; c[1][0][0][0]=1;
  for(int i=0;i<len;i++){ u64 n[2][3][3][2]; memset(n,0,sizeof n);
    for(int o=0;o<2;o++){ int h1=o?y[i]:x[i], h2=o?x[i]:y[i];
      for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ u64 v=c[o][a][b][d]; if(!v) continue;
        n[o][h1?2:1][b][d]+=v; if(b||h1) n[o][a][h2?2:1][d]+=v; if(d||h2) n[o][a][b][1]+=v; } }
    memcpy(c,n,sizeof c); }
  u64 H=0; for(int o=0;o<2;o++){ int h12=o?!e:e; for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++) if(!(a==2||b==2||(b==0&&!h12))) H+=c[o][a][b][d]; }
  return H; }
#define B 14000
static unsigned char ok[B+1];
static int ispow3(u64 v){ if(!v) return 0; while(v%3==0) v/=3; return v==1; }
int main(void){
  int xx[32],yy[32];
  for(L=0;L<=9;L++) for(u64 code=0; code<(1ULL<<(2*L)); code++){
    for(int i=0;i<L;i++){ X[i]=(code>>(2*i))&1; Y[i]=(code>>(2*i+1))&1; }
    for(int e=0;e<2;e++){ u64 m=Hof(X,Y,L,e); if(2*m>B) continue;
      for(int a=0;a<4;a++)for(int b=0;b<4;b++)for(int f=0;f<2;f++){
        // two appended letters a,b (each letter = 2 bits)
        memcpy(xx,X,sizeof(int)*L); memcpy(yy,Y,sizeof(int)*L); xx[L]=a&1; yy[L]=a>>1; xx[L+1]=b&1; yy[L+1]=b>>1;
        u64 n=Hof(xx,yy,L+2,f); if(n<=B && n>2*m && ispow3(n-2*m)) ok[n]=1;
        // two prepended
        xx[0]=a&1; yy[0]=a>>1; xx[1]=b&1; yy[1]=b>>1; for(int i=0;i<L;i++){xx[i+2]=X[i]; yy[i+2]=Y[i];}
        n=Hof(xx,yy,L+2,f); if(n<=B && n>2*m && ispow3(n-2*m)) ok[n]=1;
        // one each end
        xx[0]=a&1; yy[0]=a>>1; for(int i=0;i<L;i++){xx[i+1]=X[i]; yy[i+1]=Y[i];} xx[L+1]=b&1; yy[L+1]=b>>1;
        n=Hof(xx,yy,L+2,f); if(n<=B && n>2*m && ispow3(n-2*m)) ok[n]=1; }
      for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int f=0;f<2;f++){
        memcpy(xx,X,sizeof(int)*L); memcpy(yy,Y,sizeof(int)*L); xx[L]=a; yy[L]=b;
        u64 n=Hof(xx,yy,L+1,f); if(n<=B && n>2*m && ispow3(n-2*m)) ok[n]=1;
        xx[0]=a; yy[0]=b; for(int i=0;i<L;i++){xx[i+1]=X[i]; yy[i+1]=Y[i];}
        n=Hof(xx,yy,L+1,f); if(n<=B && n>2*m && ispow3(n-2*m)) ok[n]=1; }
    } }
  int bad=0; printf("odd n in [23,%d] (n!=21) NOT reached by the decoder step from any realizer of length<=10:",B);
  for(int n=23;n<=B;n+=2){ if(!ok[n]){ bad++; if(bad<=30) printf(" %d",n);} }
  printf("\ntotal %d\n",bad); return 0; }
