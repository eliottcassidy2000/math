// Low-block two-hub words: x = 2-block (s, s^1 at cut p), y = <=3 blocks. Record which 1+2^a3^b are hit.
#include <stdio.h>
#include <string.h>
typedef unsigned __int128 u128;
static int X[64],Y[64],E,L;
static u128 H2(void){
  u128 c[2][3][3][2]; memset(c,0,sizeof c); c[0][0][0][0]=1; c[1][0][0][0]=1;
  for(int i=0;i<L;i++){ u128 n[2][3][3][2]; memset(n,0,sizeof n);
    for(int o=0;o<2;o++){ int h1=o?Y[i]:X[i], h2=o?X[i]:Y[i];
      for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ u128 v=c[o][a][b][d]; if(!v) continue;
        n[o][h1?2:1][b][d]+=v; if(b||h1) n[o][a][h2?2:1][d]+=v; if(d||h2) n[o][a][b][1]+=v; } }
    memcpy(c,n,sizeof c); }
  u128 H=0; for(int o=0;o<2;o++){ int h12=o?!E:E;
    for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++) if(!(a==2||b==2||(b==0&&!h12))) H+=c[o][a][b][d]; }
  return H; }
static int hit[40][26]; static int ex[40][26][6];
int main(void){
  for(L=1;L<=26;L++)
   for(int sx=0;sx<2;sx++)for(int p=0;p<=L;p++)
    for(int sy=0;sy<2;sy++)for(int q=0;q<=L;q++)for(int r=q;r<=L;r++)
     for(E=0;E<2;E++){
       for(int i=0;i<L;i++){ X[i]= i<p? sx: !sx; Y[i]= i<q? sy: (i<r? !sy: sy); }
       u128 H=H2(); if(H<2) continue; u128 t=H-1; int a=0,b=0; while(t%2==0){t/=2;a++;} while(t%3==0){t/=3;b++;}
       if(t==1 && a<40 && b<26 && !hit[a][b]){ hit[a][b]=1; int*e=ex[a][b]; e[0]=L;e[1]=sx;e[2]=p;e[3]=sy;e[4]=q;e[5]=r*2+E; }
     }
  int miss=0,tot=0,mmiss=0,mtot=0;
  for(int a=3;a<40;a++)for(int b=0;b<26;b++){ double v=1; for(int i=0;i<a;i++)v*=2; for(int i=0;i<b;i++)v*=3; if(v>3e11) continue;
      tot++; if(!hit[a][b]){ miss++; if(miss<=12) printf("micro miss 1+2^%d*3^%d\n",a,b);} 
      if(a>=2&&b>=2){ mtot++; if(!hit[a][b]) mmiss++; } }
  printf("microcosm 1+2^a3^b (a>=3, value<=3e11): %d of %d realized by low-block two-hub words (len<=26)\n",tot-miss,tot);
  printf("macrocosm 1+36*2^k3^l: %d of %d realized\n",mtot-mmiss,mtot);
  printf("1+2^a with a<3: %d %d %d ; 1297=1+2^4 3^4 hit? %d (L=%d)\n",hit[0][0],hit[1][0],hit[2][0],hit[4][4],ex[4][4][0]);
  return 0; }
