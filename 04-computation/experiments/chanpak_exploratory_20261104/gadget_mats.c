// Gadget matrices M_{G,u,z} for the near-join T = T1 => G with back arc u -> v (v = port of T1), new port z in G:
//   (h, e_z, B_z)(T) = M (h1, e_v, B_v)(T1),
// M = [[h_G, -s_u, E_u], [e_z, -P_{u->z}, E_{u,z}], [B_z, -J1, J2 + K]] (formulas in the S78 note).
// Usage: gm MMAX [dump]   enumerates every labelled tournament G of order 1..MMAX and all (u,z); prints matrices if dump.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
#define MX 8
static int m; static uint32_t out[MX], inn[MX];
static ll hS[1<<MX], eS[1<<MX][MX], sS[1<<MX][MX], P[1<<MX][MX][MX]; // P[S][a][b]: paths covering S from a to b
static void tables(void){
  int full=(1<<m)-1;
  memset(P,0,sizeof(ll)*(size_t)(full+1)*MX*MX);
  for(int a=0;a<m;a++) P[1<<a][a][a]=1;
  for(int S=1;S<=full;S++) for(int a=0;a<m;a++) if(S>>a&1) for(int b=0;b<m;b++){ ll c=P[S][a][b]; if(!c) continue;
      uint32_t nb=out[b]&~(uint32_t)S; while(nb){ int w=__builtin_ctz(nb); nb&=nb-1; P[S|(1<<w)][a][w]+=c; } }
  hS[0]=1;
  for(int S=1;S<=full;S++){ ll t=0; for(int a=0;a<m;a++){ eS[S][a]=0; sS[S][a]=0; }
    for(int a=0;a<m;a++) for(int b=0;b<m;b++){ ll c=P[S][a][b]; t+=c; sS[S][a]+=c; eS[S][b]+=c; }
    hS[S]=t; }
  for(int a=0;a<m;a++){ eS[0][a]=0; sS[0][a]=0; }
}
static ll Eu_of(int S,int u){ // end response of u in G[S]: sum_{B subset S, u in B} e_u(B) h(S-B)
  ll t=0; for(int B=S; ; B=(B-1)&S){ if(B>>u&1) t+=eS[B][u]*hS[S&~B]; if(B==0) break; } return t; }
static void matrix(int u,int z,ll M[3][3]){
  int full=(1<<m)-1;
  ll hG=hS[full], su=sS[full][u], Eu=Eu_of(full,u);
  ll ez=eS[full][z], Puz=P[full][u][z], Euz=0;
  for(int B1=0;B1<=full;B1++) if((B1>>u&1) && !(B1>>z&1)) Euz+=eS[B1][u]*eS[full&~B1][z];
  ll Bz=0; for(int S=0;S<=full;S++) if(!(S>>z&1)) Bz+=hS[S]*sS[full&~S][z];
  ll J1=0,J2=0; for(int S2=0;S2<=full;S2++) if((S2>>u&1) && !(S2>>z&1)){ J1+=sS[S2][u]*sS[full&~S2][z]; J2+=Eu_of(S2,u)*sS[full&~S2][z]; }
  ll K=0; int rest=full&~(1<<z)&~(1<<u);
  for(int S2=rest; ; S2=(S2-1)&rest){ int R=full&~S2; // Y1 subset R containing z,u
     int Rr=R&~(1<<z)&~(1<<u);
     for(int X=Rr; ; X=(X-1)&Rr){ int Y1=X|(1<<z)|(1<<u); ll p=P[Y1][z][u]; if(p) K+=hS[S2]*p*hS[R&~Y1]; if(X==0) break; }
     if(S2==0) break; }
  M[0][0]=hG; M[0][1]=-su; M[0][2]=Eu;
  M[1][0]=ez; M[1][1]=-Puz; M[1][2]=Euz;
  M[2][0]=Bz; M[2][1]=-J1; M[2][2]=J2+K;
}
int main(int argc,char**argv){
  int MMAX=argc>1?atoi(argv[1]):4, dump=argc>2?atoi(argv[2]):0;
  for(m=1;m<=MMAX;m++){
    int na=m*(m-1)/2;
    for(long code=0; code<(1L<<na); code++){
      memset(out,0,sizeof out); int k=0;
      for(int i=0;i<m;i++) for(int j=i+1;j<m;j++,k++){ if(code>>k&1) out[i]|=1u<<j; else out[j]|=1u<<i; }
      tables();
      for(int u=0;u<m;u++) for(int z=0;z<m;z++){ ll M[3][3]; matrix(u,z,M);
        if(dump) printf("%d %ld %d %d %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",m,code,u,z,M[0][0],M[0][1],M[0][2],M[1][0],M[1][1],M[1][2],M[2][0],M[2][1],M[2][2]); }
    }
  }
  return 0;
}
