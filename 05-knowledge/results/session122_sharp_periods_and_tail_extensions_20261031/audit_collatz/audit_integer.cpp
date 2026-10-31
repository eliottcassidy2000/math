// Adapted from precision_lifts_ramsey_deficits_20261030/audit_adaptive.cpp
// at d345f993475dd7b5272ac6c3b4ab693899515a38; arithmetic independently z-major.
// Integer 5/8 slope, letter cap64; no producer optimization imported.
// Independent edge checker: z-major enumeration, repeated modular halving,
// no maximization or producer transition routine. Failures use explicit exits.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
using namespace std; using I=long long;
int main(int argc,char**argv){
 if(argc!=3&&argc!=4)return 2; int E=atoi(argv[1]),R=argc==4?atoi(argv[3]):1;if(E<1||E>10||R<1||R>E)return 2;
 vector<int> p(E+3,1);for(int j=1;j<E+3;j++)p[j]=5*p[j-1];
 int q=p[E];I D=8,N=5; const int cap=64;
 vector<int> base(E+1),width(E+1);I size=0;
 for(int h=0;h<=E;h++){base[h]=size;width[h]=4*p[max(R,h)-1];size+=1LL*(h==E?1:4*q/p[h+1])*width[h];if(size>numeric_limits<int>::max()){cerr<<"state count exceeds index type\n";return 8;}}
 auto val=[&](int x){if(!x)return E;int h=0;while(x%5==0){x/=5;++h;}return h;};
 auto index=[&](int z,int m){int h=val(z),j=(z==0?0:z/p[h]);return base[h]+(z==0?0:(j-1-j/5)*width[h])+m-1-m/5;};
 ifstream f(argv[2],ios::binary|ios::ate);if(!f||f.tellg()!=8*size){cerr<<"bad binary size\n";return 3;}f.seekg(0);
 vector<int16_t> V(size);for(I i=0;i<size;i++){unsigned char b[8];f.read((char*)b,8);uint64_t a=0;for(int j=7;j>=0;j--)a=256*a+b[j];I decoded=(I)a; if(decoded<=-30000||decoded>0){cerr<<"invalid encoded value\n";return 9;} V[i]=decoded;}
 I lo=*min_element(V.begin(),V.end()),hi=*max_element(V.begin(),V.end());
 if(lo<=-30000||hi>0){cerr<<"invalid potential range\n";return 4;}
 auto half=[&](int x){return (x+(x%2?q:0))/2;};
 I edges=0,entries=0,bounded_edges=0,explicit_edges=0;
 for(int z=0;z<q;z++){
  if(z%1000000==0)cerr<<"independent audit z="<<z<<" edges="<<edges<<"\n";
  int h=val(z),mod=p[max(R,h)];
  int targetz[65],factor[65];int zz=(3LL*z+q-2)%q,ff=3;
  for(int a=1;a<=cap;a++){zz=half(zz);ff=half(ff);targetz[a]=(zz+1)%q;factor[a]=ff;}
  for(int m=1;m<mod;m++)if(m%5){
   I src=V[index(z,m)];
   for(int a=1;a<=cap;a++){
    int z2=targetz[a],h2=val(z2),mod2=p[max(R,h2)],common=min(mod,mod2);
    int residue=(I)m*factor[a]%common;
    int reward_upper=min(h2,val(residue-1));
    if(src+D*reward_upper-N*a<=lo){I count=mod2/common;edges+=count;bounded_edges+=count;continue;}
    for(int m2=residue;m2<mod2;m2+=common){
     int reward=min(h2,val(m2-1));I rhs=src+D*reward-N*a;
     if(V[index(z2,m2)]<rhs){cerr<<"edge failure "<<z<<" "<<m<<" "<<a<<" "<<m2<<"\n";return 5;}
     ++edges;++explicit_edges;
    }
   }
  }
 }
 int inv=1;for(int a=1;a<=cap;a++){
  inv=half(inv);if(a==1)continue;int z=(q+1-2LL*inv%q)%q,h=val(z),mod=p[max(R,h)],T=2*p[max(R,h)-1],m=1;
  for(int k=0;k<T;k++){
   int m2=3LL*m*inv%mod;I charge=0;for(int j=1;j<=E;j++)charge+=k/(2*p[j-1]);
   charge+=min(h,val(m2-1));
   if(V[index(z,m2)]<D*charge-N*(k+a)){cerr<<"entry failure "<<a<<" "<<k<<"\n";return 6;}
   ++entries;m=(3LL*m+((3*m)%2?mod:0))/2%mod;
  }
 }
 bool tail=(hi+D*E-(cap+1)*N<=lo&&D*E-(cap+1)*N<=lo);
 if(!tail){cerr<<"tail failure\n";return 7;}
 cout<<"{\"E\":"<<E<<",\"R\":"<<R<<",\"states\":"<<size<<",\"edges\":"<<edges<<",\"explicit_edges\":"<<explicit_edges<<",\"bounded_edges\":"<<bounded_edges<<",\"entries\":"<<entries<<",\"D\":"<<D<<",\"N\":"<<N<<",\"min\":"<<lo<<",\"max\":"<<hi<<",\"terminal\":"<<V[index(0,1)]<<",\"all_positive_letters\":true,\"pass\":true}\n";
}
