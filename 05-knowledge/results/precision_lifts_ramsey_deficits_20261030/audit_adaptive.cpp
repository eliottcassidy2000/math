// Independent edge checker: z-major enumeration, repeated modular halving,
// no maximization or producer transition routine. Failures use explicit exits.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std; using I=long long;
int main(int argc,char**argv){
 if(argc!=3)return 2; int E=atoi(argv[1]);if(E<1||E>10)return 2;
 vector<int> p(E+3,1);for(int j=1;j<E+3;j++)p[j]=5*p[j-1];
 int q=p[E];I D=200LL*q,N=125LL*q-5;
 vector<int> base(E+1),width(E+1);I size=0;
 for(int h=0;h<=E;h++){base[h]=size;width[h]=4*p[max(1,h)-1];size+=(h==E?1:4*q/p[h+1])*width[h];}
 auto val=[&](int x){if(!x)return E;int h=0;while(x%5==0){x/=5;++h;}return h;};
 auto index=[&](int z,int m){int h=val(z),j=(z==0?0:z/p[h]);return base[h]+(z==0?0:(j-1-j/5)*width[h])+m-1-m/5;};
 ifstream f(argv[2],ios::binary|ios::ate);if(!f||f.tellg()!=8*size){cerr<<"bad binary size\n";return 3;}f.seekg(0);
 vector<I> V(size);for(I i=0;i<size;i++){unsigned char b[8];f.read((char*)b,8);uint64_t a=0;for(int j=7;j>=0;j--)a=256*a+b[j];V[i]=(I)a;}
 I lo=*min_element(V.begin(),V.end()),hi=*max_element(V.begin(),V.end());
 if(lo<=-(1LL<<59)||hi>0){cerr<<"invalid potential range\n";return 4;}
 auto half=[&](int x){return (x+(x%2?q:0))/2;};
 I edges=0,entries=0;
 for(int z=0;z<q;z++){
  int h=val(z),mod=p[max(1,h)];
  int targetz[65],factor[65];int zz=(3LL*z+q-2)%q,ff=3;
  for(int a=1;a<=64;a++){zz=half(zz);ff=half(ff);targetz[a]=(zz+1)%q;factor[a]=ff;}
  for(int m=1;m<mod;m++)if(m%5){
   I src=V[index(z,m)];
   for(int a=1;a<=64;a++){
    int z2=targetz[a],h2=val(z2),mod2=p[max(1,h2)],common=min(mod,mod2);
    int residue=(I)m*factor[a]%common;
    for(int m2=residue;m2<mod2;m2+=common){
     int reward=min(h2,val(m2-1));I rhs=src+D*reward-N*a;
     if(V[index(z2,m2)]<rhs){cerr<<"edge failure "<<z<<" "<<m<<" "<<a<<" "<<m2<<"\n";return 5;}
     ++edges;
    }
   }
  }
 }
 int inv=1;for(int a=1;a<=64;a++){
  inv=half(inv);if(a==1)continue;int z=(q+1-2LL*inv%q)%q,h=val(z),mod=p[max(1,h)],T=2*p[max(1,h)-1],m=1;
  for(int k=0;k<T;k++){
   int m2=3LL*m*inv%mod;I charge=0;for(int j=1;j<=E;j++)charge+=k/(2*p[j-1]);
   charge+=min(h,val(m2-1));
   if(V[index(z,m2)]<D*charge-N*(k+a)){cerr<<"entry failure "<<a<<" "<<k<<"\n";return 6;}
   ++entries;m=(3LL*m+((3*m)%2?mod:0))/2%mod;
  }
 }
 bool tail=(hi+D*E-65*N<=lo&&D*E-65*N<=lo);
 if(!tail){cerr<<"tail failure\n";return 7;}
 cout<<"{\"E\":"<<E<<",\"states\":"<<size<<",\"edges\":"<<edges<<",\"entries\":"<<entries<<",\"D\":"<<D<<",\"N\":"<<N<<",\"min\":"<<lo<<",\"max\":"<<hi<<",\"terminal\":"<<V[index(0,1)]<<",\"all_positive_letters\":true,\"pass\":true}\n";
}
