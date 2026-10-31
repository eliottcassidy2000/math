// Height-adaptive non-one-word price producer. A checker is mandatory.
#include <vector>
#include <deque>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <limits>
using namespace std;using ll=long long;
int main(int argc,char**argv){
 if(argc!=3&&argc!=4)return 2;int E=atoi(argv[1]),R=argc==4?atoi(argv[3]):1;if(E<1||E>10||R<1||R>E)return 3;
 vector<int>pw(E+3,1);for(int i=1;i<E+3;i++)pw[i]=5*pw[i-1];
 int q=pw[E],cap=64;ll D=8,N=5,NEG=-30000;
 vector<unsigned char>ht(q,E);for(int z=1;z<q;z++){int t=z;ht[z]=0;while(t%5==0)ht[z]++,t/=5;}
 vector<int>r(E+1),u(E+1),offset(E+2);for(int h=0;h<=E;h++){r[h]=pw[max(R,h)];u[h]=4*r[h]/5;int nz=h==E?1:4*(q/pw[h+1]);ll next=offset[h]+1LL*nz*u[h];if(next>numeric_limits<int>::max()){cerr<<"state count exceeds index type\n";return 7;}offset[h+1]=next;}
 auto rank=[](int x){return x-x/5-1;};auto unit=[](int j){return j+1+j/4;};
 auto id=[&](int z,int m){int h=ht[z];return offset[h]+(h==E?0:rank(z/pw[h]))*u[h]+rank(m);};
 int size=offset[E+1];vector<int16_t>V(size,NEG);vector<char>in(size);deque<int>todo;
 auto relax=[&](int t,ll value){if(value>V[t]){V[t]=value;if(value>0){cerr<<"POSITIVE abstract excess\n";exit(4);}if(!in[t])todo.push_back(t),in[t]=1;}};
 vector<int>inv(cap+1);inv[0]=1;for(int a=1;a<=cap;a++)inv[a]=1LL*inv[a-1]*((q+1)/2)%q;
 // All initial-one prefixes reduce to their least phase representative at this target precision.
 for(int a=2;a<=cap;a++){
  int z=(1-2LL*inv[a]%q+q)%q,H=ht[z],rr=r[H],T=2*pw[max(R,H)-1],m=1;
  for(int k=0;k<T;k++){
   int mm=3LL*m*inv[a]%rr;ll reward=0;for(int e=1;e<=E;e++)reward+=k/(2*pw[e-1]);
   relax(id(z,mm),D*(reward+min<int>(ht[mm-1],H))-N*(k+a));
   m=3LL*m*((rr+1)/2)%rr;
  }
 }
 ll pops=0;while(!todo.empty()){
  int s=todo.front();todo.pop_front();in[s]=0;++pops; if(pops%10000000==0) cerr<<"pops="<<pops<<" todo="<<todo.size()<<endl;
  int h=upper_bound(offset.begin(),offset.end(),s)-offset.begin()-1;
  int rem=s-offset[h],m=unit(rem%u[h]),z=h==E?0:pw[h]*unit(rem/u[h]);ll value=V[s];
  for(int a=1;a<=cap;a++){
   int zz=(3LL*z*inv[a]+1-2LL*inv[a]%q+q)%q,H=ht[zz],rr=r[H],small=min(r[h],rr);
   int first=3LL*m*inv[a]%small;
   for(int mm=first;mm<rr;mm+=small)relax(id(zz,mm),value+D*min<int>(ht[mm-1],H)-N*a);
  }
 }
 ll mn=0,mx=NEG;for(ll x:V){if(x==NEG){cerr<<"unreached\n";return 5;}mn=min(mn,x);mx=max(mx,x);}
 cout<<"E="<<E<<" R="<<R<<" states="<<size<<" terminal="<<V[id(0,1)]<<" D="<<D<<" min="<<mn<<" max="<<mx<<" pops="<<pops<<" tail="<<(mx+D*E-(cap+1)*N<=mn&&D*E-(cap+1)*N<=mn)<<endl;
 ofstream f(argv[2],ios::binary);for(ll x:V){uint64_t v=(uint64_t)x;char b[8];for(int j=0;j<8;j++)b[j]=(v>>(8*j))&255;f.write(b,8);}return f?0:6;
}
