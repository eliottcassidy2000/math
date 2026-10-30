// F_2 epoch envelope: exact first two layers, delayed high-layer payment.
// Producer only. Certificate must be checked independently.
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdint>
using namespace std;using ll=long long;
int main(int argc,char**argv){
 if(argc!=3)return 2;int E=atoi(argv[1]),q=1,cap=64;if(E<3||E>9)return 3;
 for(int i=0;i<E;i++)q*=5;ll scale=25LL*q,D=8*scale,N=5*scale-5,NEG=-(1LL<<60);
 vector<int>h(q,E),pw(E+1,1),iq(cap+1),ir(cap+1),decode,ordinary(25,-1),pending(25,-1);
 for(int m=1;m<25;m++)if(m%5)ordinary[m]=decode.size(),decode.push_back(m);
 for(int m=2;m<25;m++)if(m%5==1||m%5==4)pending[m]=decode.size(),decode.push_back(m);
 int modes=decode.size();for(int i=1;i<=E;i++)pw[i]=5*pw[i-1];
 for(int z=1;z<q;z++){int x=z;h[z]=0;while(x%5==0)h[z]++,x/=5;}
 iq[0]=ir[0]=1;for(int a=1;a<=cap;a++){iq[a]=(ll)iq[a-1]*((q+1)/2)%q;ir[a]=(ll)ir[a-1]*13%25;}
 auto step=[&](int s,int a){int mode=s/q,z=s%q,m=3LL*decode[mode]*ir[a]%25,nz=(3LL*z*iq[a]+1-2LL*iq[a]%q+q)%q,H=h[nz];
   bool active=m%5==1||m%5==4;int j=ordinary[m];ll gross=D*((H>=1&&m%5==1)+(H>=2&&m==1));
   if(a>1){if(H>2&&active){if(m==1)gross+=D*(H-2);else j=pending[m];}}
   else if(mode>=20){if(m==1)gross+=D*(H-2);else j=pending[m];}
   else if(H>2&&active)gross+=scale/5-5*scale/pw[H];
   return pair<int,ll>(j*q+nz,gross-N*a);
 };
 size_t size=(size_t)modes*q;vector<ll>V(size,NEG);vector<char>in(size);deque<int>todo;
 int m=1;do{for(int a=2;a<=cap;a++){auto[t,c]=step(ordinary[m]*q,a);if(c>V[t]){V[t]=c;if(!in[t])todo.push_back(t),in[t]=1;}}m=3*m*13%25;}while(m!=1);
 ll pops=0;while(!todo.empty()){int s=todo.front();todo.pop_front();in[s]=0;pops++;for(int a=1;a<=cap;a++){auto[t,c]=step(s,a);c+=V[s];if(c>V[t]){V[t]=c;if(c>0){cerr<<"positive envelope excess\n";return 4;}if(!in[t])todo.push_back(t),in[t]=1;}}}
 ll mn=0,mx=NEG;for(ll v:V)if(v!=NEG)mn=min(mn,v),mx=max(mx,v);
 cout<<"E="<<E<<" terminal="<<V[0]<<" D="<<D<<" min="<<mn<<" max="<<mx<<" pops="<<pops<<" tail="<<(mx+D*E-65*N<=mn&&D*E-65*N<=mn)<<endl;
 ofstream out(argv[2],ios::binary);for(ll v:V){uint64_t x=(uint64_t)v;char b[8];for(int j=0;j<8;j++)b[j]=(x>>(8*j))&255;out.write(b,8);}return out?0:5;
}
