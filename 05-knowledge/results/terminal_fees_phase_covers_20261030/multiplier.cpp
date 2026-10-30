// Session113 epoch envelope with an additional exact multiplier residue.
// This is a producer, not a proof checker.
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdint>
using namespace std;using ll=long long;
int main(int argc,char**argv){
 if(argc<4)return 2;int E=atoi(argv[1]),S=atoi(argv[2]),q=1,r=1,cap=64;
 if(E<1||E>9||S<1||S>4)return 3;
 for(int i=0;i<E;i++)q*=5;for(int i=0;i<S;i++)r*=5;
 ll scale=25LL*q,D=8*scale,N=5*scale-5,NEG=-(1LL<<60);
 int units=4*r/5,modes=r;size_t size=(size_t)modes*q;
 vector<int>h(q,E),pw(E+1,1),iq(cap+1),ir(cap+1),decode(units);
 for(int i=1;i<=E;i++)pw[i]=5*pw[i-1];
 for(int z=1;z<q;z++){int x=z;h[z]=0;while(x%5==0)h[z]++,x/=5;}
 for(int m=1;m<r;m++)if(m%5)decode[m-m/5-1]=m;
 iq[0]=ir[0]=1;for(int a=1;a<=cap;a++){iq[a]=(ll)iq[a-1]*((q+1)/2)%q;ir[a]=(ll)ir[a-1]*((r+1)/2)%r;}
 auto next=[&](int s,int a){int mode=s/q,u=mode>=units?4+5*(mode-units):decode[mode],z=s%q;int nu=3LL*u*ir[a]%r,nz=(3LL*z*iq[a]+1-2LL*iq[a]%q+q)%q;int j=nu-nu/5-1;if(a>1&&nu%5==4&&h[nz])j=units+(nu-4)/5;return j*q+nz;};
 auto cost=[&](int s,int a,int t){int k=s/q,j=t/q,H=h[t%q],nu=j>=units?4+5*(j-units):decode[j];if(a>1)return (nu%5==1?D*H:0LL)-N*a;if(k>=units)return D*H-N;if(H&&(nu%5==1||nu%5==4))return (nu%5==1?D:0LL)+scale-5*scale/pw[H]-N;return -N;};
 vector<ll>V(size,NEG);vector<char>in(size);deque<int>todo;
 vector<int>origins;int u=1;do{origins.push_back(u);u=3LL*u*ir[1]%r;}while(u!=1);
 for(int u:origins)for(int a=2;a<=cap;a++){int src=(u-u/5-1)*q,t=next(src,a);ll c=cost(src,a,t);if(c>V[t]){V[t]=c;if(!in[t])todo.push_back(t),in[t]=1;}}
 ll pops=0;while(!todo.empty()){int s=todo.front();todo.pop_front();in[s]=0;pops++;for(int a=1;a<=cap;a++){int t=next(s,a);ll c=V[s]+cost(s,a,t);if(c>V[t]){V[t]=c;if(c>0){cerr<<"positive envelope excess\n";return 4;}if(!in[t])todo.push_back(t),in[t]=1;}}}
 ll mn=0,mx=NEG;for(ll v:V)if(v!=NEG)mn=min(mn,v),mx=max(mx,v);
 cout<<"E="<<E<<" S="<<S<<" terminal="<<V[0]<<" D="<<D<<" min="<<mn<<" max="<<mx<<" pops="<<pops<<" tail="<<(mx+D*E-65*N<=mn&&D*E-65*N<=mn)<<endl;
 ofstream out(argv[3],ios::binary);for(ll v:V){uint64_t x=(uint64_t)v;char b[8];for(int j=0;j<8;j++)b[j]=(x>>(8*j))&255;out.write(b,8);}return out?0:5;
}
