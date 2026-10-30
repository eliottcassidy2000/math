// Adapted from Session113 uniform_content/repaired.cpp. State precision E, slope precision E+2.
#include <vector>
#include <deque>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdint>
using namespace std;using ll=long long;
int main(int argc,char**argv){int E=argc>1?atoi(argv[1]):3,cap=64,q=1;if(E<1||E>9){cerr<<"supported precision is1..9"<<endl;return 4;}for(int i=0;i<E;i++)q*=5;int R=E+2,scale=q*25;ll D=8LL*scale,N=5LL*scale-5;vector<int>h(q,E),pw(E+1,1),inv(cap+1),iv5(cap+1);for(int i=1;i<=E;i++)pw[i]=5*pw[i-1];for(int z=1;z<q;z++){int t=z;h[z]=0;while(t%5==0)h[z]++,t/=5;}inv[0]=iv5[0]=1;for(int a=1;a<=cap;a++){inv[a]=(ll)inv[a-1]*((q+1)/2)%q;iv5[a]=iv5[a-1]*3%5;}
auto next=[&](int s,int a){int k=s/q,u=k==4?4:k+1,z=s%q;int v=u*3*iv5[a]%5,y=(3LL*z*inv[a]+1-2LL*inv[a]%q+q)%q;int j=v-1;if(a>1&&v==4&&h[y]>0)j=4;return j*q+y;};
auto cost=[&](int s,int a,int t)->ll{int k=s/q,j=t/q,H=h[t%q];if(a>1)return (j==0?D*H:0)-N*a;if(k==4)return D*H-N;if(H>0&&(j==0||j==3))return (j==0?D:0)+scale-5LL*scale/pw[H]-N;return -N;};
const ll NEG=-(1LL<<60);vector<ll>V(5*q,NEG);vector<int>pred(5*q,-1),letter(5*q),origin(5*q);vector<char>in(5*q);deque<int>todo;
for(int o:{0,3*q})for(int a=2;a<=cap;a++){int t=next(o,a);ll c=cost(o,a,t);if(c>V[t]){V[t]=c;pred[t]=-1;origin[t]=o/q+1;letter[t]=a;if(!in[t])todo.push_back(t),in[t]=1;}}
ll pops=0;while(!todo.empty()){int s=todo.front();todo.pop_front();in[s]=0;pops++;for(int a=1;a<=cap;a++){int t=next(s,a);ll c=V[s]+cost(s,a,t);if(c>V[t]){V[t]=c;pred[t]=s;letter[t]=a;origin[t]=origin[s];if(c>0){cerr<<"FAIL E="<<E<<" score="<<c<<" D="<<D<<" origin="<<origin[t]<<" REVERSE";int x=t;for(int j=0;x>=0&&j<100000;j++){cerr<<" "<<letter[x];x=pred[x];}cerr<<endl;return 1;}if(!in[t])todo.push_back(t),in[t]=1;}}}
// Pending states with h=0 are unreachable. Ignore them in certificate extrema.
ll mn=0,mx=NEG;int unreachable=0;for(int s=0;s<5*q;s++){if(V[s]==NEG){unreachable++;continue;}mn=min(mn,V[s]);mx=max(mx,V[s]);}
bool tail=mx+D*E-N*(cap+1)<=mn&&D*E-N*(cap+1)<=mn;
cout<<"terminal="<<V[0]<<" D="<<D<<" slopeE="<<R<<" E="<<E<<" min="<<mn<<" max="<<mx<<" pops="<<pops<<" unreachable="<<unreachable<<" tail="<<tail<<endl;
if(argc>2){ofstream f(argv[2],ios::binary);for(ll value:V){uint64_t bits=static_cast<uint64_t>(value);char bytes[8];for(int j=0;j<8;j++)bytes[j]=static_cast<char>((bits>>(8*j))&255);f.write(bytes,8);}if(!f){cerr<<"binary output failed"<<endl;return 5;}}return tail?0:2;}
