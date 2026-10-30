#include <vector>
#include <deque>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <fstream>
using namespace std; using ll=long long;
int main(int argc,char**argv){int E=argc>1?atoi(argv[1]):4,cap=argc>2?atoi(argv[2]):64;if(E<1||E>9||cap<2||cap>4096){cerr<<"supported inputs: 1<=E<=9, 2<=cap<=4096\n";return 4;}int q=1;for(int i=0;i<E;i++)q*=5;ll scale=8LL*q,coef=5LL*q-5;vector<int> h(q,E),pw(E+1,1);for(int i=1;i<=E;i++)pw[i]=5*pw[i-1];for(int z=1;z<q;z++){int t=z;h[z]=0;while(t%5==0){h[z]++;t/=5;}}
vector<int> inv(cap+1),k(cap+1);inv[0]=1;for(int a=1;a<=cap;a++){inv[a]=(ll)inv[a-1]*((q+1)/2)%q;int z=(1-2LL*inv[a]%q+q)%q;k[a]=h[z];}
auto next=[&](int s,int a){int e=s/q,z=s%q;return ((e+(a+1)%2)%2)*q+(3LL*z*inv[a]+1-2LL*inv[a]%q+q)%q;};
auto delta=[&](int s,int a,int t)->ll{if(a==1)return s/q==0?-(5LL*q/pw[h[s%q]]-5):-coef;return (t/q==0?scale*max(h[t%q],k[a]):0)-coef*a;};
const ll neg=-(1LL<<60);vector<ll> V(2*q,neg);vector<int> pred(2*q,-1),letter(2*q);vector<bool> in(2*q,false);deque<int> todo;for(int a=2;a<=cap;a++){int t=next(0,a);ll v=delta(0,a,t);if(v>V[t]){V[t]=v;pred[t]=-1;letter[t]=a;if(!in[t])todo.push_back(t),in[t]=true;}}
long long pops=0,relax=0;while(!todo.empty()){int s=todo.front();todo.pop_front();in[s]=false;pops++;for(int a=1;a<=cap;a++){int t=next(s,a);ll val=V[s]+delta(s,a,t);if(val>V[t]){V[t]=val;pred[t]=s;letter[t]=a;relax++;if(val>0){cerr<<"FAIL "<<E<<" "<<val<<" state "<<t<<" word backwards:";int x=t;for(int j=0;x>=0&&j<1000;j++){cerr<<" "<<letter[x];x=pred[x];}cerr<<endl;return 1;}if(!in[t])todo.push_back(t),in[t]=true;}}}
ll mn=*min_element(V.begin(),V.end()),mx=*max_element(V.begin(),V.end()); bool tail=mx+scale*E-coef*(cap+1)<=mn;cout<<"E="<<E<<" q="<<q<<" states="<<2*q<<" min="<<mn<<" max="<<mx<<" pops="<<pops<<" relax="<<relax<<" tail="<<tail<<" entry_tail="<<(scale*E-coef*(cap+1)<=mn)<<endl;
if(argc>3){ofstream out(argv[3],ios::binary);out.write(reinterpret_cast<char*>(V.data()),V.size()*sizeof(ll));}
for(int s=0;s<2*q;s++)for(int a=1;a<=cap;a++){int t=next(s,a);if(V[t]<V[s]+delta(s,a,t))return 2;}return tail?0:3;
}
