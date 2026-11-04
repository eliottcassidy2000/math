#include <algorithm>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;using ll=long long;
bool strong(const vector<int>&r,int S){if(!S)return false;int b=S&-S;for(int rev=0;rev<2;++rev){int reach=b;while(1){int nxt=reach;for(int i=0;i<(int)r.size();++i)if(!rev&&(reach>>i&1))nxt|=r[i]&S;else if(rev&&(r[i]&reach))nxt|=(1<<i)&S;if(nxt==reach)break;reach=nxt;}if(reach!=S)return false;}return true;}
pair<ll,ll> counts(const vector<int>&r){int n=r.size(),N=1<<n;vector<ll>d(N*n);for(int i=0;i<n;++i)d[(1<<i)*n+i]=1;for(int s=1;s<N;++s)for(int i=0;i<n;++i){ll z=d[s*n+i];if(!z)continue;int rem=r[i]&~s;while(rem){int b=rem&-rem;rem-=b;d[(s|b)*n+__builtin_ctz(b)]+=z;}}ll h=0,h0=0;for(int i=0;i<n;++i){h+=d[(N-1)*n+i];h0+=d[(N-2)*n+i];}return {h,h0};}
int main(int argc,char**argv){int n=atoi(argv[1]),cap=atoi(argv[2]),K;cin>>K;ll tested=0,valid=0,kept=0,best=LLONG_MAX;vector<int>bw;pair<ll,ll>bc;for(int k=0;k<K;++k){vector<int>r(n);for(auto &x:r)cin>>x;for(int m=1;m<(1<<n)-1;++m){vector<int>t=r;for(int i=0;i<n;++i)if(m>>i&1)t[i]|=1<<n;t.push_back(((1<<n)-1)^m);++tested;if(argc<4 && !strong(t,((1<<(n+1))-1)^1))continue;++valid;auto [h,h0]=counts(t);ll delta=h-h0;if(delta<best){best=delta;bw=t;bc={h,h0};}if(delta>cap)continue;++kept;cout<<h<<' '<<h0;for(auto x:t)cout<<' '<<x;cout<<'\n';}}cerr<<"parents "<<K<<" tested "<<tested<<" valid "<<valid<<" kept "<<kept<<" min "<<best<<" h "<<bc.first<<" parent_h "<<bc.second<<" witness";for(auto x:bw)cerr<<' '<<x;cerr<<endl;}
