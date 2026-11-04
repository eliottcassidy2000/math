#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using W=__int128_t;
std::string dec(W x){if(!x)return "0";std::string s;while(x){s+=char('0'+x%10);x/=10;}std::reverse(s.begin(),s.end());return s;}
struct Counter{int n,states=1;std::vector<unsigned>R,available;std::vector<int>cap,stride;std::vector<unsigned char>complete;std::vector<W>memo;Counter(std::vector<unsigned>r,std::vector<int>k):n(r.size()),R(r),cap(k){stride.resize(n);for(int i=0;i<n;++i){stride[i]=states;states*=cap[i]+1;}available.resize(states);complete.resize(states);for(int s=1;s<states;++s){bool all=true;for(int i=0;i<n;++i){int used=s/stride[i]%(cap[i]+1);if(!used)all=false;if(used<cap[i])available[s]|=1u<<i;}complete[s]=all;}memo.assign(size_t(states)*n,-1);}W suffix(int state,int last){W&ans=memo[size_t(state)*n+last];if(ans>=0)return ans;ans=complete[state];unsigned todo=R[last]&available[state];while(todo){int v=__builtin_ctz(todo);todo&=todo-1;int used=state/stride[v]%(cap[v]+1);ans+=(used==1&&cap[v]==2?2:1)*suffix(state+stride[v],v);}return ans;}W total(){W ans=0;for(int i=0;i<n;++i)ans+=suffix(stride[i],i);return ans;}};
int main(){int cases;std::cin>>cases;while(cases--){int n;std::cin>>n;std::vector<unsigned>R(n);std::vector<int>cap(n);for(auto&r:R)std::cin>>r;for(auto&k:cap)std::cin>>k;std::cout<<dec(Counter(R,cap).total())<<'\n';}}
