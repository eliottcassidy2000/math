#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
using namespace std;
using I=unsigned __int128;
string text(I x){if(!x)return "0";string s;while(x){s.push_back(char('0'+x%10));x/=10;}return string(s.rbegin(),s.rend());}
int main(){unsigned tests;if(!(cin>>tests))return 1;while(tests--){unsigned n,lo,hi;cin>>n>>lo>>hi;if(n>11||lo<1||hi<1||lo>8||hi>8)return 2;vector<unsigned>A(n+2);for(auto&x:A)cin>>x;unsigned N=1u<<n,ks=lo+1,ls=hi+1;vector<I> dp(size_t(N)*ks*ls*(n+2));auto at=[&](unsigned mask,unsigned k,unsigned l,unsigned v)->I&{return dp[(((size_t(mask)*ks+k)*ls+l)*(n+2))+v];};for(unsigned v=0;v<n;v++)at(1u<<v,0,0,v)=1;at(0,1,0,n)=1;at(0,0,1,n+1)=1;for(unsigned mask=0;mask<N;mask++)for(unsigned k=0;k<=lo;k++)for(unsigned l=0;l<=hi;l++)for(unsigned u=0;u<n+2;u++){I val=at(mask,k,l,u);if(!val)continue;for(unsigned v=0;v<n+2;v++){if(!(A[u]>>v&1))continue;if(v<n){if(!(mask>>v&1))at(mask|(1u<<v),k,l,v)+=val;}else if(v==n){if(k<lo)at(mask,k+1,l,v)+=val;}else if(l<hi)at(mask,k,l+1,v)+=val;}}vector<vector<I>>S(1+max(lo,hi),vector<I>(1+max(lo,hi)));S[0][0]=1;for(unsigned t=1;t<S.size();t++)for(unsigned h=1;h<=t;h++)S[t][h]=h*(S[t-1][h]+S[t-1][h-1]);I ans=0;for(unsigned k=1;k<=lo;k++)for(unsigned l=1;l<=hi;l++)for(unsigned u=0;u<n+2;u++)ans+=at(N-1,k,l,u)*S[lo][k]*S[hi][l];cout<<text(ans)<<'\n';}}
