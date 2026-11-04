#include <cstdint>
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;using I=__int128_t;
string dec(I v){if(v==0)return "0";string s;while(v){s+=char('0'+v%10);v/=10;}reverse(s.begin(),s.end());return s;}
int main(){int count;cin>>count;for(int g=0;g<count;g++){int k;cin>>k;vector<int>r(k),cap(k),stride(k);for(int&x:r)cin>>x;for(int&x:cap)cin>>x;int N=1;for(int i=0;i<k;i++){stride[i]=N;N*=cap[i]+1;}vector<I>d((size_t)N*k);for(int i=0;i<k;i++)d[(size_t)stride[i]*k+i]=1;I ans=0;for(int state=1;state<N;state++){int avail=0,twos=0;bool complete=true;for(int i=0;i<k;i++){int used=state/stride[i]%(cap[i]+1);if(used<cap[i])avail|=1<<i;if(!used)complete=false;if(cap[i]==2&&used==2)twos++;}for(int last=0;last<k;last++){I val=d[(size_t)state*k+last];if(!val)continue;if(complete)ans+=val*(I(1)<<twos);int next=avail&r[last];while(next){int v=__builtin_ctz((unsigned)next);next&=next-1;d[(size_t)(state+stride[v])*k+v]+=val;}}}cout<<dec(ans)<<'\n';if(g%10==0)cerr<<"query "<<g<<" done\n";}}
