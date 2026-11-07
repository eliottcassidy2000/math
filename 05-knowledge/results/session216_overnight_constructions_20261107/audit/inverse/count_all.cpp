#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
using namespace std;
int main(){unsigned long long want;string s;int count=0,maxn=0;while(cin>>want>>s){if(s=="-")s="";int n=s.size()+2;if(n>12)return 10;vector<unsigned> in(n);in[1]|=1;for(int v=2;v<n;v++){int mask=s[v-2]-'0';if(mask<0||mask>3)return 11;for(int a=0;a<2;a++){if(mask>>a&1)in[v]|=1u<<a;else in[a]|=1u<<v;}for(int u=2;u<v;u++)in[v]|=1u<<u;}
// Pull recurrence: count paths ending at v from predecessor-ending paths.
vector<uint64_t> f((1<<n)*n);for(unsigned set=1;set<(1u<<n);set++)for(int v=0;v<n;v++)if(set>>v&1){unsigned rest=set^(1u<<v);if(!rest){f[set*n+v]=1;continue;}unsigned p=rest&in[v];while(p){int u=__builtin_ctz(p);p&=p-1;f[set*n+v]+=f[rest*n+u];}}
uint64_t got=0;for(int v=0;v<n;v++)got+=f[((1<<n)-1)*n+v];if(got!=want){cerr<<want<<" != "<<got<<endl;return 12;}count++;maxn=max(maxn,n);}if(count!=4998)return 13;cout<<"{\"checked\":"<<count<<",\"max_vertices\":"<<maxn<<",\"overflow_bound\":\"12! < 2^64\"}"<<endl;}
