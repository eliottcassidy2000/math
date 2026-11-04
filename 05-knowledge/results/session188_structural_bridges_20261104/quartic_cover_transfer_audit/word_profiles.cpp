#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
int main(){unsigned tests;if(!(cin>>tests))return 1;while(tests--){unsigned n;cin>>n;if(n<1||n>6)return 2;vector<unsigned>T(n),w(n);for(auto&r:T)cin>>r;iota(w.begin(),w.end(),0);uint64_t C[8][8]={};for(unsigned i=0;i<=n;i++){C[i][0]=1;for(unsigned j=1;j<=i;j++)C[i][j]=C[i-1][j-1]+C[i-1][j];}vector<uint64_t>q(n+1);do{unsigned bad=0;for(unsigned j=0;j+1<n;j++)bad+=!(T[w[j]]>>w[j+1]&1);for(unsigned t=0;t<=n-1-bad;t++)q[1+bad+t]+=C[n-1-bad][t];}while(next_permutation(w.begin(),w.end()));for(auto v:q)cout<<v<<' ';cout<<'\n';}}
