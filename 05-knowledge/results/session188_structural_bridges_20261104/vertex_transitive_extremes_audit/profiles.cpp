#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
int main(){unsigned tests;if(!(cin>>tests))return 1;while(tests--){unsigned n;cin>>n;if(n==0||n>9)return 2;vector<unsigned>T(n),w(n);for(auto&x:T)cin>>x;iota(w.begin(),w.end(),0);uint64_t C[10][10]={};for(unsigned i=0;i<=n;i++){C[i][0]=1;for(unsigned j=1;j<=i;j++)C[i][j]=C[i-1][j-1]+C[i-1][j];}vector<uint64_t>q(n+1),E(n*(n+1)),S(n*(n+1));do{vector<unsigned>bad(n);unsigned b=0;for(unsigned j=0;j+1<n;j++){bad[j]=!(T[w[j]]>>w[j+1]&1);b+=bad[j];}for(unsigned z=0;z<=n-1-b;z++){unsigned h=1+b+z;q[h]+=C[n-1-b][z];S[w[0]*(n+1)+h]+=C[n-1-b][z];}unsigned suffix=b;bool prefix=true;for(unsigned t=0;t<n&&prefix;t++){if(t+1<n)suffix-=bad[t];unsigned mandatory=(t+1<n?1:0)+suffix;unsigned free=(t+1<n?n-t-2:0)-suffix;for(unsigned z=0;z<=free;z++)E[w[t]*(n+1)+1+mandatory+z]+=C[free][z];if(t+1<n&&bad[t])prefix=false;}}while(next_permutation(w.begin(),w.end()));for(auto x:q)cout<<x<<' ';for(auto x:E)cout<<x<<' ';for(auto x:S)cout<<x<<' ';cout<<'\n';}}
