#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <cstdint>
using namespace std;using U=uint64_t;
struct Result{U low,high,merged;int aut;bool sc;};
// Row-major encoding deliberately differs from the producer's column-major encoding.
Result canon(const vector<int>&r){int n=r.size();vector<int>p(n);iota(p.begin(),p.end(),0);U lo=~U(0),hi=0;int count=0;do{U v=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)v=(v<<1)|((r[p[i]]>>p[j])&1);if(v<lo){lo=v;count=1;}else if(v==lo)count++;hi=max(hi,v);}while(next_permutation(p.begin(),p.end()));U opp=((U(1)<<(n*(n-1)/2))-1)^hi;return{lo,hi,min(lo,opp),count,lo==opp};}
int main(){int cases;cin>>cases;cout<<'[';for(int id=0;id<cases;id++){int n;cin>>n;vector<int>r(n);for(int&i:r)cin>>i;auto z=canon(r);if(id)cout<<',';cout<<"{\"merged\":"<<z.merged<<",\"aut\":"<<z.aut<<",\"self_converse\":"<<z.sc<<",\"cards\":[";for(int a=0;a<n;a++){vector<int>g;for(int i=0;i<n;i++)if(i!=a){int row=0,jj=0;for(int j=0;j<n;j++)if(j!=a){if(r[i]>>j&1)row|=1<<jj;jj++;}g.push_back(row);}auto c=canon(g);if(a)cout<<',';cout<<"{\"deleted\":"<<a<<",\"merged\":"<<c.merged<<'}';}cout<<"]}";}cout<<"]\n";}
