#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;using ll=long long;using Graph=vector<int>;
void need(bool t,const char*s){if(!t)throw runtime_error(s);}
pair<Graph,int> build(vector<int>js){Graph a{0};int root=0;for(int j:js){int n=a.size();for(int&i:a)i|=((1<<(j+1))-1)<<n;a.resize(n+j+1);a[root]^=1<<n;a[n]|=1<<root;for(int x=1;x<=j;x++){a[n+x]|=1<<n;for(int y=x+1;y<=j;y++)a[n+x]|=1<<(n+y);}root=n;}return{a,root};}
array<ll,5> count(const Graph&a,int r){int n=a.size(),N=1<<n;vector<array<ll,16>>dp(N),st(N);vector<ll>h(N);h[0]=1;for(int S=1;S<N;S++){for(int v=0;v<n;v++)if(S>>v&1){int U=S^(1<<v);if(!U)dp[S][v]=st[S][v]=1;else for(int w=0;w<n;w++)if(U>>w&1){if(a[w]>>v&1)dp[S][v]+=dp[U][w];else st[S][v]+=st[U][w];}h[S]+=dp[S][v];}}ll R=0;for(int S=1;S<N;S++)if(S>>r&1)R+=st[S][r]*h[(N-1)^S];return{h[N-1],dp[N-1][r],R,st[N-1][r],h[(N-1)^(1<<r)]};}
ll cycles(const Graph&a){int n=a.size();function<ll(int,int)>walk=[&](int S,int v)->ll{if(S==(1<<n)-1)return a[v]&1;ll c=0;for(int bits=a[v]&~S;bits;bits&=bits-1){int b=bits&-bits;c+=walk(S|b,__builtin_ctz((unsigned)b));}return c;};return walk(1,0);}
int main(){int controls=0;cout<<"{\n  \"status\":\"PASS\",\n  \"witnesses\":[\n";bool first=true;for(int m=1;m<=4;m++){ll P=1;for(int j=1;j<=m;j++)P*=((1<<j)+1);for(int last=1;last<=m;last++){vector<int>js;for(int j=1;j<=m;j++)if(j!=last)js.push_back(j);js.push_back(last);auto [a,r]=build(js);need(cycles(a)==1,"unique core cycle");int n=a.size();for(int&i:a)i|=1<<n;a.push_back(0);auto z=count(a,r);ll d=P/((1<<last)+1);need(z==array<ll,5>{P,0,2*P,1,d},"weak ambiguity profile");need((int)a.size()==2+m*(m+3)/2,"same order");if(!first)cout<<",\n";first=false;cout<<"    {\"m\":"<<m<<",\"last\":"<<last<<",\"order\":"<<a.size()<<",\"H\":"<<P<<",\"R\":"<<2*P<<",\"d\":"<<d<<",\"root\":"<<r<<",\"rows\":[";for(int i=0;i<(int)a.size();i++)cout<<(i?",":"")<<a[i];cout<<"]}";controls++;}}cout<<"\n  ],\n  \"direct_graph_controls\":"<<controls<<"\n}\n";}
