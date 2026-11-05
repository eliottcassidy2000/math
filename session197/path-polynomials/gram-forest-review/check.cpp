#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <stdexcept>
using namespace std;using G=vector<int>;using ll=long long;
void need(bool x,const char*s){if(!x)throw runtime_error(s);}G tr(int n){G a(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)a[i]|=1<<j;return a;}
G tour(int n,int c){G a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){bool f=c>>k&1;a[f?i:j]|=1<<(f?j:i);}return a;}
G cyclic(const G&a,const G&b,const G&c){int n=a.size(),m=b.size(),l=c.size();G out(n+m+l);for(int i=0;i<n;i++)out[i]=a[i]|(((1<<m)-1)<<n);for(int i=0;i<m;i++)out[n+i]=(b[i]<<n)|(((1<<l)-1)<<(n+m));for(int i=0;i<l;i++)out[n+m+i]=(c[i]<<(n+m))|((1<<n)-1);return out;}
vector<ll>ordered(const G&a){int n=a.size();vector<ll>v(n+1);vector<int>p(n);iota(p.begin(),p.end(),0);do{int bad=0;for(int i=0;i<n-1;i++)if(!(a[p[i]]>>p[i+1]&1))bad|=1<<i;for(int cut=0;cut<(1<<(n-1));cut++)if((cut&bad)==bad)v[1+__builtin_popcount((unsigned)cut)]++;}while(next_permutation(p.begin(),p.end()));return v;}
ll walks(const G&a){int n=a.size(),N=1<<n;ll H=0;for(int S=1;S<N;S++){vector<ll>d(n);for(int i=0;i<n;i++)d[i]=(S>>i)&1;for(int k=1;k<n;k++){vector<ll>e(n);for(int i=0;i<n;i++)if(d[i])for(int q=a[i]&S;q;q&=q-1)e[__builtin_ctz((unsigned)q)]+=d[i];d=e;}ll W=accumulate(d.begin(),d.end(),0LL);H+=((n-__builtin_popcount((unsigned)S))%2?-W:W);}return H;}
vector<vector<ll>>jacobi(vector<ll>f){int m=f.size()-1,L=m+1;f.resize(m+3);vector<vector<ll>>J(L,vector<ll>(L));for(int r=1;r<=L;r++){J[r-1][r-1]=f[r-1]+3*f[r]+f[r+1];if(r<L)J[r-1][r]=J[r][r-1]=f[r]+f[r+1];}return J;}
ll form(const vector<vector<ll>>&J,vector<ll>x,vector<ll>y){int L=J.size();x.resize(L+1);y.resize(L+1);ll s=0;for(int i=0;i<L;i++)for(int j=0;j<L;j++)s+=x[i+1]*J[i][j]*y[j+1];return s;}
int main(){vector<G>small;for(int n=1;n<=3;n++)for(int c=0;c<(1<<(n*(n-1)/2));c++)small.push_back(tour(n,c));int physical=0,gram=0,extra=0;for(auto&V:small){auto f=ordered(V);auto J=jacobi(f);for(auto&T:small)for(auto&U:small){auto x=ordered(T),y=ordered(U);ll H=walks(cyclic(T,U,V));need(H==form(J,x,y),"literal walk IE versus cover kernel");need(H*H<=form(J,x,x)*form(J,y,y),"physical CS");x.resize(J.size()+1);y.resize(J.size()+1);bool proportional=1;for(int r=1;r<=(int)J.size();r++)proportional&=x[r]*y[1]==y[r]*x[1];need((H*H==form(J,x,x)*form(J,y,y))==proportional,"CS equality condition");physical++;}for(int a=1;a<=(int)J.size();a++)for(int b=1;b<=(int)J.size();b++){need(walks(cyclic(tr(a),tr(b),V))==form(J,ordered(tr(a)),ordered(tr(b))),"actual full rank Gram entries");gram++;}}
for(int n:{4,5})for(int seed=0;seed<8;seed++){G V=tour(n,(seed*197+13)&((1<<(n*(n-1)/2))-1));auto f=ordered(V);auto J=jacobi(f);for(auto [a,b]:vector<array<int,2>>{{1,1},{1,2},{2,1},{2,2}}){need(walks(cyclic(tr(a),tr(b),V))==form(J,ordered(tr(a)),ordered(tr(b))),"larger arbitrary V kernel");extra++;}}
cout<<"PASS literal_physical_triples "<<physical<<" full_transitive_Gram_entries "<<gram<<" extra_order4_5_third_blocks "<<extra<<"\n";
}
