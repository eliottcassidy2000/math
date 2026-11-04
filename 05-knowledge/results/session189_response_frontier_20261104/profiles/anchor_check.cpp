#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
vector<int64_t> ends(const vector<unsigned>&r,int repeated){int n=r.size(),states=1;vector<int>k(n,1),stride(n);if(repeated>=0)k[repeated]=2;for(int i=0;i<n;i++){stride[i]=states;states*=k[i]+1;}vector<int64_t>d(states*n);for(int i=0;i<n;i++)d[stride[i]*n+i]=1;for(int st=1;st<states;st++){unsigned avail=0;for(int j=0;j<n;j++)if(st/stride[j]%(k[j]+1)<k[j])avail|=1<<j;for(int i=0;i<n;i++){auto v=d[st*n+i];if(!v)continue;unsigned z=r[i]&avail;while(z){int j=__builtin_ctz(z);z&=z-1;d[(st+stride[j])*n+j]+=v;}}}vector<int64_t> out(n);for(int i=0;i<n;i++)out[i]=d[(states-1)*n+i];return out;}
int64_t total(const vector<unsigned>&r,int j){int64_t z=0;for(auto v:ends(r,j))z+=v;return z;}
vector<unsigned> ext(vector<unsigned>r,int c){int n=r.size();r.push_back(((1<<n)-1)^c);for(int i=0;i<n;i++)if(c>>i&1)r[i]|=1<<n;return r;}
int main(){int nq;cin>>nq;for(int z=0;z<nq;z++){vector<unsigned>Q(7),U(7);for(auto&v:Q)cin>>v;for(int i=0;i<7;i++)U[i]=127^(1<<i)^Q[i];vector<vector<int64_t>> sk(7,vector<int64_t>(7));for(int j=0;j<7;j++){auto a=ends(Q,j),b=ends(U,j);for(int i=0;i<7;i++)sk[i][j]=b[i]-a[i];}for(auto a:sk)for(auto v:a)cout<<v<<" ";for(int i=0;i<7;i++)cout<<total(ext(Q,1<<i),-1)-total(ext(U,1<<i),-1)<<" ";for(int i=0;i<7;i++)for(int j=0;j<7;j++)cout<<total(ext(Q,1<<i),j)-total(ext(U,1<<i),j)<<" ";cout<<"\n";}}
