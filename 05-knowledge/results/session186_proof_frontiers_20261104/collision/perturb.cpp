#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using U=uint64_t;
vector<U> ending(const vector<int>& a){int n=a.size(),N=1<<n;vector<U>d(N*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(d[m*n+v]){int opts=a[v]&~m;while(opts){int b=opts&-opts;opts-=b;d[(m|b)*n+__builtin_ctz(unsigned(b))]+=d[m*n+v];}}return vector<U>(d.end()-n,d.end());}
vector<int> op(vector<int>a){for(int i=0;i<int(a.size());i++)a[i]^=((1<<a.size())-1)^(1<<i);return a;}
bool pairmod(const vector<int>&a){int n=a.size();for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(!((a[i]^a[j])&~((1<<i)|(1<<j))))return true;return false;}
int main(){
 array<array<int,6>,3> seeds={{{16,5,33,39,46,3},{16,5,49,39,42,3},{18,4,49,39,42,3}}};
 for(int r=6;r<=8;r++)for(int sj=0;sj<3;sj++){
 int n=5+r;vector<int>a(n);for(int i=0;i<5;i++){a[i]=seeds[sj][i]&31;for(int j=5;j<n;j++)if(seeds[sj][i]>>5&1)a[i]|=1<<j;else a[j]|=1<<i;}for(int i=5;i<n;i++)for(int j=i+1;j<n;j++)a[i]|=1<<j;
 auto E=ending(a),S=ending(op(a));U H=accumulate(E.begin(),E.end(),U(0));vector<pair<int,int>>edges;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.push_back({i,j});int m=edges.size();long long total[4]={},nopair[4]={},sameH[4]={},sameE[4]={},sameSE[4]={};
 auto flip=[&](int e){auto [i,j]=edges[e];a[i]^=1<<j;a[j]^=1<<i;};
 auto test=[&](int k){total[k]++;if(pairmod(a))return;nopair[k]++;auto e=ending(a);if(accumulate(e.begin(),e.end(),U(0))!=H)return;sameH[k]++;cout<<"SAMEH n="<<n<<" seed="<<sj<<" radius="<<k<<" rows=";for(auto x:a)cout<<x<<",";cout<<endl;if(e!=E)return;sameE[k]++;if(ending(op(a))!=S)return;sameSE[k]++;cout<<"SURVIVOR n="<<n<<" seed="<<sj<<" radius="<<k<<" rows=";for(auto x:a)cout<<x<<",";cout<<endl;};
 for(int i=0;i<m;i++){flip(i);test(1);for(int j=i+1;j<m;j++){flip(j);test(2);for(int k=j+1;k<m;k++){flip(k);test(3);flip(k);}flip(j);}flip(i);}
 for(int k=1;k<=3;k++)cout<<"n="<<n<<" seed="<<sj<<" H="<<H<<" radius="<<k<<" total="<<total[k]<<" no_pair="<<nopair[k]<<" same_H="<<sameH[k]<<" same_E="<<sameE[k]<<" same_SE="<<sameSE[k]<<endl;
 }
}
