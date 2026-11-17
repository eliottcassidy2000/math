#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
struct Calc {
 int n,f; array<int,10>A{}; array<int,512>C{};
 void dfs(int root,int v,int S,int len,int avail){
    if(len>=3 && len%2 && ((A[v]>>root)&1))++C[S];
    int rem=A[v]&avail;
    while(rem){int t=rem&-rem;rem^=t;dfs(root,__builtin_ctz(t),S|t,len+1,avail^t);}
 }
 bool strong(){
    int full=(1<<n)-1;
    for(int i=0;i<n;i++){
      int reach=1<<i;
      for(int j=0;j<n;j++){
        int last=reach;
        for(int k=0;k<n;k++)if((reach>>k)&1)reach|=A[k];
        if(reach==last)break;
      }
      if(reach!=full)return false;
    }
    return true;
 }
 vector<int64_t> eval(){
    f=(1<<n)-1; C.fill(0); for(int u=0;u<n;u++)dfs(u,u,1<<u,1,f^((2<<u)-1));
    int K=(n-1)/2;
    vector<int64_t> g((K+1)*4);
    vector<int> V;
    for(int s=1;s<=f;s++)if(C[s]){
      int k=(__builtin_popcount(s)-1)/2;
      g[k*4+1]+=C[s];
      V.push_back(s);
    }
    for(unsigned a=0;a<V.size();a++)for(unsigned b=a+1;b<V.size();b++) {
      int x=V[a],y=V[b];
      if(x&y)continue;
      int z=x|y;
      g[((__builtin_popcount(z)-2)/2)*4+2]+=(int64_t)C[x]*C[y];
      if(n>=9 && __builtin_popcount(z)==6){
        int w=f^z;
        if(w>y && C[w])g[3*4+3]+=(int64_t)C[x]*C[y]*C[w];
      }
    }
    return g;
 }
};
int main(int argc,char**argv){
 if(argc!=3)return 2;
 int n=stoi(argv[1]);string mode=argv[2],line; uint64_t id=0;
 if(n<1||n>9)return 2;
 while(getline(cin,line)){
    if(line.empty()||line[0]<'0'||line[0]>'9')continue;
    Calc a; a.n=n; istringstream is(line);
    if(mode=="code"){
      is>>id; int k=n*(n-1)/2-1;
      for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k--)
        if((id>>k)&1)a.A[i]|=1<<j;else a.A[j]|=1<<i;
    } else if(mode=="rows"){
      ++id;
      for(int i=0;i<n;i++)if(!(is>>a.A[i]))return 2;
    } else return 2;
    for(int i=0;i<n;i++) {
        if(a.A[i]>>i&1)return 3;
        for(int j=0;j<i;j++)if(((a.A[i]>>j&1)+(a.A[j]>>i&1))!=1)return 3;
    }
    auto g=a.eval(); cout<<id<<' '<<a.strong();
    for(int k=1;2*k<n;k++)for(int c=1;c<=min(k,n-2*k);c++)cout<<' '<<g[4*k+c];
    cout<<'\n';
 }
}
