#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using I=long long;
using Adj=std::vector<std::vector<int>>;
using Row=std::array<I,5>; // h,q,E-L,J
Adj make(int n,int code){
 Adj a(n,std::vector<int>(n));
 for(int i=0,k=0;i<n;i++)for(int j=i+1;j<n;j++){
   a[i][j]=code>>k++&1;a[j][i]=1-a[i][j];
 }
 return a;
}
bool strong(const Adj&a){
 int n=a.size();
 for(int fun=0;fun<2;fun++){
  int got=1,old=0;
  while(old!=got){old=got;for(int i=0;i<n;i++)if(got>>i&1)for(int j=0;j<n;j++)if(fun?a[j][i]:a[i][j])got|=1<<j;}
  if(got!=(1<<n)-1)return false;
 }
 return true;
}
std::vector<Row> profile(const Adj&a,int roots){
 const int n=a.size(),V=1<<n,F=V-1;
 I en[4096][12]={},st[4096][12]={},all[4096]={},E[12]={},L[12]={},J[12][12]={},q=0;
 all[0]=1;
 for(int i=0;i<n;i++)en[1<<i][i]=st[1<<i][i]=1;
 for(int s=1;s<V;s++)for(int i=0;i<n;i++)if(s>>i&1){
  all[s]+=en[s][i];for(int j=0;j<n;j++)if(!(s>>j&1)){
   if(a[i][j])en[s|1<<j][j]+=en[s][i];else st[s|1<<j][j]+=st[s][i];
  }
 }
 for(int s=0;s<V;s++){
  q+=all[s]*all[F^s];for(int u=0;u<roots;u++){
   E[u]+=en[s][u]*all[F^s];L[u]+=st[s][u]*all[F^s];
   for(int v=0;v<roots;v++)J[u][v]+=en[s][u]*st[F^s][v];
  }
 }
 std::vector<Row> ans(1<<roots);
 for(int mask=0;mask<(1<<roots);mask++){
  I ee=0,ll=0,jj=0;
  for(int u=0;u<roots;u++)if(mask>>u&1){
   ee+=E[u];ll+=L[u];for(int v=0;v<roots;v++)if(mask>>v&1)jj+=J[u][v];
  }
  ans[mask]={all[F],q,ee,ll,jj};
 }
 return ans;
}
Adj blow(const Adj&a,int w,int m){
 const int n=a.size(); std::vector<int> own;
 for(int u=0;u<n;u++)if(u!=w)own.push_back(u);
 for(int i=0;i<m;i++)own.push_back(w);
 Adj b(own.size(),std::vector<int>(own.size()));
 for(int i=0;i<(int)own.size();i++)for(int j=0;j<(int)own.size();j++)b[i][j]=own[i]==own[j]?i<j:a[own[i]][own[j]];
 return b;
}
int main(){
 int longest=0;
 for(int n=3;n<=7;n++){
  std::vector<bool> used(1<<(n*(n-1)/2));
  int count=0;
  for(int code=0;code<(int)used.size();code++)if(!used[code]){
   auto a=make(n,code);
   std::vector<int> p(n);std::iota(p.begin(),p.end(),0);
   do{int c=0,k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)c|=a[p[i]][p[j]]<<k++;used[c]=true;}while(std::next_permutation(p.begin(),p.end()));
   if(!strong(a))continue;count++;
   const int M=2+(n-1)/2;
   for(int w=0;w<n;w++){
    std::vector<std::vector<Row>> lib;
    for(int m=1;m<=M;m++)lib.push_back(profile(blow(a,w,m),n-1));

    std::map<std::vector<I>,std::map<I,std::array<int,3>>> core;
    for(int v=0;v<(1<<(n-1));v++)for(int a=1;a<=4;a++)for(int b=1;b<=4;b++){
     std::vector<I> key;I base=0; 
     for(int m=0;m<M;m++){
      auto [h,q,E,L,J]=lib[m][v];
      I H=h+q*(1LL<<(a+b))+E*(1LL<<a)-L*(1LL<<b)-J;
      if(m)key.push_back(H-base);else base=H;
     }
     core[key].emplace(base,std::array<int,3>{v,a,b});
    }
    for(auto&[key,table]:core)for(auto i=table.begin();i!=table.end();i++){
     auto j=i;int L=1;for(auto j=std::next(i);j!=table.end()&&j->first==i->first+2*L;j++)L++;
     if(L>longest){longest=L;printf("n=%d code=%d port=%d len=%d\n",n,code,w,L);
      for(int l=0;l<L;l++,j++){auto [v,a,b]=j->second;printf("mask=%d a=%d b=%d seed=%lld ",v,a,b,j->first);
       for(int m=0;m<M;m++){auto [h,q,E,L,D]=lib[m][v];printf("(%lld,%lld,%lld,%lld,%lld)",h,q,E,L,D);}puts("");}
      fflush(stdout);
     }
    }
   }
  }
  fprintf(stderr,"n=%d strong classes=%d\n",n,count);
 }
}
