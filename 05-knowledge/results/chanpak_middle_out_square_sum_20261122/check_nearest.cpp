#include<bits/stdc++.h>
using namespace std;using I=long long;using R=vector<int>;
R flip(R a,int u,int v){a[u]^=1<<v;a[v]^=1<<u;return a;}
R extend(R a,int s){for(int i=0;i<(int)a.size();i++)if(!(s>>i&1))a[i]|=1<<a.size();a.push_back(s);return a;}
void dp(const R&a,array<array<I,10>,1024>& f){
 int n=a.size(),m=1<<n;
 for(int s=0;s<m;s++)f[s].fill(0);
 for(int i=0;i<n;i++)f[1<<i][i]=1;
 for(int s=1;s<m-1;s++)for(int i=0;i<n;i++)if(f[s][i])for(int t=a[i]&~s;t;t&=t-1){int j=__builtin_ctz(t);f[s|1<<j][j]+=f[s][i];}
}
I hp(const R&a){
 static array<array<I,10>,1024> f;
 dp(a,f); return accumulate(f[(1<<a.size())-1].begin(),f[(1<<a.size())-1].end(),0LL);
}
vector<I> all_apex(const R&a){
 int n=a.size(),m=1<<n;R b(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i]>>j&1)b[j]|=1<<i;
 static array<array<I,10>,1024> E,B;
 dp(a,E); dp(b,B);
 I w[10][10]{},L[10]{};
 for(int s=1;s<m-1;s++)for(int u=s;u;u&=u-1){
  int i=__builtin_ctz(u);if(!E[s][i])continue;
  for(int v=(m-1)^s;v;v&=v-1){int j=__builtin_ctz(v);w[i][j]+=E[s][i]*B[(m-1)^s][j];}
 }
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){w[i][j]+=w[j][i];w[j][i]=w[i][j];L[(a[i]>>j&1)?i:j]+=w[i][j];}
 vector<I> ans(m);
 ans[0]=accumulate(E[m-1].begin(),E[m-1].end(),0LL);
 for(int s=1;s<m;s++){
  int i=__builtin_ctz(s),q=s&(s-1);I v=ans[q]+L[i];
  for(int t=q;t;t&=t-1)v-=w[i][__builtin_ctz(t)];
  ans[s]=v;
 }
 return ans;
}
struct Best {I mid=-1; vector<int> path;R rows;int i=-1,j=-1;};
int main(int argc,char** argv){
 int p=atoi(argv[1]),k=atoi(argv[2]);I q=1LL<<p,T=q*(2*q-1);Best best[2][2];I checked=0;
 set<R> seen;
 int pp,i,j;
 while(cin>>pp>>i>>j){
  int n=argc>3?atoi(argv[3]):(pp==3?6:8);R a(n);for(int&v:a)cin>>v;
  if(pp!=p)continue;
  if(p==3&&n==6&&seen.find(a)!=seen.end())continue;
  seen.insert(a);
  auto f=[&](const R& b, vector<int> prefix,auto &&f)->void{
   if((int)prefix.size()+1<k){
    for(int s=0;s<(1<<(int)b.size());s++){
     auto r=extend(b,s);
     auto d=prefix;d.push_back(s);f(r,d,f);
    }
    return;
   }
   auto low=all_apex(b),high=all_apex(flip(b,i,j));checked+=low.size();
   for(int s=0;s<(int)low.size();s++){
    I dif=high[s]-low[s];if(abs(dif)!=2)continue;
    I mid=(high[s]+low[s])/2;
    int side=mid>T, rev=dif<0;
    if(best[side][rev].mid==-1||abs(mid-T)<abs(best[side][rev].mid-T)){
     auto path=prefix;path.push_back(s);
     auto c=extend(b,s);assert(hp(c)==low[s]&&hp(flip(c,i,j))==high[s]);
     best[side][rev]={mid,path,a,i,j};
    }
   }
  };
  f(a,{},f);
 }
 for(int side=0;side<2;side++)for(int rev=0;rev<2;rev++){
  auto z=best[side][rev];
  cout<<"side "<<side<<" reverse "<<rev<<" center "<<z.mid<<" arc "<<z.i<<","<<z.j<<"  tags ";for(int w:z.path)cout<<w<<" ";cout<<" rows ";for(int w:z.rows)cout<<w<<" ";cout<<"\n";
 }cout<<"checks "<<checked<<"\n";
}
