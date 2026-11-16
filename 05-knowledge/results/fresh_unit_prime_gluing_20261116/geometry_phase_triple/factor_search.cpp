#include<bits/stdc++.h>
using namespace std;
using V=vector<int>;
struct Data{
 int n,F;V P,h; Data(const V&a):n(a.size()),F((1<<n)-1),P((F+1)*n),h(F+1){
 h[0]=1;
 for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){
   int i=__builtin_ctz(b),R=S^(1<<i),v=!R;
   for(int x=R&a[i];x;x&=x-1)v+=P[R*n+__builtin_ctz(x)];
   P[S*n+i]=v;h[S]+=v;
 }}
};
V join(V a,const V&b){
 int n=a.size(),added=((1<<b.size())-1)<<n;
 for(auto&x:a)x|=added;
 for(int x:b)a.push_back(x<<n);
 return a;
}
int odd(int a){return a>>__builtin_ctz(a);}
int sums(const Data &d){int p=0;for(int i=0;i<=d.F;i++)p+=d.h[i]*d.h[d.F^i];return p;}
struct Plan{int t,g,h;V pre,post;int xmult;};
using Key=pair<int,int>;
map<int,vector<Plan>> plans;
map<Key,array<pair<V,int>,3>> found;
int success=0,N=0,accept=0;
void examine(V a,Data &d){
 auto j=plans.find(d.h[d.F]);if(j==plans.end())return; ++accept;
 int n=d.n,F=d.F, h=d.h[F];V X(F+1),w(n),ev(n*n);X[0]=h;
 for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){int i=__builtin_ctz(b);w[i]+=d.P[S*n+i]*d.h[F^S];}
 V p((F+1)*n);
 for(int st=0;st<n;st++){
  fill(p.begin(),p.end(),0);p[(1<<st)*n+st]=1;
  for(int S=1<<st;S<=F;S++)if(S>>st&1)for(int b=S;b;b&=b-1){
   int k=__builtin_ctz(b),v=p[S*n+k];if(!v)continue;
   if((__builtin_popcount(S)&1)==0)ev[st*n+k]+=v*d.h[F^S];
   for(int u=(F^S)&a[k];u;u&=u-1)p[(S|(u&-u))*n+__builtin_ctz(u)]+=v;
  }
 }
 V cut(F+1),tot(n);
 for(int i=0;i<n;i++)for(int k=0;k<n;k++)tot[i]+=ev[i*n+k];
 for(int q=1;q<=F;q++){
   int i=__builtin_ctz(q),s=q&(q-1);X[q]=X[s]+w[i];cut[q]=cut[s]+tot[i];
   for(int b=s;b;b&=b-1){int k=__builtin_ctz(b);cut[q]-=ev[i*n+k]+ev[k*n+i];}
   int H=h+2*cut[q];
   for(auto pl:j->second)if(H==pl.h){
    auto &slot=found[{pl.g,odd(X[q]*pl.xmult)}];if(!slot[pl.t].first.empty())continue;
    V c=a;for(int k=0;k<n;k++)if(!(q>>k&1))c[k]|=1<<n;c.push_back(q);
    auto b=join(join(pl.pre,c),pl.post);slot[pl.t]={b,n+pl.pre.size()};
    if(!slot[0].first.empty()&&!slot[1].first.empty()&&!slot[2].first.empty()){
     cout<<"TRIPLE g "<<pl.g<<" odd "<<odd(X[q]*pl.xmult)<<" nprocessed "<<N<<"\n";
     for(auto &[rr,r]:slot){cout<<"root "<<r<<" rows";for(int x:rr)cout<<" "<<x;cout<<"\n";} success=1;return;
    }
   }
 }
}
int main(int argc,char **argv){
 int rounds=argc>1?atoi(argv[1]):100000;
 map<int,map<int,V>> palette;palette[1][1]={};
 for(int n=3;n<=5;n++){
  vector<pair<int,int>> e;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)e.emplace_back(i,j);
  for(int z=0;z<(1<<(int)e.size());z++){
   V a(n);int i=0;for(auto [u,v]:e){if(z>>i&1)a[u]|=1<<v;else a[v]|=1<<u;++i;}
   Data d(a);int h=d.h[d.F];if(h==3||h==5||h==11||h==13||h==15)palette[h].try_emplace(odd(sums(d)),a);
  }
 }
 for(auto [k,a]:palette[3])for(auto [l,b]:palette[5])palette[15].try_emplace(k*l,join(a,b));
 vector<V> factors={{3,5,15},{11},{13}};V P={6015,6017,6019};
 for(int g=1;g<=6015;g+=2)if(g%429==0||g%715==0){
  for(int t=0;t<3;t++)for(int f:factors[t])if(g%f==0)for(int l:{1,3,5,11,13,15})if(f%l==0){
   for(auto &[x,a]:palette[l])for(auto &[y,b]:palette[f/l])plans[g/f].push_back({t,g,P[t]/f,a,b,l*y});
  }
 }
 mt19937 rng(17009);
 map<int,vector<V>> pool;
 for(int it=0;it<rounds&&!success;it++){
  int n=6+rng()%4;if(it%10==0)n=10;
  V a(n);
  auto &bank=pool[n];
  if(!bank.empty()&&it%3!=0){
   a=bank[rng()%bank.size()];int edits=1+rng()%6;
   while(edits--){int i=rng()%n,j=rng()%(n-1);j+=j>=i;a[i]^=1<<j;a[j]^=1<<i;}
  }else for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(rng()&1)a[i]|=1<<j;else a[j]|=1<<i;
  Data d(a);++N;
  if(plans.count(d.h[d.F])&&bank.size()<3000)bank.push_back(a);
  examine(a,d);
 }
 int pairs=0;map<int,array<int,3>> hits;
 for(auto&[k,a]:found){int s=0;for(int i=0;i<3;i++)if(!a[i].first.empty()){++s;hits[k.first][i]++;}pairs+=s>=2;}
 cout<<"processed "<<N<<" joined "<<accept<<" keys "<<found.size()<<" pairs "<<pairs<<"\n";
 for(auto &[g,a]:hits)cout<<g<<" "<<a[0]<<" "<<a[1]<<" "<<a[2]<<"\n";
}
