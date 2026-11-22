#include <bits/stdc++.h>
using namespace std; using ll=long long; using I=__int128;
using R=vector<int>; struct K{ll c[10][10]={};};
K profile(const R& g,int u,int v){
 int n=g.size(),P=0, bit[10]={}, full=(1<<(n-2))-1, d=n;
 for(int i=0;i<n;i++) if(i!=u&&i!=v) bit[i]=1<<P++;
 vector<ll> a((full+1)*d*d*n);
 auto D=[&](int S,int i,int j,int a0)->ll&{return a[(((S*d)+i)*d+j)*n+a0];};
 for(int i=0;i<n;i++) D(bit[i],i==u,i==v,i)=1;
 K k;
 for(int S=0;S<=full;S++) for(int i=0;i<d;i++) for(int j=0;j<d;j++) for(int s=0;s<n;s++){
   ll N=D(S,i,j,s); if(!N)continue;
   if(S==full && i&&j)k.c[i][j]+=N;
   for(int z=g[s];z;z&=z-1){
    int t=__builtin_ctz(z);
    if(bit[t]){if(!(S&bit[t]))D(S|bit[t],i,j,t)+=N;}
    else if(t==u){if(i<d-1)D(S,i+1,j,t)+=N;}
    else if(j<d-1)D(S,i,j+1,t)+=N;
   }
 }
 return k;
}
bool strong(const R&g){int n=g.size(),all=(1<<n)-1;for(int i=0;i<n;i++){int s=1<<i;for(int j=0;j<n;j++)for(int k=0;k<n;k++)if(s>>k&1)s|=g[k];if(s!=all)return false;}return true;}
struct f{int a,b,d;};
vector<f> ray(const K&k){
 ll a=k.c[2][2], b=k.c[2][1]-2*a,c=k.c[1][2]-2*a,d=k.c[1][1]-2*k.c[2][1]-2*k.c[1][2]+4*a;
 if(!a || (a&(a-1)) || abs(d)>25)return {};
 for(int i=1;i<10;i++)for(int j=1;j<10;j++)if(max(i,j)>=3&&k.c[i][j])return {};
 vector<f> ans;int sum=1+__builtin_ctzll(a);
 for(int x=-30;x<=30;x++){int y=sum-x,h=max({0,-x,-y}); if(((I)b*(I(1)<<(y+h)))+((I)c*(I(1)<<(x+h)))==-(I(a)<<h))ans.push_back({x,y,(int)d});}
 return ans;
}
void examine(const R&g, bool near,int onlyu=-1,int onlyv=-1){
 int n=g.size(); if(!strong(g))return;
 for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
  if((onlyu>=0&&u!=onlyu)||(onlyv>=0&&v!=onlyv))continue;
  K k=profile(g,u,v);
  auto xs=ray(k);
  for(auto [a,b,d]:xs) {
   if(near || abs(d)==1) {
    cout<<"ray "<<n<<" "<<u<<" "<<v<<" "<<a<<" "<<b<<" "<<d<<" rows";
    for(int g0:g)cout<<" "<<g0;
    cout<<" signature "<<k.c[1][1]<<","<<k.c[1][2]<<","<<k.c[2][1]<<","<<k.c[2][2];
    cout<<"\n";
   }
   if(d!=-1&&d!=1)continue;
   for(int x=0;x<n;x++)for(int y=x+1;y<n;y++){
    if(x==u||x==v||y==u||y==v)continue;
    R h=g; h[x]^=1<<y;h[y]^=1<<x; if(!strong(h))continue;
    for(auto [i,j,e]:ray(profile(h,u,v)))if(i==a&&j==b&&e==-d)
     cout<<"STRADDLE "<<u<<","<<v<<" "<<a<<","<<b<<" flip "<<x<<","<<y<<"\n";
   }
  }
 }
}
int main(int argc,char**argv){
 int limit=argc>1?atoi(argv[1]):8;
 string deck=argc>2?argv[2]:"session212/cycle-levels";
 if(limit==0){
   for(R g: {R{990,1020,1016,976,928,137,688,768,608,32}, R{998,1020,1016,977,929,456,656,768,576,32}}){
    examine(g,true,1,9);
    int n=g.size();
    vector<pair<int,int>> E;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)E.emplace_back(i,j);
    for (int i=0;i<int(E.size());i++) {auto [a,b]=E[i];g[a]^=1<<b;g[b]^=1<<a;examine(g,true,1,9);
      for(int j=i+1;j<int(E.size());j++){auto [c,d]=E[j];g[c]^=1<<d;g[d]^=1<<c;examine(g,true,1,9);g[c]^=1<<d;g[d]^=1<<c;}
      g[a]^=1<<b;g[b]^=1<<a;
    }
   }
   return 0;
 }
 for(int n=3;n<=limit;n++){
  ifstream in(deck+"/classes"+to_string(n)+".rows");R g(n);
  if(!in)throw runtime_error("cannot read source deck");
  long long count=0;
  while(in>>g[0]){for(int i=1;i<n;i++)if(!(in>>g[i]))throw runtime_error("truncated source deck");examine(g,false);++count;}
  const int expected[]={0,1,1,2,4,12,56,456,6880};
  if(n<=8&&count!=expected[n])throw runtime_error("incorrect deck cardinality");
  cerr<<"order "<<n<<" source classes "<<count<<"\n";
 }
}
