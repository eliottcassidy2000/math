#include<bits/stdc++.h>
using namespace std;using U=uint64_t;using V=vector<int>;
int nout=8;bool witnesses=false;
U canon(V a){
 int n=a.size();vector<V> b(n);
 for(int i=0;i<n;i++)b[__builtin_popcount(a[i])].push_back(i);
 V p;vector<pair<int,int>> seg;
 for(auto& x:b)if(x.size()){
  int L=p.size();p.insert(p.end(),x.begin(),x.end());
  if(x.size()>1)seg.push_back({L,(int)p.size()});
 }
 U best=~0ull;
 auto visit=[&](auto&&self,size_t i)->void{
  if(i==seg.size()){
   U z=0;int k=0;
   for(int j=0;j<n;j++)for(int l=j+1;l<n;l++,k++)z|=U(a[p[j]]>>p[l]&1)<<k;
   best=min(best,z);return;
  }
  auto [a,b]=seg[i];do{self(self,i+1);}while(next_permutation(p.begin()+a,p.begin()+b));
 };
 visit(visit,0);return best;
}
V rows(U z,int n){
 V a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){
  if(z>>k&1)a[i]|=1<<j;else a[j]|=1<<i;
 }return a;
}
V add(V a,int w){int n=a.size();for(int i=0;i<n;i++)if(!(w>>i&1))a[i]|=1<<n;a.push_back(w);return a;}
V remove(V a,int u){
 V d;for(int i=0;i<(int)a.size();i++)if(i!=u)d.push_back((a[i]&((1<<u)-1))|((a[i]>>(u+1))<<u));
 return d;
}
int H(const V&a){
 int n=a.size(),M=1<<n;int d[512][9]={};
 for(int i=0;i<n;i++)d[1<<i][i]=1;
 for(int s=1;s<M;s++)for(int t=s;t;t&=t-1){
  int i=__builtin_ctz(t),v=d[s][i];if(!v)continue;
  for(int k=a[i]&~s;k;k&=k-1)d[s|(k&-k)][__builtin_ctz(k)]+=v;
 }
 return accumulate(d[M-1],d[M-1]+n,0);
}
bool strong(const V&a){
 int n=a.size(),w=1,y=1;for(int j=0;j<n;j++)for(int i=0;i<n;i++){
  if(w>>i&1)w|=a[i];if(a[i]&y)y|=1<<i;
 }
 return w==((1<<n)-1)&&w==y;
}
int main(int argc,char**argv){
 if(argc>1)nout=atoi(argv[1]);if(nout<1 || nout>9)return 2; if(argc>2)witnesses=(string(argv[2])=="--witnesses");
 set<U> old={0};map<int,int> birth={{1,1}};
 map<U,int> Ho={{0,1}};set<U> Boo={0};
 if(witnesses)cout<<"order\tcount\toutgoing_masks\n1\t1\t0\n";
 for(int n=2;n<=nout;n++){
  set<U> cur;for(auto z:old)for(int u=0;u<(1<<(n-1));u++)cur.insert(canon(add(rows(z,n-1),u)));
  map<U,int> Ha;
  set<int> news;
  for(auto z:cur){int h=H(rows(z,n));Ha[z]=h;if(!birth.count(h))news.insert(h);}
  for(auto h:news)birth[h]=n;
  if(witnesses){
   map<int,U> reps;for(auto [z,h]:Ha)if(news.count(h))reps.try_emplace(h,z);
   for(auto [h,z]:reps){
    cout<<n<<"\t"<<h<<"\t";auto a=rows(z,n);
    for(int i=0;i<n;i++)cout<<(i?",":"")<<a[i];cout<<"\n";
   }
   old=cur;continue;
  }

  set<U> Bo, esc;
  for(auto z:cur)if(news.count(Ha[z])){
   Bo.insert(z);
   bool hit=false;auto a=rows(z,n);
   for(int i=0;i<n;i++)if(Boo.count(canon(remove(a,i))))hit=true;
   if(!hit)esc.insert(z);
  }
  int strongly=0, bbflip=0, bnflip=0, bdeath=0, dsame=0, nf=0; set<pair<U,U>> flips;
  map<U,set<U>> parents;
  for(auto z:cur){
   auto a=rows(z,n);
   if(Bo.count(z)&&strong(a))strongly++;
   for(int i=0;i<n;i++)parents[z].insert(canon(remove(a,i)));
   for(auto y:parents[z]){
    if(Bo.count(z)&&Boo.count(y))bdeath++;
    if(Ho.at(y)==Ha[z])dsame++;
   }
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
    a[i]^=1<<j;a[j]^=1<<i;U x=canon(a);
    if(x!=z)flips.insert(minmax(z,x));else nf++;
    a[i]^=1<<j;a[j]^=1<<i;
   }
  }
  for(auto [x,y]:flips){int l=Bo.count(x)+Bo.count(y);if(l==2)bbflip++;else if(l==1)bnflip++;}
  cout<<"LEVEL "<<n<<" classes "<<cur.size()<<" new_values "<<news.size()<<" birth_classes "<<Bo.size()<<" strong_birth "<<strongly<<" parentless_birth "<<esc.size()<<" cross_edges "<<bdeath<<" same_count_deletions "<<dsame<<" simple_flips "<<flips.size()<<" birth_birth_flips "<<bbflip<<" birth_old_flips "<<bnflip<<" neutral_isom_flips_labeled "<<nf<<"\n";
  map<U,set<U>> adj;
  set<int> uni; int high=0, unit=0, directionmax=0, oldmax=0,minnew=INT_MAX;
  for(auto [y,h]:Ho)oldmax=max(h,oldmax);
  for(auto z:Bo){if(Ha[z]>oldmax)high++;minnew=min(minnew,Ha[z]);}
  for(auto [x,y]:flips)if(Bo.count(x)&&Bo.count(y)){adj[x].insert(y);adj[y].insert(x);if(abs(Ha[x]-Ha[y])==2){unit++;uni.insert(min(Ha[x],Ha[y]));}}
  set<U> alive=Bo;vector<int> comps;map<U,int> anchor,dist;
  queue<U> wave;
  for(auto a:Bo)if(!esc.count(a)){dist[a]=0;wave.push(a);}
  while(!wave.empty()){auto x=wave.front();wave.pop();for(auto y:adj[x])if(!dist.count(y)){dist[y]=dist[x]+1;wave.push(y);}}
  while(!alive.empty()){
   vector<U> visit={*alive.begin()}; alive.erase(visit[0]);
   for(size_t at=0;at<visit.size();at++)for(auto y:adj[visit[at]])if(alive.erase(y))visit.push_back(y);
   comps.push_back(visit.size());
  }
  for(auto x:Bo){bool high=false;for(auto y:adj[x])if(Ha[y]>Ha[x])high=true;if(!high)directionmax++;}
  sort(comps.begin(),comps.end(),greater<int>());
  cout<<"SHAPE "<<n<<" old_max "<<oldmax<<" high_birth "<<high<<" unit_birth_flips "<<unit<<" weak_local_maxima "<<directionmax<<" components ";
  for(auto x:comps)cout<<x<<",";
  cout<<" relay_dist ";
  for(auto x:esc)if(dist.count(x))cout<<dist[x]<<",";else cout<<"inf,";
  cout<<"\n";
  int valuepairs=0; vector<int> missed;
  for(auto q:news)if(news.count(q+2)){valuepairs++;if(!uni.count(q))missed.push_back(q);}
  cout<<"UNITS "<<n<<" adjacent_values "<<valuepairs<<" physically_joined "<<uni.size()<<" missing_first ";
  for(size_t k=0;k<min<size_t>(10,missed.size());k++)cout<<missed[k]<<",";
  cout<<"\n";
  for(auto x:esc){
   cout<<"ORPH "<<n<<" "<<x<<" H "<<Ha[x]<<" old counts ";
   set<int> p;for(auto y:parents[x])p.insert(Ho.at(y));for(auto y:p)cout<<y<<",";
   int near=0;set<int> nearh;for(auto [z,w]:flips)if(z==x||w==x){auto y=z==x?w:z;if(Bo.count(y)&&!esc.count(y)){near++;nearh.insert(Ha[y]);}}
   cout<<" anchored birth neighbors "<<near<<" counts ";for(auto h:nearh)cout<<h<<",";
   cout<<" rows ";for(auto r:rows(x,n))cout<<r<<",";cout<<"\n";
  }
  if(n<=8)for(auto z:Bo){
   cout<<"ATLAS "<<n<<" "<<z<<" "<<Ha[z]<<" ";for(int x:rows(z,n))cout<<x<<",";cout<<" ";
   for(auto p:parents[z])cout<<p<<":"<<Ho.at(p)<<":"<<Boo.count(p)<<",";cout<<"\n";
  }
  old=cur;Ho=Ha;Boo=Bo;
 }
}
