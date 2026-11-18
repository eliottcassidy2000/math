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
  for(auto [z,h]:Ha)if(news.count(h)){
   auto a=rows(z,n);
   cout<<n<<"\t"<<h<<"\t";
   for(int i=0;i<n;i++)cout<<(i?",":"")<<a[i];
   cout<<"\n";
  }
  old=cur;
 }
}
