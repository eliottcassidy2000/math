#include <bits/stdc++.h>
using namespace std;using V=vector<int>;
int H(V const& a){
 int n=a.size(),N=1<<n;V dp(N*n);
 for(int i=0;i<n;i++)dp[(1<<i)*n+i]=1;
 for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(int h=dp[s*n+i]){
  for(int w=a[i]&~s;w;w&=w-1){int b=w&-w;dp[(s|b)*n+__builtin_ctz(b)]+=h;}
 }
 return accumulate(dp.end()-n,dp.end(),0);
}
bool strong(V const&a){
 for(int i=0;i<(int)a.size();i++){
  int f=1<<i;
  for(;;){int g=f;for(int v=f;v;v&=v-1)g|=a[__builtin_ctz(v)];if(f==g)break;f=g;}
  if(f!=(1<<(int)a.size())-1)return false;
 }
 return true;
}
V canon(V const&a){
 V p(a.size());iota(p.begin(),p.end(),0);V ans=a;
 do{V x(a.size());
 for(int i=0;i<(int)p.size();i++)for(int j=0;j<(int)p.size();j++)if(a[p[i]]>>p[j]&1)x[i]|=1<<j;
 ans=min(ans,x);
 }while(next_permutation(p.begin(),p.end()));
 return ans;
}
map<int,set<V>> factors(int g,const char* reps,const char* reps8){
 map<int,set<V>> C;
 for(int n=3;n<=6;n++)for(int bits=0;bits<(1<<(n*(n-1)/2));bits++){
  V a(n);int k=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(bits>>k++&1)a[i]|=1<<j;else a[j]|=1<<i;
  int h=H(a);if(h>1&&g%h==0&&strong(a))C[h].insert(canon(a));
 }
 ifstream in(reps);
 V a(7);int count=0;
 while(in>>a[0]){for(int j=1;j<7;j++)if(!(in>>a[j]))abort();count++;int h=H(a);if(h>1&&g%h==0&&strong(a))C[h].insert(canon(a));}
 if(count!=456)abort();
 if(g>=45){
  ifstream in(reps8); V a(8); count=0;
  while(in>>a[0]){for(int j=1;j<8;j++)if(!(in>>a[j]))abort();count++;int h=H(a);if(h>1&&g%h==0&&strong(a))C[h].insert(a);}
  if(count!=6880)abort();
 }
 return C;
}
V ordinal(vector<V> const& Cs){
 int n=0;for(auto&c:Cs)n+=c.size();
 V a;int m=0;
 for(auto&c:Cs){for(int x:c)a.push_back(x<<m|((1<<n)-(1<<(m+c.size()))));m+=c.size();}
 return a;
}
struct W{int t,h,U;V D;};
using Bank=map<int,map<int,W>>;

struct P { int e=0,o=0; P operator*(P b)const{return {e*b.e+o*b.o,e*b.o+o*b.e};}};
struct Local{
 int h,n; vector<P> s,t,M;P opt;
 Local(V const&a):n(a.size()),s(n),t(n),M(n*n){
  int N=1<<n,nn=n*n;
  V H(N),D(N*nn);H[0]=1;
  for(int i=0;i<n;i++)D[(1<<i)*nn+i*n+i]=1;
  for(int z=1;z<N;z++)for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   int c=D[z*nn+i*n+j];H[z]+=c;if(!c)continue;
   for(int w=a[j]&~z;w;w&=w-1){int b=w&-w;D[(z|b)*nn+i*n+__builtin_ctz(b)]+=c;}
  }
  h=H.back();opt.e=h;
  for(int z=1;z<N;z++)for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   int c=H[(N-1)^z]*D[z*nn+i*n+j],k=i*n+j;
   if(__builtin_parity(z))s[i].o+=c,t[j].o+=c,M[k].o+=c,opt.o+=c;
   else s[i].e+=c,t[j].e+=c,M[k].e+=c,opt.e+=c;
  }
 }
};
void scan(vector<V> const& parts,int g,Bank&B){
 V a=ordinal(parts);
 int n=a.size(),m=parts.size();
 static map<V,Local>C;
 vector<Local*>L;
 for(auto& T:parts)L.push_back(&(C.try_emplace(T,T).first->second));
 V start(1,0),prefix(1,1);
 for(auto r:L)start.push_back(start.back()+r->n),prefix.push_back(prefix.back()*r->h);
 if(prefix.back()!=g)abort();
 V K(n*n),E(n),RO(n);
 for(int j=0;j<m;j++){
  auto A=L[j];int outer=g/A->h;
  for(int u=0;u<A->n;u++)for(int v=0;v<A->n;v++){
   int i=start[j]+u,k=start[j]+v;
   auto w=A->M[u*A->n+v];K[i*n+k]=w.e*outer;RO[i]+=w.e*outer;E[i]+=w.o*outer;
  }
  P mid={1,0};
  for(int k=j+1;k<m;k++){
   auto Z=L[k];int out=prefix[j]*(g/prefix[k+1]);
   for(int u=0;u<A->n;u++)for(int v=0;v<Z->n;v++){
    auto w=A->s[u]*mid*Z->t[v];
    int i=start[j]+u,q=start[k]+v;K[i*n+q]=w.e*out;RO[i]+=w.e*out;E[i]+=w.o*out;
   }
   mid=mid*Z->opt;
  }
 }
 auto dfs=[&](auto&&self,int pos,int U,int k,int e,int o)->void{
  if(abs(k-g)<=1&&e&&o==e){
   int t=__builtin_ctz(2*e),b=(2*e)>>t,d=2*(g-k),h=2*e-g+d;
   auto& b0=B[b];if(!b0.count(d)||b0[d].t>t)b0[d]={t,h,U,a};
  }
  for(int i=pos;i<n;i++){
   int add=0;
   for(int v=U;v;v&=v-1){int j=__builtin_ctz(v);add+=K[i*n+j]+K[j*n+i];}
   if(k+add<=g+1)self(self,i+1,U|(1<<i),k+add,e+E[i],o+RO[i]);
  }
 };
 dfs(dfs,0,0,0,0,0);
}
void outw(W const&w){
 cout<<" t="<<w.t<<" h="<<w.h<<" U="<<w.U<<" D";
 for(int x:w.D)cout<<" "<<x;
 cout<<endl;
}
int main(int argc,char**argv){
 for(int arg=3;arg<argc;arg++){
  int g=atoi(argv[arg]);if(g<5||g>73||!(g&1))abort();auto C=factors(g,argv[1],argv[2]);
  int cap=4+(31-__builtin_clz(g+1)),seen=0;
  cout<<"g="<<g<<" singletons cap "<<cap<<endl;
  for(auto&[h,S]:C)cout<<"strong "<<h<<" "<<S.size()<<endl;
  vector<V> cs;Bank B;
  auto rec=[&](auto&&self,int rem)->void{
   if(rem!=1){
    for(auto& [h,S]:C)if(rem%h==0)for(auto& a:S){cs.push_back(a);self(self,rem/h);cs.pop_back();}
    return;
   }
   vector<V> tmpl;
   auto locations=[&](auto&&self,int at,int left)->void{
    int old=tmpl.size();
    for(int k=0;k<=left;k++){
     if(k)tmpl.push_back({0});
     if(at==(int)cs.size()){scan(tmpl,g,B);seen++;}
     else{tmpl.push_back(cs[at]);self(self,at+1,left-k);tmpl.pop_back();}
    }
    tmpl.resize(old);
   };
   locations(locations,0,cap);
  };
  rec(rec,g);
  cout<<"shapes "<<seen<<endl;
  for(auto &[a,ds]:B)if(ds.count(0)){
   cout<<"A "<<a;outw(ds.at(0));
   for(auto& [d,w]:ds)if(d){cout<<" d "<<d;outw(w);}
  }
 }
}
