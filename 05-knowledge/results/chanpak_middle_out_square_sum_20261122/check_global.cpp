#include<bits/stdc++.h>
using namespace std;
using I=long long;using A=vector<int>;
struct Cut{
 I h; vector<vector<I>> w;
 Cut(A const&a){
  int n=a.size(),M=1<<n,U=M-1;
  vector<I> hs(M),D(M*n);
  for(int j=0;j<n;j++)D[(1<<j)*n+j]=1;
  for(int S=1;S<M;S++){
   for(int j=0;j<n;j++){
    I p=D[S*n+j];if(!p)continue;
    hs[S]+=p;
    for(int t=a[j]&~S;t;t&=t-1)D[(S|(t&-t))*n+__builtin_ctz(t)]+=p;
   }
  }
  h=hs[U];hs[0]=1;
  w.assign(n,vector<I>(n,0));
  for(int u=0;u<n;u++){
   fill(D.begin(),D.end(),0);
   int S0=1<<u;D[S0*n+u]=1;
   for(int S=S0;S<M;S=(S+1)|S0){
    bool even=__builtin_parity(S)==0;
    for(int j=0;j<n;j++){
     I p=D[S*n+j];if(!p)continue;
     if(even)w[u][j]+=2*p*hs[U^S];
     for(int t=a[j]&~S;t;t&=t-1)D[(S|(t&-t))*n+__builtin_ctz(t)]+=p;
    }
   }
  }
 }
 vector<I> values()const{
  int n=w.size(),M=1<<n;
  vector<I> r(M,h);
  for(int S=1;S<M;S++){
   int j=__builtin_ctz(S),X=S&(S-1);
   I d=0;
   for(int t=0;t<n;t++){
    if((S>>t&1)==0)d+=w[j][t];
    if(X>>t&1)d-=w[t][j];
   }
   r[S]=r[X]+d;
  }
  return r;
 }
};
int n,p,i0,jflip,K;
I targ;
A base;
I visits[12]{}, prunes[12]{};
vector<int> viable;
vector<int> addr;
I hc(A const&a){
 int n=a.size(),M=1<<n;static I dp[65536][16];
 memset(dp,0,(M)*sizeof dp[0]);for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int S=0;S<M-1;S++)for(int i=0;i<n;i++)if(dp[S][i]){
  for(int e=a[i]&~S;e;e&=e-1){int j=__builtin_ctz(e);dp[S|1<<j][j]+=dp[S][i];}
 }
 return accumulate(dp[M-1],dp[M-1]+n,0LL);
}
A flip(A a){a[i0]^=1<<jflip;a[jflip]^=1<<i0;return a;}
A add(A a,int m){int k=a.size();for(int i=0;i<k;i++)if(!(m>>i&1))a[i]|=1<<k;a.push_back(m);return a;}
void search(A a,int frommask){
 int depth=a.size()-n;
 visits[depth]++;
 if(depth==K-1){
  auto C=Cut(a).values();
  bool found=false;
  for(int m:viable)if(m>=frommask)for(int hi=0;hi<(1<<depth);hi++){
   int x=m+(hi<<n);
   if(C[x]!=targ-1)continue;
   found=true;
   A aa=add(a,x);
   if(hc(flip(aa))!=targ+1)continue;
   cout<<"HIT depth="<<depth+1<<" rows";for(int x:aa)cout<<" "<<x;
   cout<<" raw addr";for(int x:addr)cout<<" "<<x;cout<<" "<<x<<endl;
   exit(0);
  }
  prunes[depth]+=!found;return;
 }
 for(int m:viable)if(m>=frommask)for(int hi=0;hi<(1<<depth);hi++){
  int x=m+(hi<<n);A aa=add(a,x);
  if(hc(aa)>targ-1||hc(flip(aa))>targ+1){prunes[depth]++;continue;}
  addr.push_back(x);search(aa,m);addr.pop_back();
 }
}

bool strong(A const&a){
 int n=a.size(),all=(1<<n)-1;
 for(int i=0;i<n;i++){
  int seen=1<<i;
  while(true){
   int z=seen;
   for(int j=0;j<n;j++)if(seen>>j&1)z|=a[j];
   if(z==seen)break;
   seen=z;
  }
  if(seen!=all)return false;
 }
 return true;
}
A canon(A const&a,int x,int y){
 int n=a.size(),all=(1<<n)-1;
 A best(n+2,INT_MAX),path;
 function<void(int)> go=[&](int rest){
  if(!rest){
   A out(n+2,0);
   for(int i=0;i<n;i++){
    if(path[i]==x)out[n]=i;
    if(path[i]==y)out[n+1]=i;
    for(int j=0;j<n;j++)if(a[path[i]]>>path[j]&1)out[i]|=1<<j;
   }
   if(out[n]>out[n+1])swap(out[n],out[n+1]);
   best=min(best,out);
   return;
  }
  int choices=path.empty()?rest:(a[path.back()]&rest);
  for(int q=choices;q;q&=q-1){
   int k=q&-q;path.push_back(__builtin_ctz(k));
   go(rest^k);path.pop_back();
  }
 };
 go(all);return best;
}
set<A> marked_seeds(int n){
 vector<pair<int,int>> free;
 for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)free.emplace_back(i,j);
 A a(n),b(n);set<A> seen;
 for(int m=0;m<(1<<(int)free.size());m++){
  fill(a.begin(),a.end(),0);
  for(int i=0;i+1<n;i++)a[i]|=1<<(i+1);
  for(int k=0;k<(int)free.size();k++){
   auto [i,j]=free[k];
   if(m>>k&1)a[i]|=1<<j;else a[j]|=1<<i;
  }
  if(!strong(a)||hc(a)!=27)continue;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   b=a;b[i]^=1<<j;b[j]^=1<<i;
   if(hc(b)==29&&strong(b))seen.insert(canon(a,i,j));
  }
 }
 return seen;
}
void run(){
 viable.clear();addr.clear();
 fill(begin(visits),end(visits),0);fill(begin(prunes),end(prunes),0);
 auto C=Cut(base).values(),D=Cut(flip(base)).values();
 for(int m=0;m<(1<<n);m++)if(C[m]<=targ-1&&D[m]<=targ+1)viable.push_back(m);
 search(base,0);
 cout<<"MISS n="<<n<<" extras="<<K<<" eligible "<<viable.size();
 for(int i=0;i<K;i++)cout<<" "<<visits[i]<<"/"<<prunes[i];
 cout<<"\n";
}
int main(){
 p=3;targ=120;
 for(n=3;n<=7;n++){
  auto all=marked_seeds(n);
  cout<<"MARKED_SOURCE_CLASSES "<<n<<" "<<all.size()<<"\n";
  for(A a:all){
   i0=a[n];jflip=a[n+1];a.resize(n);
   cout<<"SEED "<<i0<<" "<<jflip<<" rows";
   for(int r:a)cout<<" "<<r;
   cout<<"\n";
   for(K=1;K<=9-n;K++){
    base=a;cout<<"ORDER_PRESERVING ";run();
    base=flip(a);cout<<"ORDER_SWAPPING ";run();
   }
  }
 }
}
