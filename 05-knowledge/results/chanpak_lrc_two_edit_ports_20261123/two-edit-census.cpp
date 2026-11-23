#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>
using namespace std;
using Rows=vector<int>;

bool eligible(int p){
 if(p<5||p==7)return false;
 for(int j=2;j*j<=p;j++)if(p%j==0)return false;
 return true;
}
Rows anchor(int n,int mask){
 Rows r(n);int bit=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)
  if(j==i+1||((mask>>bit++)&1))r[i]|=1<<j;
  else r[j]|=1<<i;
 return r;
}
void flip(Rows&r,int u,int v){r[u]^=1<<v;r[v]^=1<<u;}
bool strong(const Rows&r){
 int s=1,p=1,z=0;
 while(s!=z){z=s;for(int i=0;i<(int)r.size();i++)if(s>>i&1)s|=r[i];}
 z=0;
 while(p!=z){z=p;for(int i=0;i<(int)r.size();i++)if(r[i]&p)p|=1<<i;}
 return s==((1<<r.size())-1)&&p==s;
}
struct DP{
 array<array<int,10>,1024> start;
 array<int,1024> h;
 int compute(const Rows&r){
  const int N=1<<r.size();
  h[0]=1;
  for(int s=1;s<N;s++){
   h[s]=0;
   for(int t=s;t;t&=t-1){
    int i=__builtin_ctz(t),m=s^(1<<i),c=!m;
    for(int j=r[i]&m;j;j&=j-1)c+=start[m][__builtin_ctz(j)];
    start[s][i]=c;h[s]+=c;
   }
  }
  return h[N-1];
 }
};
Rows pad(const Rows&c,int k){
 Rows r=c;int n=c.size(),F=(1<<(n+k))-1;
 for(int i=0;i<n;i++)r[i]|=F^((1<<n)-1);
 for(int j=n;j<n+k;j++)r.push_back(F^((1<<(j+1))-1));
 return r;
}

int main(){
 DP dp,aux;
 const set<int> seven{5,11,13,19};
 const set<int> eight{5,11,13,17,19,23,29,31,37,41,43,47,53,59,67};
 for(int n=4;n<=8;n++){
  vector<pair<int,int>> E;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)E.push_back({i,j});
  set<int> hits;
  long long targetFrames=0,arcs=0;
  const int frames=1<<((n-1)*(n-2)/2);
  for(int mask=0;mask<frames;mask++){
   Rows r=anchor(n,mask);
   int h=dp.compute(r),p=h/7;
   if(h!=7*p||!eligible(p))continue;
   targetFrames++;
   if(hits.count(p))continue; // a witness already establishes presence, never absence
   for(size_t i=0;i<E.size()&&!hits.count(p);i++){
    auto [a,b]=E[i];flip(r,a,b);
    for(size_t j=i+1;j<E.size()&&!hits.count(p);j++){
     auto [c,d]=E[j];flip(r,c,d);arcs++;
     if(aux.compute(r)==p)hits.insert(p);
     flip(r,c,d);
    }
    flip(r,a,b);
   }
  }
  const set<int> empty;
  assert(hits==(n<7?empty:n==7?seven:eight));
  cout<<"order "<<n<<": "<<targetFrames<<" eligible upper frames, "<<arcs<<" pair tests; p =";
  for(int p:hits)cout<<" "<<p;
  cout<<endl;
 }
}
