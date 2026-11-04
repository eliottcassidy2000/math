// Independent audit: permutation orbits, literal split words, full attachments.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using I=long long;
using T=vector<unsigned>;
void need(bool ok,const char* message){if(!ok){cerr<<message<<endl;exit(2);}}

T decode(unsigned code,int n){
 T t(n);int k=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)
  if(code>>k&1)t[i]|=1u<<j;else t[j]|=1u<<i;
 return t;
}
unsigned encode(const T&t,const vector<int>&p){
 unsigned code=0;int k=0;
 for(int i=0;i<(int)t.size();i++)for(int j=i+1;j<(int)t.size();j++,k++)
  if(t[p[i]]>>p[j]&1)code|=1u<<k;
 return code;
}
vector<T> orbit_bank(int n){
 vector<bool>seen(1u<<(n*(n-1)/2));vector<T>out;vector<int>p(n);
 for(unsigned code=0;code<seen.size();code++)if(!seen[code]){
  T t=decode(code,n);out.push_back(t);iota(p.begin(),p.end(),0);
  do{seen[encode(t,p)]=true;}while(next_permutation(p.begin(),p.end()));
 }
 return out;
}
struct Cover{I h=0,p=0;vector<I>f;vector<vector<I>>j;};
Cover literal_covers(const T&t){
 int n=t.size();Cover q;q.f.assign(n,0);q.j.assign(n,vector<I>(n));
 vector<int>w(n);iota(w.begin(),w.end(),0);
 do{
  int bad=-1,count=0;
  for(int i=0;i+1<n;i++)if(!(t[w[i]]>>w[i+1]&1)){
   bad=i;if(++count>1)break;
  }
  if(count>1)continue;
  if(count==0){
   q.h++;q.p+=n+1;
   for(int k=1;k<=n;k++)q.f[w[k-1]]++;
   for(int k=1;k<n;k++)q.j[w[k-1]][w[k]]++;
  }else{
   q.p++;q.f[w[bad]]++;q.j[w[bad]][w[bad+1]]++;
  }
 }while(next_permutation(w.begin(),w.end()));
 return q;
}
T attach(const T&b,unsigned W,int clock){
 int n=b.size(),root=n+clock;T t(root+1);
 for(int i=0;i<clock;i++){
  for(int j=i+1;j<root;j++)t[i]|=1u<<j;
  t[root]|=1u<<i;
 }
 for(int i=0;i<n;i++){
  t[clock+i]=b[i]<<clock;
  if(W>>i&1)t[clock+i]|=1u<<root;else t[root]|=1u<<(clock+i);
 }
 return t;
}
bool strong(const T&t){
 auto reach=[&](bool reverse){
  unsigned have=1,todo=1;
  while(todo){
   unsigned bit=todo&-todo;todo-=bit;int v=countr_zero(bit);unsigned next=0;
   if(!reverse)next=t[v];
   else for(int u=0;u<(int)t.size();u++)if(t[u]>>v&1)next|=1u<<u;
   next&=~have;have|=next;todo|=next;
  }
  return have==((1u<<t.size())-1);
 };
 return reach(false)&&reach(true);
}
I full_h(const T&t){
 int n=t.size(),N=1<<n;vector<I>d(N*n);
 for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;
 for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(d[m*n+v]){
  unsigned can=t[v]&~m;
  while(can){
   unsigned bit=can&-can;can-=bit;
   d[(m|bit)*n+countr_zero(bit)]+=d[m*n+v];
  }
 }
 I h=0;for(int v=0;v<n;v++)h+=d[(N-1)*n+v];return h;
}
T born(T t,unsigned cut){
 int n=t.size();for(int i=0;i<n;i++)if(cut>>i&1)t[i]|=1u<<n;
 t.push_back(((1u<<n)-1)^cut);return t;
}
struct Stats{I graphs=0,portsets=0,strongports=0,f63=0;array<I,4>hits{};};
Stats scan(const vector<T>&bank){
 Stats s;array<I,4>goals{-37,-35,-31,-29};
 for(const auto&t:bank){
  s.graphs++;auto q=literal_covers(t);
  for(unsigned W=0;W<(1u<<t.size());W++){
   s.portsets++;auto full=attach(t,W,1);if(!strong(full))continue;
   s.strongports++;I f=0,j=0;
   for(int u=0;u<(int)t.size();u++)if(W>>u&1){
    f+=q.f[u];
    for(int v=0;v<(int)t.size();v++)if(W>>v&1)j+=q.j[u][v];
   }
   if(f!=63)continue;s.f63++;
   I H1=full_h(full),H2=full_h(attach(t,W,2));
   need(H2-H1==126,"literal slope disagrees");
   need(2*H1-H2==q.h-j,"literal intercept disagrees");
   for(int k=0;k<4;k++)if(q.h-j==goals[k])s.hits[k]++;
  }
 }
 return s;
}
void print(int n,const Stats&s){
 cout<<"order="<<n<<" graphs="<<s.graphs<<" all_portsets="<<s.portsets
     <<" strong_portsets="<<s.strongports<<" f63="<<s.f63<<" hits22_24_42_44=";
 for(auto x:s.hits)cout<<x<<",";cout<<"\n";
}
int main(){
 array<int,7>expected{0,1,1,2,4,12,56};vector<T>six;I f63=0;array<I,4>total{};
 for(int n=1;n<=6;n++){
  auto bank=orbit_bank(n);need(bank.size()==expected[n],"independent class inventory");
  if(n==6)six=bank;auto s=scan(bank);print(n,s);f63+=s.f63;
  for(int k=0;k<4;k++)total[k]+=s.hits[k];
 }
 vector<T>seven;
 for(auto&t:six)for(unsigned cut=0;cut<64;cut++)seven.push_back(born(t,cut));
 need(seven.size()==3584,"birth coverage");auto s=scan(seven);print(7,s);f63+=s.f63;
 for(int k=0;k<4;k++)total[k]+=s.hits[k];
 need(total[0]>0&&total[2]>0,"positive22/42 control");
 need(total[1]==0&&total[3]==0,"unexpected24/44 hit");
 need(f63==185,"fixed-slope universe changed");
 cout<<"PASS independent orbit-marking, literal permutation-split covers, actual attachment strongness, and two full Hamilton counts for all "<<f63<<" fixed-slope candidates\n";
}
