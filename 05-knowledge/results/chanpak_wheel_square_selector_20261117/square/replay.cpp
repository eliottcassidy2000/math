#include<bits/stdc++.h>
using namespace std;
using M=uint64_t;
vector<M> adj;
unsigned a[40],n; M ALL,odd;
map<M,vector<int>> table;
long long countp,keypaths;
set<pair<int,int>> forced;
const M KEY=1073479743;
void visit(int v,M rem,int j,M parity){
 a[j]=v;
 if(!(j&1))parity|=1ull<<v;
 if(!rem){
  countp++;M key=parity^odd;
  if(key==KEY){set<pair<int,int>> edges;for(unsigned t=1;t<n;t++)edges.insert(minmax(a[t-1]+1,a[t]+1));
   if(!keypaths++)forced=edges;else {set<pair<int,int>> cur;set_intersection(edges.begin(),edges.end(),forced.begin(),forced.end(),inserter(cur,cur.end()));forced.swap(cur);}}
  if(!table.count(key))for(int i=0;i<n;i++)table[key].push_back(a[i]+1);
  return;
 }
 M next=rem&adj[v];if(!next)return;
 unsigned leaves=0;
 for(M t=rem;t;t&=t-1){
  int u=__builtin_ctzll(t);
  auto c=__builtin_popcountll(adj[u]&(rem|(1ull<<v)));
  if(c==0)return;
  if(c==1&&++leaves>1)return;
 }
 M seen=0,front=next;
 while(front){auto u=__builtin_ctzll(front);front&=front-1;seen|=1ull<<u;front|=adj[u]&rem&~seen;}
 if(seen!=rem)return;
 for(M t=next;t;t&=t-1){int u=__builtin_ctzll(t);visit(u,rem^(1ull<<u),j+1,parity);}
}
int main(){
 map<M,int> last; set<M> intersection; set<pair<int,int>> common; bool has29=false;
 const int expected[]={10,12,35,52,19,20,349,392,669,4041,17175,12960};
 for(n=25;n<=36;n++){
  ALL=(1ull<<n)-1;odd=0;adj.assign(n,0);table.clear();countp=keypaths=0;forced.clear();
  for(int i=0;i<(int)n;i++){
   if(!(i&1))odd|=1ull<<i;
   for(int j=i+1;j<(int)n;j++){int q=sqrt(i+j+2);if(q*q==i+j+2){adj[i]|=1ull<<j;adj[j]|=1ull<<i;}}
  }
  for(int i=0;i<(int)n;i++)visit(i,ALL^(1ull<<i),0,0);
  map<M,int> cur;int streak=0;M exemplar=0;
  for(auto &[x,w]:table){cur[x]=last.count(x)?last[x]+1:1;if(cur[x]>streak){streak=cur[x];exemplar=x;}}
  cout<<n<<" "<<countp/2<<" "<<table.size()<<" "<<streak<<" "<<exemplar;
  if(streak){cout<<" ";for(int i:table[exemplar])cout<<i<<",";}
  cout<<endl;
  if(countp/2!=expected[n-25])return 2;
  if(n==29)has29=table.count(KEY);
  if(n==30)for(auto &[x,w]:table)intersection.insert(x);
  if(n>30){set<M> z;for(auto &[x,w]:table)if(intersection.count(x))z.insert(x);intersection.swap(z);}
  if(n>=30){
    if(!table.count(KEY))return 3;
    if(n==30)common=forced;else{set<pair<int,int>> z;set_intersection(common.begin(),common.end(),forced.begin(),forced.end(),inserter(z,z.end()));common.swap(z);}
    cout<<"KEY "<<n<<" "<<keypaths<<": ";for(auto x:table[KEY])cout<<x<<",";
    cout<<" forced ";for(auto [x,y]:forced)cout<<x<<"-"<<y<<",";cout<<"\n";
  }
  last.swap(cur);
 }
 if(intersection!=set<M>{KEY} || has29 || common!=set<pair<int,int>>{{7,18},{9,27},{19,30}})return 4;
 cout<<"PASS unique 30–36 parity key\n";
}
