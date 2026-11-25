#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <array>
#include <algorithm>
using namespace std;
struct Node {array<unsigned char,56> h;array<unsigned char,21> p;};
int n; vector<Node>A;vector<unsigned char>c;
struct Option{int u;unsigned mask;};
vector<int> path,best;
long long visits[8]={};
const Node* anchor=nullptr;
const int U=1<<21;
vector<unsigned> subline;
vector<unsigned char> fano;
array<unsigned,21> through, lines;
bool geometry(){
 if(!anchor||path.size()<2)return true;
 for(int j=0;j<21;j++){
  unsigned M=1u<<j,dup=0;
  auto add=[&](int x){if(!(M&(1u<<x)))M|=1u<<x;else dup|=1u<<x;};
  add(anchor->p[j]);
  for(int k:path)add(A[k].p[j]);
  if(dup){
   if(__builtin_popcount(dup)!=1||!(subline[M^dup]&~through[__builtin_ctz(dup)]))return false;
  } else if(!fano[M]&&!subline[M]){
   bool possible=false;
   for(unsigned U=M;U;U&=U-1)if(subline[M^(U&-U)]){possible=true;break;}
   if(!possible)return false;
  }
 }
 return true;
}
bool dfs(const vector<Option>& avail,unsigned used,int goal) {
 ++visits[path.size()];
 if(path.size()>best.size())best=path;
 if(path.size()==goal)return true;
 for(size_t j=0;j<avail.size();j++){
  if(path.size()+avail.size()-j < goal)break;
  auto [u,m]=avail[j];
  assert(!(m&used));
  unsigned used1=used|m;
  vector<Option> next; next.reserve(avail.size()-j-1);
  for(size_t k=j+1;k<avail.size();k++){
    auto [v,b]=avail[k];int color=c[u*n+v];
    if(!color||(b&used1))continue;
    unsigned z=1u<<(color-1);
    if(z&(b|used1))continue;
    next.push_back({v,b|z});
  }
  path.push_back(u);
  if(geometry()&&dfs(next,used1,goal))return true;
  path.pop_back();
 }
 return false;
}
int main(int argc,char**argv){
 ifstream f(argc>1?argv[1]:"rainbow-data.bin",ios::binary);
 uint32_t a,b; f.read((char*)&a,4);f.read((char*)&b,4);n=a;
 assert(n==7560);vector<uint32_t> reps(b);f.read((char*)reps.data(),4*b);
 A.resize(n);
 auto read=[&](Node& k){f.read((char*)k.h.data(),56);f.read((char*)k.p.data(),21);};
 for(auto& k:A)read(k);
 uint32_t m;f.read((char*)&m,4);vector<Node> zero(m);for(auto& k:zero)read(k);
 f.read((char*)lines.data(),84);assert(f);
 subline.resize(U);fano.resize(U);through.fill(0);
 for(int i=0;i<21;i++){
  unsigned l=lines[i];assert(__builtin_popcount(l)==5);
  for(unsigned U=l;U;U&=U-1)through[__builtin_ctz(U)]|=1u<<i;
  for(unsigned s=l;;s=(s-1)&l){subline[s]|=1u<<i;if(s==0)break;}
 }
 int fcount=0;
 for(unsigned m=127;m<U;){
  bool ok=true;for(unsigned l:lines)if(!(__builtin_popcount(l&m)&1)){ok=false;break;}
  if(ok){fcount++;for(unsigned s=m;;s=(s-1)&m){fano[s]=1;if(!s)break;}}
  unsigned b=m&-m,c=m+b;m=(((c^m)>>2)/b)|c;
 }
 assert(fcount==360);
 c.resize((size_t)n*n);
 vector<unsigned> id(n,0);
 long long edges=0;
 for(int i=0;i<n;i++){
  for(int k=0;k<21;k++)if(A[i].p[k]==k)id[i]=1u<<k;
  for(int j=0;j<i;j++){
   if(id[i]==id[j])continue; // not possible in a rainbow set with identity
   int k=-1;
   for(int w=0;w<21;w++)if(A[i].p[w]==A[j].p[w]){if(k>=0){k=-1;break;}k=w;}
   if(k<0)continue;
   unsigned bit=1u<<k;if(bit&(id[i]|id[j]))continue;
   bool good=true;
   for(int w=0;w<56;w++)if(A[i].h[w]==A[j].h[w]){good=false;break;}
   if(good){c[i*n+j]=c[j*n+i]=k+1;edges++;}
  }
 }
 assert(edges==4999680);
 cout<<"rainbow candidate graph "<<n<<" "<<edges<<endl;
 int goal=6;
 for(auto r:reps){
  vector<Option> opts;
  for(int v=0;v<n;v++)if(int k=c[r*n+v])opts.push_back({v,id[v]|(1u<<(k-1))});
  path={int(r)};
  if(dfs(opts,id[r],goal)){
    cout<<"witness";for(int j:path)cout<<" "<<j;cout<<endl;
    return 1;
  }
 }
 cout<<"no size 7 all colored; visits";for(auto b:visits)cout<<" "<<b;cout<<endl;
 fill(begin(visits),end(visits),0);best.clear();
 for(auto& r:zero){
  anchor=&r;
  vector<Option> opts;
  for(int i=0;i<n;i++){
   int k=-1;
   for(int x=0;x<21;x++)if(r.p[x]==A[i].p[x]){if(k>=0){k=-1;break;}k=x;}
   if(k<0)continue;
   unsigned bit=1u<<k;if(bit&id[i])continue;
   bool good=true;for(int w=0;w<56;w++)if(r.h[w]==A[i].h[w]){good=false;break;}
   if(good)opts.push_back({i,id[i]|bit});
  }
  path.clear();
  if(dfs(opts,0,5)){
    cout<<"one missing color witness";for(int j:path)cout<<" "<<j;cout<<endl;return 1;
  }
 }
 cout<<"no size 7 with one uncolored edge and the pencil/Fano conditions; visits";for(auto b:visits)cout<<" "<<b;cout<<endl;
}