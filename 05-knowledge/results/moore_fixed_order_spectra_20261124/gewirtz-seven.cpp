#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
using u32=uint32_t;using u64=uint64_t;
const int Z=1<<21;
struct Perm {array<unsigned char,56> h;array<unsigned char,42> pl;};
struct State {u32 used;unsigned char dup,c;}; // dup ignored when c=0
using Image=array<State,42>;
vector<Perm>P;
vector<u64>adj,spot;
bool no_coinc=false;
vector<u32> next0(Z),Lunion(Z),Fnext(Z);
array<u32,21> lines;
int N,W;
u64 visited[8]={}, counter=0;
vector<int> seq;
bool found=false;
bool invert(int i,int j){
 for(int k=0;k<56;k+=8){
  u64 u,v;
  memcpy(&u,P[i].h.data()+k,8);memcpy(&v,P[j].h.data()+k,8);
  u^=v;
  if((u-0x0101010101010101ULL)&~u&0x8080808080808080ULL)return false;
 }
 return true;
}
bool compat(int i,int j){
 if(!((adj[size_t(i)*W+j/64]>>(j%64))&1))return false;
 if(no_coinc){for(int k=0;k<42;k++)if(P[i].pl[k]==P[j].pl[k])return false;}
 return true;
}
const u64* spots(int k,int j,int m){return &spot[((size_t(k)*3+j)*128+m)*W];}
bool check_geom(const array<u32,42>& ok,int v){
 for(int k=0;k<42;k++)if(!(ok[k]&(1u<<P[v].pl[k])))return false;
 return true;
}
u32 allowed(State a){
 if(no_coinc){assert(!a.c);return Fnext[a.used];}
 if(!a.c)return next0[a.used];
 const u32 bit=1u<<a.dup,s=a.used^bit;
 if(a.c==2)return Lunion[s]&~s;
 return Lunion[a.used]&~a.used;
}
void put(State& a,int x){
 u32 bit=1u<<x;
 if(!(a.used&bit)){a.used|=bit;return;}
 if(a.c){assert(a.dup==x&&a.c==2);a.c=3;}
 else {a.dup=x;a.c=2;}
}
bool search(const Image& I,const vector<int>& opt){
 ++visited[seq.size()];
 if(seq.size()==6){
 for(int k=0;k<42;k++){
   for(u32 l:lines){
     int c=(l>>(k%21))&1;
     for(int i:seq)c+=(l>>P[i].pl[k])&1;
     assert(c&1);
   }
 }
 found=true;++counter;cout<<"compatible seven-set:";for(int i:seq)cout<<" "<<i;cout<<'\n';return false;}
 if(opt.size()<6-seq.size())return false;
 vector<u64> optbits(W);
 if(seq.size()>=2&&opt.size()>256)
   for(int x:opt)optbits[x>>6]|=u64(1)<<(x&63);
 for(size_t j=0;j+6-seq.size()<=opt.size();j++){
  const int u=opt[j];
  Image J=I;
  array<u32,42> ok;
  bool possible=true;
  for(int k=0;k<42;k++){
   put(J[k],P[u].pl[k]);
   ok[k]=allowed(J[k]);
   if(seq.size()<5&&!ok[k]){possible=false;break;}
  }
  if(!possible)continue;
  seq.push_back(u);
  if(seq.size()==6){if(search(J,{}))return true;}
  else {
   vector<int>next;next.reserve(opt.size()-j-1);
   if(seq.size()>=3&&opt.size()>256){
     array<int,42> order;for(int k=0;k<42;k++)order[k]=k;
     partial_sort(order.begin(),order.begin()+8,order.end(),[&](int a,int b){
       return __builtin_popcount(ok[a])<__builtin_popcount(ok[b]);
     });
     array<array<const u64*,3>,8> A;
     for(int k=0;k<8;k++)for(int z=0;z<3;z++)
       A[k][z]=spots(order[k],z,(ok[order[k]]>>(z*7))&127);
     for(int w=u>>6;w<W;w++){
       u64 S=optbits[w]&adj[size_t(u)*W+w];
       if(w==(u>>6)){
         const int b=(u&63)+1;
         S&=(b==64?0:~u64(0)<<b);
       }
       if(!S)continue;
       for(int k=0;k<8&&S;k++)S&=A[k][0][w]|A[k][1][w]|A[k][2][w];
       while(S){
         const int v=w*64+__builtin_ctzll(S);S&=S-1;
         if((!no_coinc||compat(u,v))&&check_geom(ok,v))next.push_back(v);
       }
     }
   }else{
     for(size_t x=j+1;x<opt.size();x++){
       const int v=opt[x];
       if(compat(u,v)&&check_geom(ok,v))next.push_back(v);
     }
   }
   if(search(J,next))return true;
  }
  seq.pop_back();
 }
 return false;
}
int main(int argc,char** argv){
 ifstream f(argc>1?argv[1]:"gewirtz-seven.bin",ios::binary);
 u32 n,m;f.read((char*)&n,4);f.read((char*)&m,4);assert(n==19080&&m==6);
 N=n;W=(N+63)/64;
 vector<u32> reps(m);
 f.read((char*)reps.data(),m*4);f.read((char*)lines.data(),84);
 P.resize(n);
 for(auto& r:P){f.read((char*)r.h.data(),56);f.read((char*)r.pl.data(),42);}assert(f);
 for(u32 l:lines){
  for(u32 s=l;;s=(s-1)&l){Lunion[s]|=l;if(!s)break;}
  for(int x=0;x<21;x++){
   const u32 bit=1u<<x,T=l|bit;
   for(u32 s=T;;s=(s-1)&T){next0[s]|=(T^s)|bit;if(!s)break;}
  }
 }
 int fcount=0;
 for(u32 m=127;m<Z;){
  bool ok=true;for(u32 l:lines)if(!(__builtin_popcount(l&m)&1)){ok=false;break;}
  if(ok){fcount++;for(u32 s=m;;s=(s-1)&m){next0[s]|=m^s;Fnext[s]|=m^s;if(!s)break;}}
  u32 b=m&-m,c=m+b;m=(((c^m)>>2)/b)|c;
 }
 assert(fcount==360);
 adj.resize(size_t(N)*W);
 u64 edges=0;
 for(int i=0;i<N;i++)for(int j=0;j<i;j++)if(invert(i,j)){
   adj[size_t(i)*W+j/64]|=u64(1)<<(j%64);adj[size_t(j)*W+i/64]|=u64(1)<<(i%64);edges++;
 }
 cout<<"derangement compatibility edges "<<edges<<endl;
 spot.resize(size_t(42)*3*128*W);
 for(int i=0;i<N;i++)for(int k=0;k<42;k++){
   const int v=P[i].pl[k];
   const_cast<u64*>(spots(k,v/7,1<<(v%7)))[i/64]|=u64(1)<<(i%64);
 }
 for(int k=0;k<42;k++)for(int j=0;j<3;j++)for(int m=1;m<128;m++)if(m&(m-1)){
   auto A=spots(k,j,m&-m),B=spots(k,j,m&(m-1));auto T=const_cast<u64*>(spots(k,j,m));
   for(int i=0;i<W;i++)T[i]=A[i]|B[i];
 }
 for(size_t ri=0;ri<reps.size();ri++){
   const u32 r=reps[ri];
   if(argc>2&&size_t(stoi(argv[2]))!=ri)continue;
   // If no pair fixed one point in the already exhausted mixed search,
   // every pair must be coincidence-free.
   int fixed=0;for(int k=0;k<21;k++)fixed+=P[r].pl[k]==k;
   no_coinc=(fixed==0);
   Image I{};
   array<u32,42> ok;
   for(int k=0;k<42;k++){
     I[k].used=1u<<(k%21);
     put(I[k],P[r].pl[k]);
     ok[k]=allowed(I[k]);
   }
   vector<int>opts;
   for(int v=0;v<N;v++)if(compat(r,v)){
     if(no_coinc){bool good=true;for(int k=0;k<42;k++)if(P[v].pl[k]==k%21){good=false;break;}if(!good)continue;}
     bool good=true;for(int k=0;k<42;k++)if(!(ok[k]&(1u<<P[v].pl[k]))){good=false;break;}
     if(good)opts.push_back(v);
   }
   seq={int(r)};
   cout<<"class "<<r<<" first pool "<<opts.size()<<endl;
   if(search(I,opts))break;
   cout<<" cumulative:";for(auto b:visited)cout<<" "<<b;cout<<endl;
 }
 cout<<"normalized seven-set leaves "<<counter<<endl;
}