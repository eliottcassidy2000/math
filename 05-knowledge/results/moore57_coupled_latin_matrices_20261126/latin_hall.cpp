
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using U=uint64_t;
constexpr int D=56,NMAX=3136,AX=4;
using Bits=bitset<NMAX>;
int m,n;
U S;
vector<int> adj[NMAX],cand[NMAX];
Bits chosen[NMAX],possible[NMAX];
array<int,4> tag[NMAX];
array<U,4> missing[NMAX];
long long checks[12]{},drops[12]{},counts[5]{};
int forced=0;
U bit(int x){return 1ull<<x;}
int pop(U x){return __builtin_popcountll(x);}
template<class F> void each(U b,F cb){while(b){int x=__builtin_ctzll(b);b&=b-1;cb(x);}}
bool support(array<U,D>& allowed,U I,U J){
 assert(pop(I)==pop(J));
 int mate[D];fill(mate,mate+D,-1);
 vector<int> ix;each(I,[&](int x){ix.push_back(x);});
 sort(ix.begin(),ix.end(),[&](int x,int y){return pop(allowed[x])<pop(allowed[y]);});
 auto augment=[&](auto&& A,int i,U& seen)->bool{
  U choices=allowed[i]&~seen;
  while(choices){
   int j=__builtin_ctzll(choices);choices&=choices-1;
   if(seen&bit(j))continue;
   seen|=bit(j);
   if(mate[j]<0||A(A,mate[j],seen)){mate[j]=i;return true;}
  }return false;
 };
 for(int i:ix){U seen=0;if(!augment(augment,i,seen))return false;}
 array<U,D> reach{};
 for(int i:ix)each(allowed[i],[&](int j){assert(mate[j]>=0);reach[i]|=bit(mate[j]);});
 for(int j:ix)for(int i:ix)if(reach[i]&bit(j))reach[i]|=reach[j];
 for(int i:ix){
  U old=allowed[i];
  each(old,[&](int j){if(!(reach[mate[j]]&bit(i)))allowed[i]&=~bit(j);});
 }
 return true;
}
bool drop(int x,int y,int k){
 if(!possible[x][y])return false;
 possible[x][y]=possible[y][x]=false;
 drops[k]++;
 return true;
}
bool local(){
 int p=0;
 for(int a=0;a<AX;a++)for(int b=a+1;b<AX;b++,p++)
 for(int x=0;x<n;x++)if(missing[x][a]){
  array<U,D> tab{};
  for(int y:cand[x])if(possible[x][y])tab[tag[y][a]]|=bit(tag[y][b]);
  checks[p]++;
  if(!support(tab,missing[x][a],missing[x][b])){
   cout<<"vertex_hall_failure "<<x<<" "<<tag[x][0]<<" "<<tag[x][1]<<" "<<a<<" "<<b<<endl;return false;}
  for(int y:cand[x])if(possible[x][y]&&!(tab[tag[y][a]]&bit(tag[y][b])))drop(x,y,p);
 }
 return true;
}
int vert(int side,int a,int b){return side==0?a*m+b:b*m+a;}
bool line(){
 for(int side=0;side<2;side++)for(int a=0;a<m;a++)for(int b=a+1;b<m;b++){
  U I=0,J=0;array<U,D> tab{};
  for(int i=0;i<m;i++){
   const int x=vert(side,a,i);
   if(missing[x][side]&bit(b)){
    I|=bit(i);
    for(int y:cand[x])if(tag[y][side]==b&&possible[x][y])tab[i]|=bit(tag[y][1-side]);
   }
   const int z=vert(side,b,i);
   if(missing[z][side]&bit(a))J|=bit(i);
  }
  if(!I)continue;
  checks[6+side]++;
  if(!support(tab,I,J)){cout<<"line_hall_failure "<<side<<" "<<a<<" "<<b<<endl;return false;}
  for(int i=0;i<m;i++)if(I&bit(i)){
   const int x=vert(side,a,i);
   for(int y:cand[x])if(tag[y][side]==b&&possible[x][y]&&!(tab[i]&bit(tag[y][1-side])))drop(x,y,6+side);
  }
 }
 return true;
}

int ext[NMAX]; U extras;
bool fifth(){
 for(int x=0;x<n;x++)if(x/m!=1&&missing[x][0]){
   U J=extras&~bit(x%m);
   for(int y:adj[x])if(ext[y]>=0)J&=~bit(ext[y]);
   const int capacity=pop(missing[x][0])-pop(J);
   if(capacity<0)return false;
   U B=0;for(int i=0;i<m&&pop(B)<capacity;i++)if(!(J&bit(i)))B|=bit(i);
   assert(pop(B)==capacity);
   for(int a=0;a<4;a++){
    array<U,D> tab{};
    for(int y:cand[x])if(possible[x][y]){
      const auto k=ext[y]; if(k>=0)assert(J&bit(k));
      tab[tag[y][a]]|=k<0?B:bit(k);
    }
    checks[8+a]++;
    if(!support(tab,missing[x][a],J|B)){cout<<"partial_latin_failure "<<x<<" "<<a<<endl;return false;}
    for(int y:cand[x])if(possible[x][y]&&!(tab[tag[y][a]]&(ext[y]<0?B:bit(ext[y]))))drop(x,y,8+a);
   }
 }
 return true;
}
int main(int argc,char**argv){
 int e;cin>>m>>e;n=m*m;assert(m<=D);S=(1ull<<m)-1;
 for(int i=0;i<e;i++){int x,y;cin>>x>>y;assert(!chosen[x][y]);chosen[x][y]=chosen[y][x]=true;adj[x].push_back(y);adj[y].push_back(x);}
 int a=0,b=28;
 if(argc>2){a=stoi(argv[1]);b=stoi(argv[2]);}
 array<int,2> anch{a,b};int prim[2][D][D]{};
 for(int c=0;c<D;c++)for(int x=0;x<D;x++)for(int j=0;j<2;j++)prim[j][c][x]=-1;
 int full=0;
 for(int x=0;x<n;x++){
  tag[x][0]=x/m;tag[x][1]=x%m;
  full+=adj[x].size()==m-1;
  for(int k=0;k<2;k++){
   if(tag[x][0]==anch[k])tag[x][2+k]=tag[x][1];
   else{
    int z=-1;for(int y:adj[x])if(y/m==anch[k]){assert(z<0);z=y%m;}
    assert(z>=0);tag[x][2+k]=z;
   }
   prim[k][tag[x][0]][tag[x][1]]=tag[x][2+k];
  }
 }
 for(int k=0;k<2;k++)for(int i=0;i<m;i++){
  U r=0,c=0;
  for(int j=0;j<m;j++){r|=bit(prim[k][i][j]);c|=bit(prim[k][j][i]);}
  assert(r==S&&c==S);
 }
 int gram[D][D]{},h[3]{};
 for(int x=0;x<n;x++)gram[tag[x][2]][tag[x][3]]++;
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){
  auto z=gram[i][j];assert(z==int(i!=j)+int(prim[1][a][i]==j));h[z]++;
 }
 cout<<"known "<<e<<" full vertices "<<full<<" superposition "<<h[0]<<","<<h[1]<<","<<h[2]<<endl;
 for(int x=0;x<n;x++){
  if(adj[x].size()==m-1){missing[x].fill(0);continue;}
  missing[x].fill(S);
  for(int k=0;k<4;k++)missing[x][k]&=~bit(x%m);
  missing[x][0]=S^bit(x/m);
  for(int y:adj[x])for(int k=0;k<4;k++){
   assert(missing[x][k]&bit(tag[y][k]));missing[x][k]&=~bit(tag[y][k]);
  }
 }
 fill(ext,ext+n,-1);extras=0;
 for(int t=0;t<m;t++)if(adj[m+t].size()==m-1){
   extras|=bit(t);int z=m+t;
   assert(ext[z]<0);ext[z]=t;
   for(int y:adj[z]){assert(ext[y]<0);ext[y]=t;}
 }
 cout<<"partial row 1 anchors "<<pop(extras)<<endl;
 vector<Bits> p2(n),p3(n);
 for(int x=0;x<n;x++){
  p2[x]=chosen[x];p2[x].set(x);
  for(int y:adj[x])p2[x]|=chosen[y];
 }
 for(int x=0;x<n;x++){
  p3[x]=p2[x];for(int y:adj[x])p3[x]|=p2[y];
 }
 for(int x=0;x<n;x++)for(int y=x+1;y<n;y++){
  if(!(missing[x][0]&bit(tag[y][0]))||!(missing[y][0]&bit(tag[x][0])))continue;
  if(!(missing[x][1]&bit(tag[y][1]))||!(missing[y][1]&bit(tag[x][1])))continue;
  counts[0]++;
  if(p3[x][y])continue;
  counts[1]++;
  bool ok=true;
  for(int k=2;k<4;k++)ok &= bool(missing[x][k]&bit(tag[y][k]))&&bool(missing[y][k]&bit(tag[x][k]));
  if(!ok)continue;
  counts[2]++;possible[x][y]=possible[y][x]=true;cand[x].push_back(y);cand[y].push_back(x);
 }
 cout<<"candidate staged "<<counts[0]<<" "<<counts[1]<<" "<<counts[2]<<endl;
 for(int i=0;i<8;i++){
  auto s=accumulate(drops,drops+12,0ll);
  if(!local()||!line()||!fifth())return 0;
  if(s==accumulate(drops,drops+12,0ll))break;
 }
 cout<<"checks";for(long long z:checks)cout<<" "<<z;cout<<endl;
 cout<<"drops";for(long long z:drops)cout<<" "<<z;cout<<endl;
 for(int x=0;x<n;x++)for(int k=0;k<4;k++)each(missing[x][k],[&](int i){
  int t=0;for(int y:cand[x])if(possible[x][y]&&tag[y][k]==i)t++;
  if(t==1)forced++;
 });
 cout<<"remaining "<<counts[2]-accumulate(drops,drops+12,0ll)<<" forced slots "<<forced<<endl;
 if(argc>3){
   int x=stoi(argv[3]);
   cerr<<x<<" "<<possible[x].count()<<endl;
   for(int y:cand[x])if(possible[x][y])cerr<<y<<" "<<tag[y][0]<<" "<<tag[y][1]<<" "<<tag[y][2]<<" "<<tag[y][3]<<endl;
 }

}
