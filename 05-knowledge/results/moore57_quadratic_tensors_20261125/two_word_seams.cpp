// Exact 2-SAT on the union of two complete projected words.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using U=uint64_t;
constexpr U ALL=(1ULL<<52)-1;
struct Word {
 array<int,52> y;
 array<array<int,3>,52> a;
 U routes[3][56]{};
 int drop;
};
struct Seam {
 const Word *w[2];
 array<int,104> idx{},low{},comp{};
 int stack[104],top=0, cc=0,t=0; bool on[104]{};
 U conf(int v,int k) const {
  int i=v/2,j=v%2;
  return (w[k]->routes[0][w[j]->a[i][0]]|w[k]->routes[1][w[j]->a[i][1]]|
          w[k]->routes[2][w[j]->a[i][2]])&~(1ULL<<i);
 }
 bool propagate(U &a,U &b,U A,U B) const {
  while(A||B){
   if(A&B)return false;
   A&=~a;B&=~b;
   if((A&b)||(B&a))return false;
   a|=A;b|=B;
   U X=0,Y=0;
   for(U R=A;R;R&=R-1){
    int v=2*__builtin_ctzll(R); X|=conf(v,1); Y|=conf(v,0);
   }
   for(U R=B;R;R&=R-1){
    int v=2*__builtin_ctzll(R)+1; X|=conf(v,1); Y|=conf(v,0);
   }
   A=X&~a;B=Y&~b;
  }
  return true;
 }
 void dfs(int v,U free) {
  idx[v]=low[v]=++t;stack[top++]=v;on[v]=1;
  for(int j=0;j<2;j++)for(U R=conf(v,j)&free;R;R&=R-1){
   int m=2*__builtin_ctzll(R)+(j^1);
   if(!idx[m]){dfs(m,free);low[v]=min(low[v],low[m]);}
   else if(on[m])low[v]=min(low[v],idx[m]);
  }
  if(low[v]==idx[v]){
   ++cc;
   int m;do{m=stack[--top];on[m]=0;comp[m]=cc;}while(m!=v);
  }
 }
 bool solve(U equal,U& A,U& B){
  if(!propagate(A,B,equal,0))return false;
  U F=ALL&~(A|B);
  for(U h=F;h;h&=h-1){
   int i=2*__builtin_ctzll(h);
   if(!idx[i])dfs(i,F);
   if(!idx[i+1])dfs(i+1,F);
   if(comp[i]==comp[i+1])return false;
  }
  for(U h=F;h;h&=h-1){
   int i=2*__builtin_ctzll(h);
   (comp[i]<comp[i+1]?A:B)|=h&-h;
  }
  return true;
 }
};
int main(int argc,char**argv){
 int n;cin>>n;
 if(n!=53||argc<2)return 2;
 array<int,3>C[53][56];bool valid[53][56]{};
 for(int i=0;i<n;i++){
  int old,k;cin>>old>>k;
  while(k--){int y;cin>>y;for(int &z:C[i][y])cin>>z;valid[i][y]=true;}
 }
 ifstream f(argv[1]);
 vector<Word> W;
 int g,k,y;
 while(f>>g>>k>>y){
  if(y!=7||g<0||g>=3)return 3;
  Word w{};w.drop=g;
  for(int i=0;i<52;i++){
   f>>y;if(y<0||y>=56||!valid[i+1][y])return 4;
   w.y[i]=y;w.a[i]=C[i+1][y];
   for(int j=0;j<3;j++){
    int z=w.a[i][j];if(z==C[0][7][j])return 4;
    w.routes[j][z]|=1ULL<<i;
   }
  }
  W.push_back(w);
 }
 if(W.size()!=5223)return 5;
 uint64_t tested[3][3]{};
 for(int a=0;a<(int)W.size();a++)for(int b=a+1;b<(int)W.size();b++){
  auto &x=W[a], &y=W[b];
  tested[x.drop][y.drop]++;
  U eq=0;
  for(int i=0;i<52;i++)if(x.y[i]==y.y[i])eq|=1ULL<<i;
  // Treat an identical physical contact as forced to its first copy.
  Seam S{{&x,&y}};U A=0,B=0;
  if(!S.solve(eq,A,B))continue;
  cout<<"full splice "<<a<<" "<<b<<endl;
  cout<<"7";
  U seen[3]{};
  for(int i=0;i<52;i++){
   auto &w=(A>>i&1)?x:y;
   cout<<" "<<w.y[i];
   for(int j=0;j<3;j++){U z=1ULL<<w.a[i][j];if(seen[j]&z)return 6;seen[j]|=z;}
  }
  cout<<endl;
  return 0;
 }
 for(int i=0;i<3;i++)for(int j=i;j<3;j++)cout<<i<<" "<<j<<" "<<tested[i][j]<<endl;
}
