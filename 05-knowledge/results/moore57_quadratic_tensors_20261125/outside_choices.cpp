// Exact lower bound for imports from outside the union of nearest projections.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
using M=array<U,3>;
constexpr U ALL=(1ULL<<52)-1;
U bit(int x){return 1ULL<<x;}
struct E {array<int,3> r;int y;};
vector<E> pin[52], ext[52];
uint64_t nodes,solutions[5];
bool free(M used,E e){
 for(int j=0;j<3;j++)if(used[j]&bit(e.r[j]))return false;
 return true;
}
M add(M used,E e){for(int j=0;j<3;j++)used[j]|=bit(e.r[j]);return used;}
bool import_test(U open,M used,int k){
 if(k<=0)return true;
 if(__builtin_popcountll(open)<k)return false;
 int i=__builtin_ctzll(open);open^=bit(i);
 for(E e:ext[i])if(free(used,e)&&import_test(open,add(used,e),k-1))return true;
 return import_test(open,used,k);
}
void dfs(U left,U holes,M used,int budget){
 ++nodes;
 if(!left){
  int k=__builtin_popcountll(holes);
  ++solutions[k];
  if(import_test(holes,used,k-1)){
   cerr<<"unexpected 51-branch patch with at most three imports\n";exit(5);
  }
  return;
 }
 int h=-1,n=1000,obligate=0;
 for(U R=left;R;R&=R-1){
  int i=__builtin_ctzll(R),m=0;
  for(E e:pin[i])if(free(used,e))++m;
  if(!m&&++obligate>budget)return;
  if(m<n){h=i;n=m;}
 }
 left^=bit(h);
 if(budget)dfs(left,holes|bit(h),used,budget-1);
 for(E e:pin[h])if(free(used,e))dfs(left,holes,add(used,e),budget);
}
int main(int argc,char**argv){
 if(argc!=2)return 2;
 ifstream f(argv[1]);
 int g,k,y,n,count=0;
 U P[52]{};
 while(f>>g>>k){
  ++count;
  for(U&t:P){f>>y;if(y<0||y>=56)return 2;t|=bit(y);}
 }
 if(count!=25)return 2;
 cin>>n;
 E m;cin>>k>>count>>m.y>>m.r[0]>>m.r[1]>>m.r[2];
 if(n!=53||k!=7||count!=1||m.y!=7)return 2;
 M used=add({},m);
 for(int i=0;i<52;i++){
  cin>>k>>count;
  while(count--){
   E e;cin>>e.y>>e.r[0]>>e.r[1]>>e.r[2];
   if(!free(used,e))continue;
   (P[i]&bit(e.y)?pin[i]:ext[i]).push_back(e);
  }
 }
 dfs(ALL,0,used,4);
 for(int i=0;i<=4;i++)cout<<(52-i)<<" "<<solutions[i]<<endl;
 cout<<"nodes "<<nodes<<endl;
 if(solutions[0]!=0||solutions[1]!=0||solutions[2]!=3||solutions[3]!=173||solutions[4]!=5172)return 6;
}
