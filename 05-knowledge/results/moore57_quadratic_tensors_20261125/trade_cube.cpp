// Enumerate all smallest trades in a linear four-partite contact board.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
using E=array<int,4>;
using Trade=array<E,2>;
int main() {
 int n;cin>>n;
 if(n!=53)return 2;
 int lookup[4][4][56][56];
 for(auto &a:lookup)for(auto&b:a)for(auto&c:b)for(int &x:c)x=-1;
 vector<E>C;
 vector<int>Y,ref(n);
 E marker;
 for(int i=0;i<n;i++){
  int k;cin>>ref[i]>>k;
  while(k--){
   int y;E t{i,0,0,0};cin>>y>>t[1]>>t[2]>>t[3];
   if(!i){marker=t;continue;}
   bool ok=true;
   for(int j=1;j<4;j++)if(t[j]==marker[j])ok=false;
   if(ok){C.push_back(t);Y.push_back(y);}
  }
 }
 if(C.size()!=1741)return 3;
 for(int a=0;a<(int)C.size();a++)for(int i=0;i<4;i++)for(int j=i+1;j<4;j++){
  int x=C[a][i],y=C[a][j];
  if(lookup[i][j][x][y]!=-1)return 3;
  lookup[i][j][x][y]=lookup[j][i][y][x]=a;
 }
 auto conflict=[&](int a,int b){for(int i=0;i<4;i++)if(C[a][i]==C[b][i])return true;return false;};
 vector<E> perm;
 E P{0,1,2,3};
 do{bool ok=true;for(int i=0;i<4;i++)if(P[i]==i)ok=false;if(ok)perm.push_back(P);}while(next_permutation(P.begin(),P.end()));
 set<Trade>T;
 int raw=0;
 for(int a=0;a<(int)C.size();a++)for(int b=a+1;b<(int)C.size();b++)if(!conflict(a,b))for(auto p:perm){
  E R;
  bool ok=true;
  for(int i=0;i<4;i++){
   int j=p[i];R[i]=lookup[i][j][C[a][i]][C[b][j]];
   if(R[i]<0){ok=false;break;}
   for(int j=0;j<i;j++)if(conflict(R[i],R[j])){ok=false;break;}
   if(!ok)break;
  }
  if(!ok)continue;
  int holes[4][2];
  for(int j=0;j<4;j++){
   int h=0;for(int r:R){int u=C[r][j];if(u!=C[a][j]&&u!=C[b][j])holes[j][h++]=u;}
   if(h!=2)return 4;
  }
  for(int m=0;m<8;m++){
   E X,W;
   for(int j=0;j<4;j++){X[j]=holes[j][(m>>j)&1];W[j]=holes[j][1-((m>>j)&1)];}
   int x=lookup[0][1][X[0]][X[1]],w=lookup[0][1][W[0]][W[1]];
   if(x<0||w<0||C[x]!=X||C[w]!=W)continue;
   if(x==a||x==b||w==a||w==b)return 4;
   for(int r:R)if(r==x||r==w||r==a||r==b)return 4;
   E L{a,b,x,w},U=R;
   sort(L.begin(),L.end());sort(U.begin(),U.end());
   if(U<L)swap(L,U);
   ++raw;T.insert({L,U});
  }
 }
 if(raw!=996||T.size()!=83)return 5;
 cout<<"# sides: (board row, sheet) pairs; first is lexically smaller\n";
 for(const auto&t:T){
  for(const auto &s:t)for(int j:s)cout<<C[j][0]<<" "<<Y[j]<<" ";
  cout<<"\n";
 }
}
