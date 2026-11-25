#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
using namespace std;
int count(const array<int,11>& A) {
 int D[2048][11];const int n=11,full=(1<<n)-1;
 for(int s=1;s<=full;s++)for(int b=s;b;b&=b-1){
  int i=__builtin_ctz((unsigned)b),r=s^(1<<i),v=r==0;
  for(int g=r&~A[i];g;g&=g-1)v+=D[r][__builtin_ctz((unsigned)g)];
  D[s][i]=v;
 }
 int H=0;for(int i=0;i<n;i++)H+=D[full][i];return H;
}
int main(){
 vector<pair<int,int>> pairs;array<int,11>A{};
 auto R=[](int d){return d==1||d==3||d==4||d==5||d==9;};
 for(int i=0;i<11;i++)for(int j=i+1;j<11;j++){
  pairs.emplace_back(i,j);
  if(R(j-i))A[i]|=1<<j;else A[j]|=1<<i;
 }
 auto flip=[&](int i){auto [u,v]=pairs[i];A[u]^=1<<v;A[v]^=1<<u;};
 set<int> heights[6],fringe;
 heights[0].insert(count(A));
 auto run=[&](auto& self,int pos,int len,int shell)->void{
  if(!len){heights[shell].insert(count(A));return;}
  for(int i=pos;i+len<=55;i++){flip(i);self(self,i+1,len-1,shell);flip(i);}
 };
 flip(0);for(int k=1;k<=5;k++)run(run,1,k-1,k);flip(0);
 for(int k=0;k<=5;k++){
  cout<<"shell "<<k<<" distinct "<<heights[k].size()<<'\n';
  for(int x:heights[k])if(x>80555)fringe.insert(x);
 }
 cout<<"verified heights above 80555: "<<fringe.size()<<'\n';
 assert(*heights[0].begin()==95095&&fringe.size()==1308);
}
