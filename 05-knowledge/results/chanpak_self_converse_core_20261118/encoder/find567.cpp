#include <array>
#include <iostream>
#include <vector>
#include <utility>
using namespace std;
int main(){
  constexpr int n=9, N=1<<n;
  vector<pair<int,int>> o;
  for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)
    if (pair{i,j}<=pair{n-1-j,n-1-i})o.push_back({i,j});
  for(unsigned w=0;w<(1u<<o.size());w++){
    array<int,n>A{};
    for(int i=0;i<n-1;i++)A[i]=1<<(i+1);
    for(unsigned b=0;b<o.size();b++){
      auto [i,j]=o[b];auto [k,l]=pair{n-1-j,n-1-i};
      if(w>>b&1){A[i]|=1<<j;A[k]|=1<<l;}else{A[j]|=1<<i;A[l]|=1<<k;}
    }
    array<array<int,n>,N>D{};
    for(int i=0;i<n;i++)D[1<<i][i]=1;
    for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(D[s][i]){
      int b=A[i]&~s;
      while(b){int v=b&-b;b-=v;D[s|v][__builtin_ctz(v)]+=D[s][i];}
    }
    int H=0;for(int z:D[N-1])H+=z;
    bool strong=true;
    for(int i=1;i<n;i++){
      bool has=false;for(int j=i;j<n;j++)if(A[j]&((1<<i)-1))has=true;
      if(!has)strong=false;
    }
    if(H==567&&strong){cout<<w;for(int z:A)cout<<" "<<z;cout<<"\n";return 0;}
  }
  cout<<"No witness at order 9\n";
}
