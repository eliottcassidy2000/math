#include <bits/stdc++.h>
using namespace std;
int hp(int n,vector<pair<int,int>> &E,int w){
 int R[8]={},D[128][7]={};
 for(auto [i,j]:E){if(w&1)R[i]|=1<<j;else R[j]|=1<<i;w>>=1;}
 for(int i=0;i<n;i++)D[1<<i][i]=1;
 for(int s=1;s<(1<<n);s++)for(int i=0;i<n;i++)if(D[s][i]){
  int M=R[i]&~s;
  while(M){int b=M&-M;M^=b; D[s|b][__builtin_ctz(b)]+=D[s][i];}
 }
 return accumulate(D[(1<<n)-1],D[(1<<n)-1]+n,0);
}
bool prime(int x){if(x<2)return false;for(int i=2;i*i<=x;i++)if(x%i==0)return false;return true;}
int main(){
 for(int n=3;n<=7;n++){
  vector<pair<int,int>> E;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)E.push_back({i,j});
  int S=1<<E.size();
  vector<int> h(S);int ng=0,twins=0;
  for(int i=0;i<S;i++)h[i]=hp(n,E,i);
  for(int i=0;i<S;i++)if(h[i]%7==0&&h[i]%49!=0)
   for(unsigned b=0;b<E.size();b++){
    int j=i^(1<<b);if(h[j]==h[i]+14&&h[j]%49!=0) {
     ng++;
     if(prime(h[i]/7)&&prime(h[j]/7)){
      if(twins++==0)cout<<n<<' '<<i<<' '<<b<<' '<<h[i]<<' '<<h[j]<<endl;
     }
    }
   }
  cout<<"totals "<<n<<' '<<ng<<' '<<twins<<endl;
 }
}
