#include <array>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
int main(){
 int n,H;
 char sep;
 map<int,vector<int>> hits;
 int count=0;map<int,vector<int>> lows;
 while(cin>>n>>H){
  vector<int>A(n);for(int i=0;i<n;i++){if(i)cin>>sep;cin>>A[i];}
  if(n>9)return 2;
  int f=(1<<n)-1, E[512][9]{},B[512][9]{},I[9]{}, hp[512]{1};
  for(int i=0;i<n;i++){E[1<<i][i]=B[1<<i][i]=1;
   for(int j=0;j<n;j++) if(A[j]>>i&1)I[i]|=1<<j;
  }
  for(int s=1;s<=f;s++)for(int u=0;u<n;u++){
   for(int t=A[u]&~s;t&&E[s][u];t&=t-1)E[s|(t&-t)][__builtin_ctz(unsigned(t))]+=E[s][u];
   for(int t=I[u]&~s;t&&B[s][u];t&=t-1)B[s|(t&-t)][__builtin_ctz(unsigned(t))]+=B[s][u];
   hp[s]+=E[s][u];
  }
  if(hp[f]!=H)return 2;
  int R[9]{},C[9][9]{};
  for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(E[s][i]){
   R[i]+=E[s][i]*hp[f^s];
   for(int j=0;j<n;j++)if(B[f^s][j])C[i][j]+=E[s][i]*B[f^s][j];
  }
  for(int b=1;b<=f;b++){
   int G=0,D=0;
   for(int i=0;i<n;i++)if(b>>i&1){G+=R[i];for(int j=0;j<n;j++)if(b>>j&1)D+=C[i][j];}
   int sign=H-D;
   if(sign==1 || sign==-1) {
    int key=G*sign;
    if(!lows.count(key)){vector<int> a={n,H,b,G};a.insert(a.end(),A.begin(),A.end());lows[key]=a;}
   }
   if(sign!=7&&sign!=-7)continue;
   if(G%7)continue;
   int key=G/7*(sign/7);
   if(!hits.count(key)){vector<int> a={n,H,b,G};a.insert(a.end(),A.begin(),A.end());hits[key]=a;}

  }
  count++;
 }
 for(auto [c,V]:hits)if(lows.count(c)){cout<<c;for(int a:V)cout<<' '<<a;cout<<" lower";
 for(int a:lows[c])cout<<' '<<a;cout<<'\n';}
 cout<<count<<" sources "<<hits.size()<<" different signed coefficients\n";
}
