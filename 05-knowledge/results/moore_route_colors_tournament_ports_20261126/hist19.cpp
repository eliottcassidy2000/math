#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cstdint>
using namespace std;
int main(){
 string w;const int N=19,U=(1<<N)-1;
 while(cin>>w){
  if(w.size()!=171||w.find_first_not_of("01")!=string::npos)return 1;
  vector<int>a(N);int j0=0;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){if(w[j0++]=='1')a[i]|=1<<j;else a[j]|=1<<i;}
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){
   int v=(a[i]>>j)&1;
   for(int k=0;k<N;k++)v+=((a[i]>>k)&1)*((a[k]>>j)&1);
   if(v!=(i==j?0:5)){cerr<<"not DRT";return 1;}
  }
  vector<uint64_t>D((U+1)*N);
  for(int i=0;i<N;i++)D[(1<<i)*N+i]=1;
  for(int s=1;s<U;s++)for(int z=s;z;z&=z-1){
   int j=__builtin_ctz(z);const auto d=D[s*N+j];
   for(int e=a[j]&~s;e;e&=e-1){
    int i=__builtin_ctz(e);D[(s|(1<<i))*N+i]+=d;
   }
  }
  uint64_t H=0;map<uint64_t,int> ends,del;
  for(int j=0;j<N;j++){
   H+=D[U*N+j];ends[D[U*N+j]]++;
   uint64_t h=0;
   for(int i=0;i<N;i++)h+=D[(U^(1<<j))*N+i];
   del[h]++;
  }
  cout<<"H "<<H<<" ends";
  for(auto [a,b]:ends)cout<<" "<<a<<":"<<b;
  cout<<" delete";
  for(auto [a,b]:del)cout<<" "<<a<<":"<<b;
  cout<<"\n";
 }
}
