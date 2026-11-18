#include <random>
#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
int main(int argc,char**argv){
 if(argc!=3)return 2;
 const int n=stoi(argv[1]),target=stoi(argv[2]);
 if(n<3||n>11)return 2;
 const int f=(1<<n)-1;
 mt19937_64 rng(838289);
 for(int attempt=0;attempt<200000;attempt++){
  int m[11]{},S[2048]{},d[2048][11]{};
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
   bool flip=(rng()%16<(attempt%2?4:8));
   if(flip)m[j]|=1<<i;else m[i]|=1<<j;
  }
  for(int i=0;i<n;i++)d[1<<i][i]=1;
  for(int s=1;s<=f;s++)for(int q=s;q;q&=q-1){
   int i=__builtin_ctz(unsigned(q)),p=d[s][i];S[s]+=p;
   for(int t=m[i]&~s;p&&t;t&=t-1)d[s|(t&-t)][__builtin_ctz(unsigned(t))]+=p;
  }
  if(S[f]!=target)continue;
  long long q=0;for(int s=1;s<f;s++)q+=(long long)S[s]*S[f^s];
  if(q%7){cout<<"attempt "<<attempt<<" H "<<target<<" q2 "<<q<<" rows ";for(int i=0;i<n;i++)cout<<m[i]<<' ';cout<<'\n';return 0;}
 }
 return 1;
}
