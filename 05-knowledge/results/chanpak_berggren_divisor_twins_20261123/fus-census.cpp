#include <cstdio>
#include <vector>
#include <set>
#include <utility>
#include <cstdlib>
int main(int argc,char**argv){
  const int n=argc>1?atoi(argv[1]):7, N=(1<<n)-1;
  std::vector<std::pair<int,int>> E;
  for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)if(b-a!=1&&b-a!=n-1)E.push_back({a,b});
  std::set<std::pair<int,int>> S;
  int hits=0;
  for(int word=0;word<(1<<E.size());++word){
    int in[8]={};
    for(int a=0;a<n;++a)in[(a+1)%n]|=1<<a;
    for(int i=0;i<(int)E.size();++i){auto [a,b]=E[i];if((word>>i)&1)std::swap(a,b);in[b]|=1<<a;}
    int H[256]={1},dp[256][8]={};
    for(int A=1;A<=N;++A){
      for(unsigned m=A;m;m&=m-1){
        int i=__builtin_ctz(m),p=A^(1<<i);
        if(!p)dp[A][i]=1;
        for(unsigned r=p&in[i];r;r&=r-1)dp[A][i]+=dp[p][__builtin_ctz(r)];
        H[A]+=dp[A][i];
      }
    }
    if(H[N]%7)continue;
    int Q=0;for(int A=0;A<=N;++A)Q+=H[A]*H[N^A];
    if((n==7 && (Q-H[N]==2401 || (2*Q-H[N]==3969 && (H[N]==119||H[N]==175)))) || (n==8&&Q-H[N]==3969)){
      if(S.insert({H[N],Q}).second){for(int i=0;i<n;++i){int out=0;for(int j=0;j<n;++j)if((in[j]>>i)&1)out|=1<<j;std::printf("%d ",out);}std::puts("");}++hits;
    }
  }
  std::printf("%d: %zu masks; %d hits and %zu profiles\n",n,size_t(1)<<E.size(),hits,S.size());
  for(auto [h,q]:S)std::printf("%d %d\n",h,q);
}
