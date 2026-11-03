#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using U=uint64_t;
std::map<U,std::vector<int>> witnesses;
std::vector<unsigned> rows(12);
std::vector<int> ports;
std::vector<U> dp;
U paths(int n){
  const unsigned end=1u<<n; dp.assign(size_t(end)*n,0);
  for(int v=0;v<n;v++)dp[size_t(1u<<v)*n+v]=1;
  for(unsigned s=1;s<end;s++) for(int v=0;v<n;v++) if(s&(1u<<v)){
    U a=dp[size_t(s)*n+v]; if(!a)continue;
    unsigned next=rows[v]&~s;
    while(next){int w=__builtin_ctz(next);next&=next-1;dp[size_t(s|(1u<<w))*n+w]+=a;}
  }
  U total=0;for(int v=0;v<n;v++)total+=dp[size_t(end-1)*n+v];return total;
}
void dfs(int n,int goal,U &cases){
  if(n==goal){U h=paths(n);cases++;if(!witnesses.count(h))witnesses[h]=ports;return;}
  for(int x=0;x<n;x++){
    rows[n]=1u<<x;
    for(int i=0;i<n;i++)if(i!=x)rows[i]|=1u<<n;
    ports.push_back(x);dfs(n+1,goal,cases);ports.pop_back();
    for(int i=0;i<n;i++)rows[i]&=~(1u<<n);
    rows[n]=0;
  }
}
int main(int argc,char**argv){int maxn=argc>1?std::stoi(argv[1]):9;
  if(maxn<1||maxn>11)return 2;
  std::cout<<"{\"layers\":[";
  for(int n=1;n<=maxn;n++){
    U cases=0;dfs(1,n,cases);U gap=1;while(gap==7||gap==21||witnesses.count(gap))gap+=2;
    if(n>1)std::cout<<",";
    std::cout<<"{\"n\":"<<n<<",\"presentations\":"<<cases<<",\"cumulative_counts\":"<<witnesses.size()<<",\"first_missing_allowed\":"<<gap<<",\"max\":"<<witnesses.rbegin()->first<<"}"<<std::flush;
  }
  std::cout<<"],\"witnesses\":{";bool first=true;
  for(auto &[h,w]:witnesses){if(!first)std::cout<<",";first=false;std::cout<<"\""<<h<<"\":[";for(size_t i=0;i<w.size();i++){if(i)std::cout<<",";std::cout<<w[i];}std::cout<<"]";}
  std::cout<<"}}\n";
}
