#include <iostream>
#include <vector>
#include <set>
#include <cstdint>
#include <cstdlib>
#include <map>
using namespace std;
int main(int argc,char**argv){int maxn=argc>1?atoi(argv[1]):9;set<int>all={1};cout<<"[";bool first=true;
 for(int n=2;n<=maxn;n++){
  long long graphs=0;set<int>vals;map<int,pair<int,uint64_t>> witnesses;
  for(int p=1;p<=n/2;p++){
   int q=n-p,E=p*q;uint64_t codes=1ull<<E;vector<unsigned>neigh(n),dom(1u<<n);
   for(uint64_t code=0;code<codes;code++){
    for(int i=0;i<n;i++)neigh[i]=1u<<i;
    for(int i=0;i<p;i++)for(int j=0;j<q;j++)if(code>>(i*q+j)&1){neigh[i]|=1u<<(p+j);neigh[p+j]|=1u<<i;}
    int d=0;dom[0]=0;unsigned full=(1u<<n)-1;
    for(unsigned S=1;S<=full;S++){unsigned b=S&-S;dom[S]=dom[S^b]|neigh[__builtin_ctz(b)];d+=dom[S]==full;}
    vals.insert(d);witnesses.emplace(d,make_pair(p,code));graphs++;
   }
  }
  all.insert(vals.begin(),vals.end());if(!first)cout<<",";first=false;
  cout<<"{\"n\":"<<n<<",\"labelled_bipartition_presentations\":"<<graphs<<",\"values\":[";bool f=true;for(int x:vals){if(!f)cout<<",";f=false;cout<<x;}cout<<"],\"cumulative_missing_below100\":[";f=true;for(int x=1;x<100;x+=2)if(!all.count(x)){if(!f)cout<<",";f=false;cout<<x;}cout<<"],\"witnesses\":[";f=true;for(auto z:witnesses){if(!f)cout<<",";f=false;cout<<"["<<z.first<<","<<z.second.first<<","<<z.second.second<<"]";}cout<<"]}"<<flush;
 }
 cout<<"]\n";
}
