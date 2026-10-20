// Independent all-cross-graph evaluation using the ordinary monomial basis.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
int main(int argc,char**argv){if(argc!=4)return 2;int n=stoi(argv[1]);if(n<1||n>7)return 2;int e=n*(n-1)/2,B=1<<e,loc[7][7];int q=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)loc[i][j]=loc[j][i]=q++;
 vector<int> reps;ifstream ri(argv[2]);int a;while(ri>>a)reps.push_back(a);ifstream raw(argv[3],ios::binary);if(!raw)return 3;uint64_t checks=0,updates=0,prefixes=0;auto start=chrono::steady_clock::now();
 for(int a:reps){vector<int32_t> poly(B,0),reference(B);raw.read((char*)reference.data(),B*sizeof(int32_t));if(raw.gcount()!=B*sizeof(int32_t))return 4;
  vector<int> order(n);iota(order.begin(),order.end(),0);
  do{for(int chosen=0;chosen<(1<<n);++chosen){if(chosen>>order.back()&1)continue;bool valid=true;int positive=0,negative=0;
   for(int t=1;t<n;++t){int u=order[t-1],v=order[t],su=chosen>>u&1,sv=chosen>>v&1,bit=loc[u][v];
    if(su==sv){bool uv=u<v?(a>>bit&1):!(a>>bit&1);if(uv==bool(su)){valid=false;break;}}
    else if(su)negative|=1<<bit;else positive|=1<<bit;
   }
   if(!valid)continue;prefixes++;int subset=negative;
   while(true){poly[positive|subset]+=__builtin_parity((unsigned)subset)?-1:1;updates++;if(!subset)break;subset=(subset-1)&negative;}
  }}while(next_permutation(order.begin(),order.end()));
  for(int i=0;i<e;++i){int bit=1<<i;for(int mask=0;mask<B;++mask)if(mask&bit)poly[mask]+=poly[mask^bit];}
  for(int b=0;b<B;++b){if(poly[b]!=reference[b]){cerr<<"Mismatch A="<<a<<" b="<<b<<" ordinary="<<poly[b]<<" Walsh="<<reference[b]<<endl;return 5;}checks++;}
 }
 char c;if(raw.get(c))return 6;
 cout<<"{\"k\":"<<n<<",\"method\":\"ordinary monomial subset-zeta replay\",\"tournament_types\":"<<reps.size()<<",\"pointwise_comparisons\":"<<checks<<",\"all_equal\":true,\"compatible_prefixes\":"<<prefixes<<",\"coefficient_updates\":"<<updates<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
