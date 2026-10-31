// Backward endpoint DP; independent of constructor profiles and forward DP.
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include <limits>
#include <stdexcept>
using namespace std;
int main(int argc,char**argv){
 if(argc<2){cerr<<"Usage: adjacency_dp ADJACENCY_FILE...\n";return 2;}
 for(int z=1;z<argc;++z){
  ifstream in(argv[z]);int n=0;uint64_t target=0;
  if(!(in>>n>>target)||n<1||n>22){cerr<<"Invalid header/order\n";return 2;}
  vector<string> rows(n);vector<uint32_t> incoming(n);
  for(int i=0;i<n;++i){
   if(!(in>>rows[i])||rows[i].size()!=size_t(n)||rows[i].find_first_not_of("01")!=string::npos){cerr<<"Invalid adjacency row\n";return 2;}
  }
  string extra;if(in>>extra){cerr<<"Trailing adjacency data\n";return 2;}
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   if((i==j&&rows[i][j]!='0')||(i!=j&&rows[i][j]==rows[j][i])){cerr<<"Not a tournament\n";return 2;}
   if(rows[i][j]=='1')incoming[j]|=1u<<i;
  }
  // Four repair witnesses have n<=20, hence every partial count<=20!<2^64.
  // Optional reconstructed samples may have n=21 or 22: check every addition.
  auto add=[](uint64_t a,uint64_t b){if(b>numeric_limits<uint64_t>::max()-a)throw overflow_error("DP count overflow");return a+b;};
  vector<uint32_t> dp(size_t(n)*(1u<<n));
  for(uint32_t mask=1;mask<(1u<<n);++mask){uint32_t ends=mask;
   while(ends){int v=__builtin_ctz(ends);ends&=ends-1;uint32_t rest=mask^(1u<<v);uint64_t count=0;
    if(!rest)count=1;else{uint32_t pred=rest&incoming[v];while(pred){int u=__builtin_ctz(pred);pred&=pred-1;count=add(count,dp[size_t(rest)*n+u]);}}
    if(count>numeric_limits<uint32_t>::max())throw overflow_error("32-bit storage overflow");dp[size_t(mask)*n+v]=uint32_t(count);
   }
  }
  uint64_t total=0;for(int v=0;v<n;++v)total=add(total,dp[size_t((1u<<n)-1)*n+v]);
  cout<<"order "<<n<<" target "<<target<<" actual "<<total<<'\n';if(total!=target)return 1;
 }
}
