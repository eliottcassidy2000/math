// Exact endpoint subset DP on literal adjacency, with checked uint64 arithmetic.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("usage: subset_dp ADJACENCY");
 std::ifstream f(argv[1]); unsigned n; uint64_t expected;
 if(!(f>>n>>expected)||n<1||n>20)throw std::runtime_error("require order1..20");
 std::vector<uint32_t> out(n),in(n);std::string s;
 for(unsigned i=0;i<n;i++){
  if(!(f>>s)||s.size()!=n)throw std::runtime_error("row size");
  for(unsigned j=0;j<n;j++){
   if(s[j]!='0'&&s[j]!='1')throw std::runtime_error("nonbinary");
   if(s[j]=='1'){out[i]|=1u<<j;in[j]|=1u<<i;}
  }
 }
 if(f>>s)throw std::runtime_error("trailing input");
 for(unsigned i=0;i<n;i++)for(unsigned j=0;j<=i;j++)if(i==j ? (out[i]>>i&1) : ((out[i]>>j&1)+(out[j]>>i&1)!=1))throw std::runtime_error("not tournament");
 const unsigned full=(1u<<n)-1;std::vector<uint64_t> dp((uint64_t)(full+1)*n);
 auto add=[](uint64_t&a,uint64_t b){if(std::numeric_limits<uint64_t>::max()-a<b)throw std::runtime_error("overflow");a+=b;};
 for(unsigned S=1;S<=full;S++){
  unsigned vertices=S;while(vertices){unsigned v=__builtin_ctz(vertices);vertices&=vertices-1;unsigned R=S^(1u<<v);auto &value=dp[(uint64_t)S*n+v];if(!R){value=1;continue;}
   unsigned pred=R&in[v];while(pred){unsigned u=__builtin_ctz(pred);pred&=pred-1;add(value,dp[(uint64_t)R*n+u]);}
  }
 }
 uint64_t total=0;for(unsigned v=0;v<n;v++)add(total,dp[(uint64_t)full*n+v]);
 if(total!=expected)throw std::runtime_error("count mismatch");
 std::cout<<"{\"status\":\"passed\",\"order\":"<<n<<",\"exact_H\":"<<total<<",\"algorithm\":\"full subset endpoint DP\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
