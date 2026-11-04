#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using U=std::uint64_t;
using Row=std::array<U,9>;

Row ends(int n,const std::array<unsigned,9>& adj){
 std::array<unsigned,9> pred{};
 for(int u=0;u<n;++u) for(int v=0;v<n;++v) if((adj[u]>>v)&1) pred[v]|=1u<<u;
 std::vector<Row> dp(1<<n);
 for(unsigned mask=1;mask<(1u<<n);++mask) for(int v=0;v<n;++v) if(mask&(1u<<v)){
  unsigned rest=mask^(1u<<v),todo=rest&pred[v];
  if(!rest) dp[mask][v]=1;
  while(todo){int u=__builtin_ctz(todo);todo&=todo-1;dp[mask][v]+=dp[rest][u];}
 }
 return dp.back();
}

int main(int argc,char**argv){
 try{
  if(argc!=2) throw std::runtime_error("expected input path");
  std::ifstream input(argv[1]);int orders;input>>orders;
  std::cout<<"{\"orders\":[";
  for(int i=0;i<orders;++i){
   int n,count;input>>n>>count;
   if(n<1||n>9||count<1)throw std::runtime_error("invalid order/count");
   U sumH=0,parity_balanced=0;std::vector<U> balanced,even_parity_codes;
   for(int j=0;j<count;++j){
    U code;input>>code;if(!input)throw std::runtime_error("truncated input");
    std::array<unsigned,9> a{},rev{};int bit=0;
    for(int u=0;u<n;++u)for(int v=u+1;v<n;++v,++bit)
     if((code>>bit)&1)a[u]|=1u<<v;else a[v]|=1u<<u;
    for(int u=0;u<n;++u)rev[u]=((1u<<n)-1)^(1u<<u)^a[u];
    auto e=ends(n,a),s=ends(n,rev);U h=0;bool same=true,same_parity=true;
    for(int u=0;u<n;++u){h+=s[u];same &= s[u]==e[u];same_parity &= (s[u]&1)==(e[u]&1);}
    if(!(h&1))throw std::runtime_error("Hamilton parity failure");
    sumH+=h;if(same)balanced.push_back(code);parity_balanced+=same_parity;
    if(same_parity && n%2==0)even_parity_codes.push_back(code);
   }
   if(i)std::cout<<",";
   std::cout<<"{\"n\":"<<n<<",\"classes\":"<<count<<",\"sum_H\":"<<sumH
            <<",\"balanced_count\":"<<balanced.size()<<",\"parity_balanced_count\":"<<parity_balanced<<",\"balanced_codes\":[";
   for(std::size_t j=0;j<balanced.size();++j){if(j)std::cout<<",";std::cout<<balanced[j];}
   std::cout<<"],\"even_parity_balanced_codes\":[";
   for(std::size_t j=0;j<even_parity_codes.size();++j){if(j)std::cout<<",";std::cout<<even_parity_codes[j];}
   std::cout<<"]}";
   std::cerr<<"finished "<<n<<" balance "<<balanced.size()<<" / "<<count<<"\n";
  }
  std::cout<<"]}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
