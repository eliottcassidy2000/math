#define main old_census_main
#include "../../two_port_beyond_nine_20261103/census/census.cpp"
#undef main
unsigned parity_ends(const Rows&a,const std::array<unsigned char,1024>&p){
 unsigned dp[1024]{};Rows pred{};for(int v=0;v<10;++v)pred[v]=1023u^(1u<<v)^a[v];
 for(unsigned mask=1;mask<1024;++mask){unsigned vs=mask,value=0;while(vs){int v=__builtin_ctz(vs);vs&=vs-1;unsigned rest=mask^(1u<<v);if(!rest||p[dp[rest]&pred[v]])value|=1u<<v;}dp[mask]=value;}return dp[1023];
}
int main(int argc,char**argv){try{
 need(argc==3||argc==4,"parents flags [limit]");U limit=argc==4?std::stoull(argv[3]):9733056;auto parents=read64(argv[1]);need(parents.size()==191536,"parent count");std::ifstream f(argv[2],std::ios::binary);need(bool(f),"flag file");
 std::array<unsigned char,1024> parity{};for(int i=1;i<1024;++i)parity[i]=parity[i>>1]^(i&1);
 U count=0,address=0,balanced=0;std::vector<U> residual;
 for(;count<limit;++count){int c=f.get();need(c!=EOF,"short input");U delta=0;int shift=0;for(;;){need(shift<63,"varint overflow");delta|=U(c&127)<<shift;if(!(c&128))break;shift+=7;c=f.get();need(c!=EOF,"truncated input");}need(!count||delta,"increasing addresses");address+=delta;need((address>>9)<parents.size(),"address bounds");
  auto a=decode(parents[address>>9],address&511);Rows rev{};for(int i=0;i<10;++i)rev[i]=1023u^(1u<<i)^a[i];unsigned e=parity_ends(a,parity),s=parity_ends(rev,parity);
  if(e==s){++balanced;auto d=make(a,true);bool okay=true;for(int i=0;i<10;++i){need((d.s[i]&1)==((s>>i)&1)&&(d.e[i]&1)==((e>>i)&1),"independent integer endpoint check");for(int j=0;j<10;++j)okay &= ((d.c[i][j]-d.c[j][i])&3)==0;}if(okay)residual.push_back(address);}
  if((count+1)%1000000==0)std::cerr<<"tested="<<count+1<<" balanced="<<balanced<<" residual="<<residual.size()<<"\n";
 }
 if(argc==3)need(f.peek()==EOF,"full bank EOF");std::cout<<"{\"n\":10,\"classes\":"<<count<<",\"complete_bank\":"<<(argc==3?"true":"false")<<",\"endpoint_parity_balanced\":"<<balanced<<",\"combined_residual_addresses\":[";for(size_t i=0;i<residual.size();++i){if(i)std::cout<<",";std::cout<<residual[i];}std::cout<<"]}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
