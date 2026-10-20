// Independent residual census: full affine pair equations are generated separately.
// Binary assignment order, explicit directed 4-cycles, Warshall closure.
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using U=uint32_t;
struct Model { int n,q,low,f; std::vector<U> base; std::vector<std::vector<U>> delta; };
void need(bool b,const std::string&s){if(!b)throw std::runtime_error(s);}
Model read(const char*path){
 std::ifstream in(path);need(bool(in),"model unreadable");Model M;
 need(bool(in>>M.n>>M.q>>M.low>>M.f),"header missing");
 need(M.n>4&&M.n<=28&&M.f>=0&&M.f<31,"invalid model size");
 M.base.resize(M.n);for(U &x:M.base)need(bool(in>>x),"missing base arc");
 M.delta.assign(M.f,std::vector<U>(M.n));for(auto&d:M.delta)for(U&x:d)need(bool(in>>x),"missing toggle");
 std::string extra;need(!(in>>extra),"trailing model content");return M;
}
bool arc(const std::vector<U>&a,int x,int y){return (a[x]>>y)&1u;}
void validate(const Model&M,const std::vector<U>&a){
 const int n=M.n; U full=(U(1)<<n)-1;
 for(int i=0;i<n;i++){
  need((a[i]&~full)==0&&!arc(a,i,i),"bad diagonal or high bit");
  need(arc(a,i,(i+1)%n),"Hamilton cycle absent");
  if(M.low>=0&&i%M.q==M.low)need(std::popcount(a[i])==4,"low orbit degree wrong");
  for(int j=i+1;j<n;j++){
   need(arc(a,i,j)!=arc(a,j,i),"not a tournament");
   need(arc(a,i,j)==arc(a,(i+M.q)%n,(j+M.q)%n),"rotation failed");
   if(M.q==3||M.q==4)need(arc(a,i,j)!=arc(a,n-1-i,n-1-j),"anti-reflection failed");
  }
 }
}
bool fourcycle(const std::vector<U>&a,const std::array<int,4>&s){
 // Test all six Hamilton cycles with smallest vertex first, not a degree proxy.
 for(int b=1;b<4;b++)for(int c=1;c<4;c++)if(c!=b){int d=6-b-c;
  if(arc(a,s[0],s[b])&&arc(a,s[b],s[c])&&arc(a,s[c],s[d])&&arc(a,s[d],s[0]))return true;
 }return false;
}
bool closure_strong(const std::vector<U>&a,U rem){
 // Bit-parallel Warshall; no forward/backward traversal from a chosen vertex.
 std::array<U,28> reach{};int n=a.size();
 for(int i=0;i<n;i++)if((rem>>i)&1)reach[i]=(a[i]&rem)|(U(1)<<i);
 for(int k=0;k<n;k++)if((rem>>k)&1)
  for(int i=0;i<n;i++)if((reach[i]>>k)&1)reach[i]|=reach[k];
 for(int i=0;i<n;i++)if(((rem>>i)&1)&&reach[i]!=rem)return false;
 return rem!=0;
}
bool certificate_valid(const std::vector<U>&a,const std::array<int,5>&w){
 int n=a.size();for(int x:w)if(x<0||x>=n)return false;
 for(int k=0;k<3;k++)if(w[k]>=w[k+1])return false;
 std::array<int,4>s={w[0],w[1],w[2],w[3]};U rem=(U(1)<<n)-1;
 for(int x:s){if(x==w[4])return false;rem&=~(U(1)<<x);}
 int exits=0;for(int x:s)exits+=arc(a,x,w[4]);
 return exits>=3&&fourcycle(a,s)&&closure_strong(a,rem);
}
bool search(const std::vector<U>&a,std::array<int,5>&w,uint64_t&position){
 int n=a.size();U full=(U(1)<<n)-1;position=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++)for(int l=k+1;l<n;l++){
  position++;std::array<int,4>s={i,j,k,l};if(!fourcycle(a,s))continue;
  U rem=full^(U(1)<<i)^(U(1)<<j)^(U(1)<<k)^(U(1)<<l);
  int terminal=-1;
  for(int v=0;v<n;v++)if((rem>>v)&1){int exits=0;for(int x:s)exits+=arc(a,x,v);if(exits>=3){terminal=v;break;}}
  if(terminal>=0&&closure_strong(a,rem)){w={i,j,k,l,terminal};return true;}
 }return false;
}
int main(int argc,char**argv){try{
 need(argc==3||argc==4,"usage: independent_scan model output_certificate [producer_certificate]");
 Model M=read(argv[1]);std::vector<U>a=M.base;
 uint64_t total=uint64_t(1)<<M.f, lastpos=0;std::array<int,5>w{};
 std::ifstream input;if(argc==4){input.open(argv[3],std::ios::binary);need(bool(input),"certificate unreadable");}
 std::ofstream output(argv[2],std::ios::binary);need(bool(output),"output unreadable");
 for(uint64_t index=0;index<total;index++){
  if(index){uint64_t changes=argc==4?uint64_t(1)<<std::countr_zero(index):(index^(index-1));
   for(int b=0;b<M.f;b++)if((changes>>b)&1)for(int v=0;v<M.n;v++)a[v]^=M.delta[b][v];
  }
  validate(M,a); // Every tournament, not a random subset.
  if(argc==4){for(int&x:w){int c=input.get();need(c!=EOF,"early certificate EOF");x=c;}}
  else{uint64_t pos;need(search(a,w,pos),"SURVIVOR at binary assignment "+std::to_string(index));if(pos>lastpos)lastpos=pos;}
  need(certificate_valid(a,w),"invalid witness at assignment "+std::to_string(index));
  for(int x:w)output.put(static_cast<char>(x));
 }
 if(argc==4)need(input.get()==EOF,"extra certificate records");
 need(bool(output),"certificate write failed");
 std::cout<<"{\"n\":"<<M.n<<",\"q\":"<<M.q<<",\"low\":"<<M.low<<",\"free\":"<<M.f<<",\"checked\":"<<total<<",\"survivors\":0,\"mode\":\""<<(argc==4?"producer_certificate_gray_order":"independent_binary_search")<<"\",\"last_subset_position\":"<<lastpos<<"}"<<std::endl;
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
