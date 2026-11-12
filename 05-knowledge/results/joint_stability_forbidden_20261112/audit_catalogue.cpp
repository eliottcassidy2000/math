// Independent C++ catalogue audit. All accumulators are exact 64-bit integers.
// Input: ordinary tournament upper triangles in row order, one per line.
#include <array>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstdint>
using I=std::int64_t;
int main(int argc,char**argv){
 if(argc!=2)throw std::runtime_error("order required");int n=std::stoi(argv[1]);
 if(n<2||n>10)throw std::runtime_error("proved fixed order bound");
 int N=1<<n,full=N-1;I classes=0,strongs=0,units=0,masksum=0,strongunits=0,strongmasks=0;
 std::map<I,I> hist;std::string line;std::vector<I>E(N*n),S(N*n),W(n*n),weight(N);std::vector<int>r(n),inc(n);
 while(std::getline(std::cin,line)){
  if(!line.empty()&&line.back()=='\r')line.pop_back();
  if(int(line.size())!=n*(n-1)/2)throw std::runtime_error("triangle length");
  std::fill(r.begin(),r.end(),0);std::fill(inc.begin(),inc.end(),0);int pos=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   char c=line[pos++];if(c!='0'&&c!='1')throw std::runtime_error("triangle alphabet");
   int a=c=='1'?i:j,b=c=='1'?j:i;r[a]|=1<<b;inc[b]|=1<<a;
  }
  bool strong=true;
  for(int reverse=0;reverse<2;reverse++){
   int reached=1,todo=1;while(todo){int v=__builtin_ctz((unsigned)todo);todo&=todo-1;int add=(reverse?inc[v]:r[v])&~reached;reached|=add;todo|=add;}
   if(reached!=full)strong=false;
  }
  std::fill(E.begin(),E.end(),0);std::fill(S.begin(),S.end(),0);
  for(int v=0;v<n;v++)E[(1<<v)*n+v]=S[(1<<v)*n+v]=1;
  for(int m=1;m<N;m++)if(m&(m-1))for(int vs=m;vs;vs&=vs-1){
   int v=__builtin_ctz((unsigned)vs),rest=m^(1<<v);
   for(int us=rest&inc[v];us;us&=us-1)E[m*n+v]+=E[rest*n+__builtin_ctz((unsigned)us)];
   for(int us=rest&r[v];us;us&=us-1)S[m*n+v]+=S[rest*n+__builtin_ctz((unsigned)us)];
  }
  I H=0;for(int v=0;v<n;v++)H+=E[full*n+v];if(H%2!=1)throw std::runtime_error("parity");
  std::fill(W.begin(),W.end(),0);
  for(int m=1;m<full;m++)for(int is=m;is;is&=is-1){int i=__builtin_ctz((unsigned)is);I a=E[m*n+i];if(!a)continue;
   for(int js=full^m;js;js&=js-1){int j=__builtin_ctz((unsigned)js);I v=a*S[(full^m)*n+j];W[i*n+j]+=v;W[j*n+i]+=v;}
  }
  I count=0;weight[0]=0;
  for(int m=1;m<N;m++){int v=__builtin_ctz((unsigned)m),rest=m^(1<<v);I w=weight[rest];for(int us=rest;us;us&=us-1)w+=W[v*n+__builtin_ctz((unsigned)us)];weight[m]=w;count+=w==H-1;}
  classes++;strongs+=strong;hist[H]++;units+=count>0;masksum+=count;
  if(strong){strongunits+=count>0;strongmasks+=count;}
  if(classes%1000000==0)std::cerr<<classes<<" ordinary classes audited\n";
 }
 std::cout<<"{\"classes\":"<<classes<<",\"strong_classes\":"<<strongs<<",\"unit_sources\":"<<units<<",\"unit_masks\":"<<masksum<<",\"strong_unit_sources\":"<<strongunits<<",\"strong_unit_masks\":"<<strongmasks<<",\"H_histogram\":{";
 bool first=true;for(auto [h,c]:hist){if(!first)std::cout<<',';first=false;std::cout<<'"'<<h<<"\":"<<c;}std::cout<<"}}\n";
}
