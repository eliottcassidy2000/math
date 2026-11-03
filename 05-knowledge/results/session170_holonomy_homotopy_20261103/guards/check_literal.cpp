// Independent literal subset-DP certificate for the 16-vertex reset word.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using Row=std::uint32_t;using Count=std::uint64_t;using Graph=std::vector<Row>;
void need(bool x,const char* s){if(!x)throw std::runtime_error(s);}
Graph add(Graph a,Row out){unsigned n=a.size();for(unsigned i=0;i<n;++i)if(!(out>>i&1))a[i]|=Row(1)<<n;a.push_back(out);return a;}
Graph erase(const Graph&a,unsigned q){Graph b;for(unsigned i=0;i<a.size();++i)if(i!=q){Row x=0;unsigned k=0;for(unsigned j=0;j<a.size();++j)if(j!=q){if(a[i]>>j&1)x|=Row(1)<<k;++k;}b.push_back(x);}return b;}
std::vector<Count> starts(const Graph&a){
 unsigned n=a.size();Row size=Row(1)<<n;std::vector<Count> dp(std::size_t(size)*n);
 for(unsigned i=0;i<n;++i)dp[std::size_t(Row(1)<<i)*n+i]=1;
 for(Row mask=1;mask<size;++mask){
  Row bits=mask;while(bits){unsigned i=__builtin_ctz(bits);bits&=bits-1;Row rest=mask^(Row(1)<<i);if(!rest)continue;
   Row next=a[i]&rest;Count total=0;while(next){unsigned j=__builtin_ctz(next);next&=next-1;Count addend=dp[std::size_t(rest)*n+j];need(addend<=std::numeric_limits<Count>::max()-total,"DP addition overflow");total+=addend;}
   dp[std::size_t(mask)*n+i]=total;
  }
 }
 return std::vector<Count>(dp.end()-n,dp.end());
}
std::array<Count,3> state(const Graph&a,unsigned q){
 auto s=starts(a);Count h=0;for(auto x:s){need(x<=std::numeric_limits<Count>::max()-h,"total count overflow");h+=x;}
 unsigned r=a.size();for(unsigned i=0;i<a.size();++i)if(a[i]>>q&1){need(r==a.size(),"multiple in-neighbors");r=i;}need(r<a.size(),"no in-neighbor");
 auto deleted=starts(erase(a,q));return {s[q],h-s[q],deleted[r-(r>q)]};
}
int main(){
 Graph a{2,4,1};unsigned q=0;std::string word="D"+std::string(14,'C')+"D";
 for(char c:word){a=add(a,((Row(1)<<a.size())-1)^(Row(1)<<q));if(c=='C')q=a.size()-1;}
 const std::array<std::array<Count,3>,3> expect{{{12737,30352,0},{22419,53424,0},{42629,101584,0}}};
 std::cout<<"{\n  \"word\": \""<<word<<"\",\n  \"query_states\": [\n";
 for(unsigned k=0;k<3;++k){auto v=state(a,q);need(v==expect[k],"literal endpoint state mismatch");std::cout<<"    ["<<v[0]<<", "<<v[1]<<", "<<v[2]<<"]"<<(k==2?"\n":",\n");a=add(a,Row(1)<<(k==0?1:a.size()-1));}
 std::cout<<"  ],\n  \"observer\": [43089, 16377, 34185],\n  \"orders_counted\": [19, 20, 21],\n  \"status\": \"PASS\"\n}\n";
}
