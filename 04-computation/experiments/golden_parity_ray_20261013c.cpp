// Golden G on positive integral rays; strict stopping at a fixed point.
// Output exact bin sums for integer sources n and n*(1+phi); no scalar closure assumed.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
using I=int64_t;using U=uint64_t;
void ck(bool b){if(!b)throw std::runtime_error("control");}
struct State{I a,b;};
State G(State x){if((x.a&1)==0&&(x.b&1)==0)return{x.a/2,x.b/2};return{(x.b+(x.b&1))/2,(x.a+x.b+((x.a+x.b)&1))/2};}
bool eq(State x,State y){return x.a==y.a&&x.b==y.b;}
int clockG(State x){for(int t=0;t<1000;t++){State y=G(x);if(eq(x,y))return t;x=y;}throw std::runtime_error("clock cap");}
int main(int argc,char**argv){try{ck(argc==3);int E=std::stoi(argv[1]);ck(E>=10&&E<=24);U N=1ULL<<E;
 const int B=1024;std::vector<U>edges;for(int k=8;k<E;k++)for(int j=0;j<B;j++){U n=U(std::floor(std::exp2((long double)k+(long double)j/B)));if(edges.empty()||n>edges.back())edges.push_back(n);}edges.push_back(N+1);
 std::ofstream o(argv[2]);o<<"{\"N\":"<<N<<",\"bins\":[\n";size_t pos=0;I sums[2]={0,0};U count=0,controls=0,maxclock=0;bool comma=false;
 for(U n=1;n<=N;n++){int t[2]={clockG({I(n),0}),clockG({I(n),I(n)})};for(int r=0;r<2;r++){maxclock=std::max(maxclock,U(t[r]));if(n<=65536||((n&(n-1))==0)){int u=clockG({I(2*n),r?I(2*n):0});ck(u==t[r]+1);controls++;}}
  if(n>=edges[0]){sums[0]+=t[0]%2?-1:1;sums[1]+=t[1]%2?-1:1;++count;if(n+1==edges[pos+1]){if(comma)o<<",\n";comma=true;o<<'['<<edges[pos]<<','<<edges[pos+1]<<','<<count<<','<<sums[0]<<','<<sums[1]<<']';++pos;sums[0]=sums[1]=0;count=0;}}
 }o<<"],\"doubling_controls\":"<<controls<<",\"max_clock\":"<<maxclock<<"}\n";ck(bool(o));}catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}}
