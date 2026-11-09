#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdint>
using I=__int128_t; using L=long long;
struct Y {L a,b,c;};
L g(L v){L z=1;while(v){if(!(v&1))z+=(v>>1)+((v>>1)&1);v>>=1;}return z;}
L ceildiv(L a,L b){return a>=0?(a+b-1)/b:a/b;}
int main(int argc,char**argv){std::ifstream in(argv[1]);std::vector<Y> ys;Y y;while(in>>y.a>>y.b>>y.c)ys.push_back(y);std::sort(ys.begin(),ys.end(),[](Y x,Y y){return x.c<y.c || (x.c==y.c&&x.a<y.a);});for(int z=2;z<argc;z++){L N=std::stoll(argv[z]); if(N<1||N>1000000000000000LL)return 2; L branches=0;bool done=false;for(auto y:ys){L M=2;while(I(y.c)*y.c*M*M<4*I(y.a)*N||I(y.c)*M<4*y.a)M*=2;L rev=0;for(L r=0;r<M;r++){branches++;L t=(r>=M/2),E=1-t,C=M-1-rev,Lc=y.c*M+y.a*C,Q=y.a*(g(r)-1)-y.b*(1-(r&1))+y.c*(1+r),x=N-Q;if(x>=y.a){L lo=std::max(0LL,ceildiv(x-2*y.a,Lc+y.a)),hi=(x-y.a)/Lc;if(hi-lo>0){std::cerr<<"width fail\n";return 3;}if(lo<=hi){L v=M*lo+r;if(I(y.a)*g(v)-y.b*(1-(v&1))+I(y.c)*(1+v)==N){std::cout<<N<<" "<<y.a<<" "<<y.b<<" "<<y.c<<" "<<v<<" "<<branches<<"\n";done=true;break;}}}L bit=M/2;while(bit && (rev&bit)){rev^=bit;bit/=2;}rev^=bit;}if(done)break;}if(!done)std::cout<<N<<" NONE "<<branches<<"\n";}}
