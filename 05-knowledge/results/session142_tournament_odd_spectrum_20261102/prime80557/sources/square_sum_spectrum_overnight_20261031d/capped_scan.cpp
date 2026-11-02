// Bounded source bank, streamed singleton axis and a packed odd support map.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=int64_t;
int main(int argc,char**argv){try{
 if(argc<5)throw runtime_error("usage: scan limit cap bitmap profile-files...");
 I N=stoll(argv[1]),cap=stoll(argv[2]);if(N<249||N>100000000||cap<0||cap>1000)throw runtime_error("declared range");
 vector<uint32_t>f(max(2*N/5,cap)+1);f[0]=1;for(size_t y=1;y<f.size();++y)f[y]=f[y/2]+(y%2?0:y/2+(y/2)%2);
 vector<unsigned char>hit((N/2+8)/8);uint64_t emitted=0,visits=0;
 auto emit=[&](I v){if(v<=N){if(v<1||!(v&1))throw runtime_error("parity");I i=v/2;hit[i/8]|=1<<(i%8);++emitted;}};
 function<void(I,I)>singles=[&](I y,I F){if(y>2*N)return;emit(F);if(y<=N){singles(2*y,F+y+(y%2));singles(2*y+1,F);}};
 emit(1);singles(1,1);
 for(I H:{5,25})for(I y=0;y<=2*N/H;++y)emit(H*f[y]);
 vector<array<I,4>> sources;
 for(int a=4;a<argc;++a){ifstream in(argv[a]);I H,E,p,q,s,t,F,mask;
  while(in>>H>>E>>p>>q>>s>>t>>F>>mask){
   if((H!=5&&H!=25)||p<0||q<0||(q>0&&q<5)||s<0||s>H||t<0||t>H||F<H)throw runtime_error("profile domain");
   for(I x=0;x<=cap;++x){I F0=H*f[x]+p*x+(H-s)*(x%2)+F-H,Q=E*x+q;
    if(Q>0&&F0<=N)sources.push_back({F0,Q,H-t,H});}
  }
  if(!in.eof())throw runtime_error("profile read");
 }
 sort(sources.begin(),sources.end());sources.erase(unique(sources.begin(),sources.end()),sources.end());
 for(auto b:sources){auto[F,Q,g,H]=b;for(I y=0;y<=(N-F)/Q;++y){++visits;emit(F+Q*y+g*(y%2)+H*(f[y]-1));}}
 vector<I>missing;I counts=0;for(I n=1;n<=N;n+=2){I i=n/2;if(hit[i/8]>>(i%8)&1)++counts;else missing.push_back(n);}
 ofstream out(argv[3],ios::binary);out.write((const char*)hit.data(),hit.size());if(!out)throw runtime_error("bitmap write");
 cout<<"{\"status\":\"passed\",\"limit\":"<<N<<",\"cap\":"<<cap<<",\"sources\":"<<sources.size()<<",\"pairs\":"<<visits<<",\"counts\":"<<counts<<",\"missing\":[";
 for(size_t i=0;i<missing.size();++i){if(i)cout<<",";cout<<missing[i];}cout<<"]}\n";
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
