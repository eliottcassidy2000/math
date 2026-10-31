// Independent source-state enumeration, using the direct binary-inversion f.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=int64_t;
#ifndef SOURCE_CAP
#define SOURCE_CAP 10000000
#endif
int main(int argc,char**argv){try{
 if(argc<4)throw runtime_error("usage: audit ceiling bitmap profiles...");
 I N=stoll(argv[1]);if(N<21||N>10000000)throw runtime_error("declared range");
 vector<uint32_t>f(2*N+1);f[0]=1;
 for(I x=1;x<=2*N;++x){I a=x,value=1;while(a>1){if(!(a&1))value+=(a>>1)+((a>>1)&1);a>>=1;}f[x]=uint32_t(value);}
 vector<unsigned char>hit(N+1);I evaluations=0,profiles=0;
 auto emit=[&](I v){if(v>=1&&v<=N){if(!(v&1))throw runtime_error("even output");hit[v]=1;}};
 for(I x=0;x<=2*N;++x)emit(f[x]);
 for(int a=3;a<argc;++a){ifstream in(argv[a]);I H,E,p,q,s,t,F,mask;
  while(in>>H>>E>>p>>q>>s>>t>>F>>mask){
   ++profiles;for(I x=0;x<=2*N/H;++x)emit(H*f[x]);
   I last=p?(N-F)/p:(N-F-q)/E;last=min(last,I(SOURCE_CAP));
   for(I x=0;x<=last;++x){
    I initial=H*f[x]+p*x+(H-s)*(x&1)+F-H;
    I Q=E*x+q,G=H-t;if(Q==0||initial>N)continue;
    for(I y=(p==0&&x>0?1:0);y<=(N-initial)/Q;++y){
     ++evaluations;emit(initial+Q*y+G*(y&1)+H*(f[y]-1));
    }
   }
  }
  if(!in.eof())throw runtime_error("malformed profile input");
 }
 if(hit[7]||hit[21])throw runtime_error("known exclusion");
 I total=0;for(auto b:hit)total+=b;
 ofstream out(argv[2],ios::binary);out.write((const char*)hit.data(),hit.size());if(!out)throw runtime_error("output failed");
 cout<<"{\"status\":\"passed\",\"limit\":"<<N<<",\"source_cap\":"<<SOURCE_CAP<<",\"profiles\":"<<profiles<<",\"digital_evaluations\":"<<evaluations<<",\"counts\":"<<total<<"}\n";
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
