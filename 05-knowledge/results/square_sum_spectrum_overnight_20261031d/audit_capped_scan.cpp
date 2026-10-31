// Independent quadratic-form image; f is updated by a binary carry identity.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=int64_t;
I next_f(I y,I old){ // y is the new argument, so y>=1.
 int k=__builtin_ctzll((uint64_t)y);I a=y>>(k+1);
 return old+(k?y-(y>>k)+1:0)-a-(a%2);
}
I from_bits(I y){I F=1;while(y){I a=y/2;if(!(y&1))F+=a+a%2;y=a;}return F;}
int main(int argc,char**argv){try{
 if(argc<5)throw runtime_error("usage: audit limit cap bitmap profile-files...");
 I N=stoll(argv[1]),cap=stoll(argv[2]);if(N<249||N>100000000||cap<0||cap>1000)throw runtime_error("declared range");
 for(I y=1,F=1;y<=1000000;++y){F=next_f(y,F);if(F!=from_bits(y))throw runtime_error("carry control");}
 vector<unsigned char>hit((N/2+8)/8);I visits=0;
 auto emit=[&](I v){if(v<=N){if(v<1||!(v&1))throw runtime_error("parity");I i=v/2;hit[i/8]|=1<<(i%8);}};
 for(I H:{1,5,25}){I fy=1;emit(H);for(I y=1;y<=2*N/H;++y){fy=next_f(y,fy);emit(H*fy);}}
 for(int a=4;a<argc;++a){ifstream in(argv[a]);I H,E,p,q,s,t,F,mask;
  while(in>>H>>E>>p>>q>>s>>t>>F>>mask){
   for(I x=0;x<=cap;++x){I Q=E*x+q;if(Q==0)continue;
    I fx=from_bits(x),fy=1,last=(N-F-p*x)/Q;
    if(N<F+p*x)continue;
    for(I y=0;y<=last;++y){if(y)fy=next_f(y,fy);++visits;
     emit(F+p*x+Q*y+(H-s)*(x%2)+(H-t)*(y%2)+H*(fx+fy-2));}
   }
  }
  if(!in.eof())throw runtime_error("profile read");
 }
 I count=0;vector<I>missing;for(I n=1;n<=N;n+=2){I i=n/2;if(hit[i/8]>>(i%8)&1)++count;else missing.push_back(n);}
 ofstream out(argv[3],ios::binary);out.write((const char*)hit.data(),hit.size());if(!out)throw runtime_error("write");
 cout<<"{\"status\":\"passed\",\"limit\":"<<N<<",\"cap\":"<<cap<<",\"pairs\":"<<visits<<",\"counts\":"<<count<<",\"carry_controls\":1000000,\"missing\":[";
 for(size_t i=0;i<missing.size();++i){if(i)cout<<",";cout<<missing[i];}cout<<"]}\n";
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
