// Independent arithmetic image via the two-integer fixed-block formula.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=int64_t;
int main(int argc,char**argv){try{
 if(argc<4)throw runtime_error("usage: scan ceiling output-bitmap profiles...");
 I N=stoll(argv[1]);if(N<1||N>10000000)throw runtime_error("ceiling out of checked integer range");
 vector<uint32_t>D(2*N+1);for(I x=1;x<=2*N;++x)D[x]=D[x/2]+(x%2?0:(x/2+1)/2);
 vector<unsigned char>seen(N+1);vector<uint32_t> least(N+1,UINT32_MAX);I visits=0,profiles=0;
 auto emit=[&](I F,I x=0){if(F<=N){if(F<1||!(F&1))throw runtime_error("invalid count");seen[F]=1;least[F]=min(least[F],(uint32_t)x);}};
 for(I x=0;x<=2*N;++x)emit(1+2*D[x]);
 for(int arg=3;arg<argc;++arg){ifstream in(argv[arg]);I H,E,p,q,s,t,F,mask;
  while(in>>H>>E>>p>>q>>s>>t>>F>>mask){
   ++profiles;if(H<1||E<2||p<0||q<0||s<0||s>H||t<0||t>H||F<H)throw runtime_error("invalid profile");
   for(I x=0;x<=2*N/H;++x)emit(H*(1+2*D[x]));
   auto value=[&](I x,I y){++visits;return F+(p+E*y)*x+q*y+(H-s)*(x%2)+(H-t)*(y%2)+2*H*(D[x]+D[y]);};
   if(q)for(I y=0;y<=(N-F)/q;++y)emit(value(0,y));
   I maxx=p?(N-F)/p:(N-F-q)/E;
   for(I x=1;x<=maxx;++x){I ymax=(N-F-p*x)/(E*x+q);
    for(I y=(p?0:1);y<=ymax;++y)emit(value(x,y),x);
   }
  }
  if(!in.eof())throw runtime_error("invalid profile file");
 }
 I count=0;vector<I>missing;for(I F=1;F<=N;F+=2){count+=seen[F];if(!seen[F]&&F!=7&&F!=21)missing.push_back(F);}
 if((N>=7&&seen[7])||(N>=21&&seen[21]))throw runtime_error("forbidden count");
 ofstream out(argv[2],ios::binary);out.write((const char*)seen.data(),seen.size());if(!out)throw runtime_error("write failed");
 cout<<"{\"prefix_caps\":[";bool first=true;
 for(I cap:vector<I>{0,1,3,7,15,31,63,127,255,511,1023,4095,16383}){if(!first)cout<<",";first=false;I count=0;vector<I>first_missing;
 for(I z=1;z<=N;z+=2){count+=least[z]<=cap;if(least[z]>cap&&z!=7&&z!=21&&first_missing.size()<16)first_missing.push_back(z);}
 cout<<"{\"cap\":"<<cap<<",\"counts\":"<<count<<",\"first_missing\":[";for(size_t i=0;i<first_missing.size();++i){if(i)cout<<",";cout<<first_missing[i];}cout<<"]}";}
 cout<<"],\"status\":\"passed\",\"limit\":"<<N<<",\"profiles\":"<<profiles<<",\"integer_pairs\":"<<visits<<",\"counts\":"<<count<<",\"allowed_missing\":[";
 for(size_t i=0;i<missing.size();++i){if(i)cout<<",";cout<<missing[i];}cout<<"]}\n";
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
