#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;using Z=long long;
struct G{int h,q,b;};
struct S{int h,e,r;bool operator<(S const&x)const{return r!=x.r?r<x.r:h!=x.h?h<x.h:e<x.e;}};
void need(bool x,const char*s){if(!x)throw runtime_error(s);}
G literal_bits(int L){if(L==1)return{1,1,1};unsigned z=L-1;int m=32-__builtin_clz(z),g=1,Q=1,W=0,previous=-1;for(int i=0;i<m;i++){int w=1-((z>>i)&1);if(i==0)g+=w;if(!w)g+=W;if(previous==0&&w==1)g++;if(w)Q+=1<<(m-1-i);W=2*W+w;previous=w;}return{g,Q,1<<m};}
int main(int argc,char**argv){need(argc==3,"output prefix and bound required");string prefix=argv[1];int N=stoi(argv[2]);need(N>=3&&N<=10000001,"bound");vector<G>g(N+1);for(int L=1;L<=N;L++)g[L]=literal_bits(L);cout<<"{\"bound\":"<<N<<",\"images\":[";bool first=true;for(int B:{500,1000}){map<int,set<S>>byR;byR[1].insert({1,1,1});vector<S>bank;for(auto it=byR.begin();it!=byR.end();++it){for(S s:it->second){bank.push_back(s);for(int L=1;L<=B/s.r;L++){auto a=g[L];Z rr=Z(a.q)*s.h+Z(a.b)*s.r;if(rr>B)continue;S out{a.h*s.h-(L&1)*s.e+L*s.r,s.h-(L==1)*s.e,int(rr)};need(out.r>s.r&&out.h<=out.r,"monotone prefix response");need(out.h<=out.r&&out.r<=2*(out.h-out.e),"output cone");byR[out.r].insert(out);}}}need(bank.size()==(B==500?2717:7497),"prefix universe");vector<uint64_t>packed((N+64)/64);auto mark=[&](Z h){need(h>=1&&(h&1),"physical odd value");if(h<=N)packed[h/64]|=uint64_t(1)<<(h%64);};mark(1);Z trials=0;for(auto s:bank){int limit=N/s.r;for(int L=1;L<=limit;L++){trials++;auto a=g[L];mark(Z(a.h)*s.h-(L&1)*s.e+Z(L)*s.r);}}vector<int>missing;for(int h=3;h<=N;h+=2)if(!(packed[h/64]>>(h%64)&1))missing.push_back(h);vector<unsigned char>bytes((N+8)/8);for(size_t i=0;i<bytes.size();i++)bytes[i]=(packed[i/8]>>(8*(i%8)))&255;ofstream f(prefix+"-"+to_string(B)+".bin",ios::binary);f.write((char*)bytes.data(),bytes.size());need(bool(f),"image write");ofstream bfile(prefix+"-"+to_string(B)+"-prefixes.csv");for(auto s:bank)bfile<<s.h<<','<<s.e<<','<<s.r<<'\n';if(!first)cout<<',';first=false;cout<<"{\"prefix_R_bound\":"<<B<<",\"prefix_states\":"<<bank.size()<<",\"terminal_trials\":"<<trials<<",\"missing_odd\":[";for(int i=0;i<(int)missing.size();i++){if(i)cout<<',';cout<<missing[i];}cout<<"]}";}cout<<"],\"status\":\"PASS\"}\n";}
