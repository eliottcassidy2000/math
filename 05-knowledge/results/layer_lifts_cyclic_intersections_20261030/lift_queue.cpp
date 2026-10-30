// Sparse monotone relaxation; the certificate is the final edge inequalities.
#include <algorithm>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using I=std::int64_t;
void need(bool b,const char*msg){if(!b)throw std::runtime_error(msg);}
int idx(int m,int b,int q){return (m-m/5-1)*q+b;}
int reward(int m,int b,int E){int r=0;for(int j=0;j<E && (m-1)%5==0 && b%5==0;++j){++r;m=(m-1)/5+1;b/=5;}return r;}
int main(int argc,char**argv){
 int E=argc>1?std::stoi(argv[1]):5,q=1;need(E>=1 && E<=5,"this producer supports audited levels 1 through 5");for(int i=0;i<E;++i)q*=5;
 std::string out=argc>2?argv[2]:"/tmp/lift-queue.bin";
 int N=q-1,D=8*(q/5),g=std::gcd(N,D);N/=g;D/=g;
 int S=q*(q-q/5),cap=4*E;std::vector<int>inv(cap+1,1);
 for(int a=1;a<=cap;++a)inv[a]=I(inv[a-1])*((q+1)/2)%q;
 std::vector<std::int32_t>H(S,0);std::vector<unsigned char>queued(S,0);std::deque<int>Q;
 for(int m=1;m<q;++m)if(m%5)for(int b=0;b<q;++b){int r=reward(m,b,E);if(r){int j=idx(m,b,q);H[j]=D*r-N;queued[j]=1;Q.push_back(j);}}
 I popped=0,relaxations=0;
 while(!Q.empty()){
   int i=Q.front();Q.pop_front();queued[i]=0;++popped;
   int mi=i/q,m=mi+1+mi/4,b=i%q;
   for(int a=1;a<=cap;++a){
     if(H[i]+D*E-N*a<=0)break;
     int mm=3*I(inv[a])*m%q,bb=(3*I(b)+1)*inv[a]%q,j=idx(mm,bb,q);
     int candidate=H[i]+D*reward(mm,bb,E)-N*a;
     if(candidate>H[j]){H[j]=candidate;++relaxations;if(!queued[j]){queued[j]=1;Q.push_back(j);}}
     need(H[j]<1000*D*E,"potential growth guard: investigate a positive cycle");
   }
 }
 int maximum=*std::max_element(H.begin(),H.end());need(H[0]==maximum,"identity maximum");
 need(N*(cap+1)>=maximum+D*E,"unbounded letter tail");
 I checked=0;
 for(int m=1;m<q;++m)if(m%5)for(int b=0;b<q;++b){int i=idx(m,b,q);
  for(int a=1;a<=cap;++a){int mm=3*I(inv[a])*m%q,bb=(3*I(b)+1)*inv[a]%q,j=idx(mm,bb,q);
   need(H[j]>=H[i]+D*reward(mm,bb,E)-N*a,"Bellman inequality");++checked;}
 }
 std::ofstream f(out,std::ios::binary);need(bool(f),"output");
 // Explicit little-endian bytes, independently decoded by the Python auditor.
 for(auto x:H){for(int j=0;j<4;++j){unsigned char c=(std::uint32_t(x)>>(8*j))&255;f.put(char(c));}}
 need(bool(f),"complete certificate write");
 std::cout<<"{\"E\":"<<E<<",\"q\":"<<q<<",\"states\":"<<S<<",\"N\":"<<N<<",\"D\":"<<D
 <<",\"maximum\":"<<maximum<<",\"letters\":"<<cap<<",\"popped\":"<<popped<<",\"relaxations\":"<<relaxations<<",\"inequalities\":"<<checked<<"}\n";
}
