// Independent depth-first complete core traversal and DIRECT outer-word test.
// No includes from the parent's enumerator or binary inverse.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using I=long long;
using W=__int128_t;
using State=std::array<W,6>;
struct Word {I bits; int n; std::array<I,9> c;};
I N,limit; int cap; uint64_t direct_tests=0,hits=0;
std::vector<uint32_t> shell;
std::vector<std::vector<Word>> words;
std::vector<uint64_t> orders,attempts;
std::vector<std::string> address;
std::ostream *out;
I power(I a,int b){I r=1;while(b--)r*=a;return r;}
Word compile(I bits,int n){
  Word w{bits,n,{}};
  for(int m=1;m<=3;++m){
    I L=power(m-1,n),D=L;
    std::array<I,6> f{1,0,0,0,0,0};
    for(int j=0;j<n;++j){
      int digit=(bits>>j)&1;
      if(digit)L+=power(m,j)*power(m-1,n-j-1);
      else D+=power(m-1,j)*power(m,n-j-1);
      std::array<I,6> next{};
      for(int old=0;old<6;++old){
        int last=old/2,started=old%2;
        next[old]+=(m-1)*f[old];
        next[2*(digit?1:2)+started]+=f[old];
        if(started||!digit)next[2*last+1]+=f[old];
      }
      f=next;
    }
    w.c[3*(m-1)]=L;w.c[3*(m-1)+1]=D;
    w.c[3*(m-1)+2]=f[0]+f[1]+f[2]+f[3];
  }
  return w;
}
State child(State s,const Word&w){
  const auto &v=w.c; W h=s[0],e=s[1],a=s[2],p=s[3],f=s[4],b=s[5];
  return {v[3]*a+v[2]*h-v[1]*e,
          v[0]*h-(w.n?0:e),
          power(2,w.n)*a+v[4]*h,
          2*(v[6]*b+v[5]*p-v[4]*f),
          v[3]*p-f,
          power(3,w.n)*b+v[7]*p};
}
std::string spelling(I word,int n){std::string s;for(int i=0;i<n;++i)s+=char('0'+((word>>i)&1));return s;}
void visit(const State&s,int n){
  ++orders[n];
  I h=I(s[0]),e=I(s[1]),G=I(s[2]+s[3]),d=I(s[2]+s[3]-s[1]-s[4]-1);
  if(N>=h+d){
    I xmax=(N-h-d)/G;
    if(xmax>=I(shell.size())){std::cerr<<"shell-bound failure\n";std::exit(5);}
    for(I x=(n==1?1:0);x<=xmax;++x){
      ++direct_tests;
      if(W(h)*shell[x]+W(G)*x+e*(x&1)+d!=N)continue;
      if(hits++)*out<<",\n";
      *out<<"{\"core_order\":"<<n<<",\"outer_word\":"<<x<<",\"address\":[";
      for(size_t j=0;j<address.size();++j){if(j)*out<<",";*out<<"\""<<address[j]<<"\"";}
      *out<<"],\"state\":[";for(int j=0;j<6;++j){if(j)*out<<",";*out<<I(s[j]);}*out<<"]}";
    }
  }
  for(int t=0;n+t+1<=cap;++t){
    if(n==1&&t==0)continue;
    W lower=t?W(power(2,t))+power(3,t-1)+1:2;
    if(lower*s[3]>limit)break;
    if(words[t].empty()){
      I count=t?(I(1)<<(t-1)):1;
      for(I j=0;j<count;++j)words[t].push_back(compile(t?j+(I(1)<<(t-1)):0,t));
    }
    for(const auto&w:words[t]){
      ++attempts[n+t+1];State z=child(s,w);
      if(z[3]>limit)continue;
      if(z[3]<2*s[3])std::exit(6);
      address.push_back(spelling(w.bits,t));
      visit(z,n+t+1);
      address.pop_back();
    }
  }
}
int main(int argc,char**argv){
  if(argc!=3)return 2;N=std::stoll(argv[1]);if(N<3||N>200000000)return 3;
  limit=2*N;cap=0;while((I(1)<<(cap+1))<=limit)++cap;
  shell.resize(N/3+1);shell[0]=1;
  for(size_t x=1;x<shell.size();++x){size_t q=x/2;shell[x]=(x&1)?shell[q]:shell[q]+q+(q&1);}
  words.resize(cap+1);orders.assign(cap+1,0);attempts.assign(cap+1,0);
  std::ofstream receipt(argv[2]);out=&receipt;
  receipt<<"{\"target\":"<<N<<",\"status\":\"complete independent depth-first core traversal with direct outer-word enumeration\",\"hits\":[\n";
  visit({1,1,1,2,1,1},1);
  receipt<<"\n],\"layers\":[";
  uint64_t total=0;for(int n=1;n<=cap;++n){if(n>1)receipt<<",";receipt<<"{\"order\":"<<n<<",\"retained\":"<<orders[n]<<",\"attempts\":"<<attempts[n]<<"}";total+=orders[n];}
  receipt<<"],\"total_cores\":"<<total<<",\"direct_outer_tests\":"<<direct_tests<<",\"found\":"<<(hits?"true":"false")<<"}\n";
  std::cout<<"target="<<N<<" cores="<<total<<" direct_outer_tests="<<direct_tests<<" hits="<<hits<<"\n";
}
