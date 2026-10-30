// Exact Bellman certificate for truncated five-adic content, all valuation letters.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

void need(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
int index_of(int m,int b,int q){return (m-m/5-1)*q+b;}
int main(int argc,char**argv){
  const std::string prefix=argc>1?argv[1]:"/tmp/five-potential";
  std::cout<<"[\n";
  int q=1;
  for(int E=1;E<=4;++E){
    q*=5;int num=q-1,den=8*(q/5),g=std::gcd(num,den);num/=g;den/=g;
    int count=q*(q-q/5),cap=2*E+4;
    std::vector<int> reward(count),h(count),ms(count),bs(count);
    for(int m=1;m<q;++m)if(m%5)for(int b=0;b<q;++b){
      int i=index_of(m,b,q);ms[i]=m;bs[i]=b;
      for(int Q=5;Q<=q;Q*=5)reward[i]+=((m-1)%Q==0 && b%Q==0);
    }
    std::vector<std::vector<int>> edges(cap,std::vector<int>(count));
    int inv=1;
    for(int a=1;a<=cap;++a){
      inv=inv*((q+1)/2)%q;
      for(int i=0;i<count;++i)edges[a-1][i]=index_of(3*inv*ms[i]%q,(3*bs[i]+1)*inv%q,q);
    }
    int rounds=0;bool change;
    do{
      change=false;++rounds;
      for(int a=1;a<=cap;++a)for(int i=0;i<count;++i){
        int j=edges[a-1][i],v=h[i]+den*reward[j]-num*a;
        if(v>h[j]){h[j]=v;change=true;}
      }
      need(rounds<=count,"potential search did not terminate");
    }while(change);
    int maximum=*std::max_element(h.begin(),h.end());
    need(*std::min_element(h.begin(),h.end())==0,"potential minimum");
    need(h[index_of(1,0,q)]==maximum,"identity must maximize potential");
    // Beyond the checked menu the cost exceeds every possible gain, without
    // identifying or resetting the native integer clock modulo a period.
    need(num*(cap+1)>=maximum+den*E,"unbounded-letter tail certificate");
    std::uint64_t inequalities=0;
    for(int a=1;a<=cap;++a)for(int i=0;i<count;++i){
      int j=edges[a-1][i];need(h[j]>=h[i]+den*reward[j]-num*a,"Bellman inequality");
      ++inequalities;
    }
    // Separate modular trajectory verifies equality and hence optimality.
    int m=1,b=0,total=0,period=2*(q/5),iv=(q+1)/2;
    for(int t=0;t<period;++t){m=3*iv*m%q;b=(3*b+1)*iv%q;total+=reward[index_of(m,b,q)];}
    need(m==1 && b==0 && den*total==num*period,"constant-one optimality");
    std::ofstream f(prefix+"-"+std::to_string(E)+".txt");
    need(bool(f),"potential output");for(int v:h)f<<v<<'\n';
    std::cout<<"  {\"E\":"<<E<<",\"q\":"<<q<<",\"states\":"<<count
      <<",\"numerator\":"<<num<<",\"denominator\":"<<den<<",\"maximum\":"<<maximum
      <<",\"checked_letters\":"<<cap<<",\"inequalities\":"<<inequalities
      <<",\"rounds\":"<<rounds<<",\"equality_period\":"<<period<<"}"<<(E==4?"\n":",\n");
  }
  std::cout<<"]\n";
}
