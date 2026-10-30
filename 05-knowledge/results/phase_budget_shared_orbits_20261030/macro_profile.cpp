// Independent modular orbit enumeration; compares against Python arrival-address formulas.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I=std::int64_t;
void need(bool b,const char*m){if(!b)throw std::runtime_error(m);}
int main(int argc,char**argv){
 need(argc==3,"case and binary path");std::string name=argv[1];std::vector<int>u;
 int K;
 if(name=="U3"){u={3,1,3};u.insert(u.end(),6,1);u.insert(u.end(),{3,1,3});K=5;}
 else if(name=="W5"){std::vector<int>r={5};r.insert(r.end(),8,1);r.push_back(5);u=r;u.insert(u.end(),233,1);u.insert(u.end(),r.begin(),r.end());K=9;}
 else throw std::runtime_error("unknown family");
 I q=1;for(int e=0;e<K;++e)q*=5;I T=2*q/5;
 std::vector<I> inv(6,1),two(6,1);
 for(int a=1;a<=5;++a){inv[a]=inv[a-1]*((q+1)/2)%q;two[a]=two[a-1]*2%q;}
 struct State{I m,z;};
 auto step=[&](State s,int a){return State{3*s.m*inv[a]%q,((3*s.z+two[a]-2+q)*inv[a])%q};};
 auto val=[&](I x){int h=0;while(h<K&&x%5==0){++h;x/=5;}return h;};
 auto reward=[&](State s){return std::min(val((s.m+q-1)%q),val(s.z));};
 State s{1,0};I initial=0;for(int a:u){s=step(s,a);initial+=reward(s);}
 State entry=s;int h=val(s.z);I powh=1;for(int i=0;i<h;++i)powh*=5;
 I gain=(s.m%5==1||s.m%5==4)?T*5*(powh-1)/(8*powh):0;
 I middle=0,best=-(I(1)<<60),where=-1,ties=0;std::ofstream f(argv[2],std::ios::binary);need(bool(f),"open");
 for(I L=0;L<T;++L){State t=s;I total=initial+middle;for(int a:u){t=step(t,a);total+=reward(t);}
  I score=T*total-gain*L;
  if(score>best){best=score;where=L;ties=1;}else if(score==best)++ties;
  need(total>=0&&total<2147483648LL,"int32 profile");
  for(int j=0;j<4;++j)f.put(char((std::uint32_t(total)>>(8*j))&255));
  s=step(s,1);middle+=reward(s);
 }
 need(s.m==entry.m&&s.z==entry.z&&middle==gain,"exact period and gain");need(bool(f),"complete output");
 std::cout<<"{\"K\":"<<K<<",\"period\":"<<T<<",\"gain\":"<<gain<<",\"max_numerator\":"<<best<<",\"first_maximizer\":"<<where<<",\"maximizers\":"<<ties<<"}\n";
}
