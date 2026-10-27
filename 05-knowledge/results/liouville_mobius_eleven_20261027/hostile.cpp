// Extension of clock_even_guarded_composition_20261026/hostile.cpp.
// Same complete path/presentation universe; new mu(H) and lambda(H) readouts.
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using U=uint64_t;using T=std::array<unsigned,12>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
std::pair<int,int> signs(U n){int mu=1,lambda=1;
 for(U p=2;p*p<=n;p++)if(n%p==0){int e=0;while(n%p==0){n/=p;e++;lambda=-lambda;}mu=e>1?0:-mu;}
 if(n>1){mu=-mu;lambda=-lambda;}return {mu,lambda};}
U paths(const T&t){std::vector<std::array<U,12>>dp(4096);
 for(int v=0;v<12;v++)dp[1<<v][v]=1;
 for(unsigned s=1;s<4096;s++)for(int v=0;v<12;v++)if(dp[s][v]){
  unsigned avail=t[v]&(4095^s);while(avail){unsigned b=avail&-avail;avail-=b;int w=__builtin_ctz(b);dp[s|b][w]+=dp[s][v];}}
 U h=0;for(auto x:dp[4095])h+=x;return h;}
T build(bool minus){T t{};for(int b=0;b<3;b++){
 for(int i=0;i<3;i++){t[4*b+i]|=1u<<(4*b+(i+1)%3);if(minus&&b==2)t[4*b+3]|=1u<<(4*b+i);else t[4*b+i]|=1u<<(4*b+3);}
 for(int i=0;i<4;i++)t[4*b+i]|=15u<<(4*((b+1)%3));}return t;}
int main(){std::cout<<"[\n";for(int minus=0;minus<2;minus++){
 auto t=build(minus);std::array<int,12>p{};std::set<U>codes;U raw=0;
 auto dfs=[&](auto&&self,unsigned used,int len)->void{
  if(len==12){raw++;U code=0;int k=0;for(int i=0;i<12;i++)for(int j=i+2;j<12;j++)code|=U((t[p[i]]>>p[j])&1)<<k++;codes.insert(code);return;}
  unsigned avail=(len?t[p[len-1]]:4095)&(4095^used);while(avail){unsigned b=avail&-avail;avail-=b;p[len]=__builtin_ctz(b);self(self,used|b,len+1);}};
 dfs(dfs,0,0);unsigned aut=minus?27:81;need(raw==paths(t)&&raw==462915&&codes.size()*aut==raw,"source/presentation control");
 std::array<U,3>mus{};std::array<U,2>lams{};U sum=0;std::vector<std::array<U,2>>hsamples;std::map<U,U>hist;
 for(U c:codes){T partner{};int k=0;for(int i=0;i<12;i++)for(int j=i+1;j<12;j++){
  bool forward=j==i+1||!((c>>k)&1);if(j>i+1)k++;if(forward)partner[i]|=1u<<j;else partner[j]|=1u<<i;}
  U h=paths(partner);sum+=h;auto [mu,lambda]=signs(h);mus[mu+1]+=aut;lams[(lambda+1)/2]+=aut;
  hist[h]+=aut;if(hsamples.size()<8)hsamples.push_back({c,h});}
 need(mus[0]+mus[1]+mus[2]==raw&&lams[0]+lams[1]==raw,"readout mass");
 if(minus)std::cout<<",\n";
 std::cout<<"{\"word\":\""<<(minus?"001":"000")<<"\",\"Aut\":"<<aut<<",\"H\":"<<raw<<",\"presentations\":"<<codes.size()
 <<",\"L2_one\":"<<sum*aut<<",\"mu_hist_minus_zero_plus\":["<<mus[0]<<","<<mus[1]<<","<<mus[2]
 <<"],\"lambda_hist_minus_plus\":["<<lams[0]<<","<<lams[1]<<"],\"L_muH\":"<<int64_t(mus[2])-int64_t(mus[0])
 <<",\"L_lambdaH\":"<<int64_t(lams[1])-int64_t(lams[0])<<",\"partner_H_histogram\":[";
 bool first=true;for(auto [h,m]:hist){if(!first)std::cout<<",";first=false;std::cout<<"["<<h<<","<<m<<"]";}
 std::cout<<"],\"sample\":[";first=true;for(auto a:hsamples){if(!first)std::cout<<",";first=false;std::cout<<"["<<a[0]<<","<<a[1]<<"]";}
 std::cout<<"]}";std::cout.flush();}
 std::cout<<"\n]\n";}
