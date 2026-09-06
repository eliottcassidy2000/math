// Exact exploratory atlas. No all-height conclusion follows from this scan.
// Universe: primitive a<b<c<=H, positive odd, all units modulo three.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using I=long long;
using V=std::array<I,3>;
struct Rat { I p=0,q=1; bool operator<(const Rat&r)const{return (__int128)p*r.q<(__int128)r.p*q;} };
struct Best {Rat val;V w{};I n=0,d=0; void put(Rat r,V v,I nn,I dd){if(val<r){val=r;w=v;n=nn;d=dd;}}};
V primitive(V v){I g=std::gcd(std::gcd(std::abs(v[0]),std::abs(v[1])),std::abs(v[2]));for(auto &x:v)x/=g;if(v[0]<0)for(auto &x:v)x=-x;return v;}
I mod(I x,I q){return (x%q+q)%q;}
I inv(I b,I m){for(I k=0;k<m;k++)if(b*k%m==1%m)return k;std::abort();}
void show(const char*s,const Best&b){I g=std::gcd(b.val.p,b.val.q);std::cout<<s<<" value="<<b.val.p/g<<"/"<<b.val.q/g<<" w="<<b.w[0]<<","<<b.w[1]<<","<<b.w[2]<<" N="<<b.n<<" directions="<<b.d<<"\n";}
int main(int argc,char**argv){
 int H=argc>1?std::atoi(argv[1]):199;I rows=0,multi=0,noCircuit=0,aboveCount=0,aboveMax=0;
 Best maxE,meanE,nOverC,minE;std::map<int,I> classes;V firstNo{};std::vector<V> firstPoints;
 std::map<V,std::pair<I,V>> circuits;
 for(I c=5;c<=H;c+=2)if(c%3)for(I b=3;b<c;b+=2)if(b%3){I g=std::gcd(b,c),m=c/g,iv=inv(b/g,m);
 for(I a=1;a<b;a+=2)if(a%3&&std::gcd(a,g)==1){
  rows++;V w{a,b,c};I Bx=(3*(b+c)-1)/14,By=(3*(a+c)-1)/14,Bz=(3*(a+b)-1)/14;
  std::vector<V> live;std::set<V> dirs;
  for(I x=-Bx;x<=Bx;x++)if(x%3&&a*x%g==0){I residue=mod(-(a*x/g)*iv,m);I y0=-By+mod(residue+By,m);
   for(I y=y0;y<=By;y+=m){I z=-(a*x+b*y)/c;if(y%3&&z%3&&std::abs(z)<=Bz){V v{x,y,z};live.push_back(v);dirs.insert(primitive(v));}}
  }
  I d=dirs.size(),n=live.size();classes[d]++;if(d<2)continue;multi++;
  if(d==3){std::vector<V> vs(dirs.begin(),dirs.end());V coeff{};
   for(int i=0;i<3;i++){V u=vs[(i+1)%3],v=vs[(i+2)%3];coeff[i]=std::abs(u[0]*v[1]-u[1]*v[0]);}
   I gg=std::gcd(std::gcd(coeff[0],coeff[1]),coeff[2]);for(auto &x:coeff)x/=gg;std::sort(coeff.begin(),coeff.end());
   auto &entry=circuits[coeff];entry.first++;if(entry.second[2]==0)entry.second=w;
  }
  I den=14*a*b*c,cap=6*a*b;V sums{};
  for(V v:live)for(int i=0;i<3;i++)sums[i]+=std::min(cap,w[i]*(3*(a+b+c-w[i])-14*std::abs(v[i])));
  maxE.put({*std::max_element(sums.begin(),sums.end()),den},w,n,d);
  minE.put({*std::min_element(sums.begin(),sums.end()),den},w,n,d);
  meanE.put({sums[0]+sums[1]+sums[2],3*den},w,n,d);nOverC.put({n,c},w,n,d);
  aboveCount+=11*n>2*c;aboveMax+=(__int128)*std::max_element(sums.begin(),sums.end())*77>6*(__int128)den;
  // A noncollinear additive circuit retains raw addresses, not only directions.
  if(d>=3){std::set<V> points(live.begin(),live.end());bool circuit=false;
   for(V u:live){for(V v:live){if(primitive(u)==primitive(v))continue;V z{u[0]+v[0],u[1]+v[1],u[2]+v[2]};if(points.count(z)){circuit=true;break;}}if(circuit)break;}
   if(!circuit){noCircuit++;if(firstNo[2]==0){firstNo=w;firstPoints=live;}}
  }
 }}
 std::cout<<"FINITE-EXACT H="<<H<<" rows="<<rows<<" multi="<<multi<<" count_gate_misses="<<aboveCount<<" every_projection_failures="<<aboveMax<<"\n";
 show("MAX_E",maxE);show("MAX_MIN_E",minE);show("MAX_MEAN_E",meanE);show("MAX_N_OVER_C",nOverC);
 std::cout<<"DIRECTION_COUNTS";for(auto [d,n]:classes)std::cout<<" "<<d<<":"<<n;std::cout<<"\n";
 std::cout<<"AT_LEAST_THREE_WITHOUT_ADDITIVE_CIRCUIT "<<noCircuit<<" first="<<firstNo[0]<<","<<firstNo[1]<<","<<firstNo[2]<<"\n";
 for(V v:firstPoints)std::cout<<"FIRST_NO_CIRCUIT_POINT "<<v[0]<<","<<v[1]<<","<<v[2]<<"\n";
 for(auto [coeff,data]:circuits)std::cout<<"THREE_RAY_CIRCUIT "<<coeff[0]<<","<<coeff[1]<<","<<coeff[2]<<" rows="<<data.first<<" first="<<data.second[0]<<","<<data.second[1]<<","<<data.second[2]<<"\n";
}
