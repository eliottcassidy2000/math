#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <map>
#include <set>
#include <vector>
using U=uint64_t;using T=std::vector<unsigned>;
U count(const T&r){
 unsigned n=r.size(),full=(1u<<n)-1;std::vector<U>d((full+1)*n);
 for(unsigned v=0;v<n;v++)d[(1u<<v)*n+v]=1;
 for(unsigned S=1;S<=full;S++)for(unsigned v=0;v<n;v++)if(d[S*n+v]){
  unsigned next=r[v]&~S&full;while(next){unsigned w=__builtin_ctz(next);next&=next-1;d[(S|(1u<<w))*n+w]+=d[S*n+v];}
 }
 U a=0;for(unsigned v=0;v<n;v++)a+=d[full*n+v];return a;
}
U literal(const T&r){std::vector<unsigned>v(r.size());std::iota(v.begin(),v.end(),0);U total=0;do{bool good=true;for(unsigned i=1;i<v.size();i++)if(!(r[v[i-1]]>>v[i]&1)){good=false;break;}total+=good;}while(std::next_permutation(v.begin(),v.end()));return total;}
struct Repair{U H,mid,target;int x,y;T rows;};
T ear(T r,unsigned mask){unsigned n=r.size();for(unsigned u=0;u<n;u++)if(!(mask>>u&1))r[u]|=1u<<n;r.push_back(mask);return r;}
T birth(T r,int x){return ear(r,1u<<x);}
int main(){int n;U H;unsigned cases=0,children=0,perms=0,repair_checks=0,pair_checks=0,ear_checks=0;std::map<std::pair<U,int>,Repair> repairs,pair_repairs,ear_repairs;std::set<U> dead={31,43,55,61,63,73,85,87,113};while(std::cin>>n>>H){T r(n);for(auto &x:r)std::cin>>x;std::vector<long long>I(n);for(auto &x:I)std::cin>>x;
 if(count(r)!=H)return 2;if(n<=7){if(literal(r)!=H)return 3;perms++;}
 for(int x=0;x<n;x++){T b=birth(r,x);U mid=count(b);if((long long)mid-2*(long long)H!=I[x])return 4;children++;
  if(dead.count(H)&&mid<=2*H+1)for(int y=0;y<=n;y++){
   U target=count(birth(b,y));repair_checks++;
   if(target==2*H-1||target==2*H+1){int sign=target==2*H-1?-1:1;auto key=std::make_pair(H,sign);if(!repairs.count(key))repairs[key]={H,mid,target,x,y,r};}
  }
 }
 if(dead.count(H))for(int x=0;x<n;x++)for(int y=x+1;y<n;y++){
  U target=count(ear(r,(1u<<x)|(1u<<y)));pair_checks++;
  if(target==2*H-1||target==2*H+1){int sign=target==2*H-1?-1:1;auto key=std::make_pair(H,sign);if(!pair_repairs.count(key))pair_repairs[key]={H,0,target,x,y,r};}
 }
 if(dead.count(H))for(unsigned mask=1;mask<(1u<<n)-1;mask++){
  int size=__builtin_popcount(mask);if(size<3)continue;
  U target=count(ear(r,mask));ear_checks++;
  if(target==2*H-1||target==2*H+1){int sign=target==2*H-1?-1:1;auto key=std::make_pair(H,sign);if(!ear_repairs.count(key)||size<ear_repairs[key].y)ear_repairs[key]={H,0,target,int(mask),size,r};}
 }
 cases++;}
 std::cout<<"{\"status\":\"PASS\",\"source_checks\":"<<cases<<",\"literal_permutation_checks\":"<<perms<<",\"full_adjacency_child_checks\":"<<children<<",\"two_birth_checks\":"<<repair_checks<<",\"repairs\":[";
 bool comma=false;for(auto pair:repairs){auto a=pair.second;if(comma)std::cout<<',';comma=true;std::cout<<"{\"H\":"<<a.H<<",\"intermediate\":"<<a.mid<<",\"target\":"<<a.target<<",\"ports\":["<<a.x<<','<<a.y<<"],\"rows\":[";for(unsigned i=0;i<a.rows.size();i++){if(i)std::cout<<',';std::cout<<a.rows[i];}std::cout<<"]}";}
 std::cout<<"],\"two_port_checks\":"<<pair_checks<<",\"two_port_repairs\":[";comma=false;
 for(auto pair:pair_repairs){auto a=pair.second;if(comma)std::cout<<',';comma=true;std::cout<<"{\"H\":"<<a.H<<",\"target\":"<<a.target<<",\"ports\":["<<a.x<<','<<a.y<<"],\"rows\":[";for(unsigned i=0;i<a.rows.size();i++){if(i)std::cout<<',';std::cout<<a.rows[i];}std::cout<<"]}";}
 std::cout<<"],\"larger_ear_checks\":"<<ear_checks<<",\"larger_ear_repairs\":[";comma=false;
 for(auto pair:ear_repairs){auto a=pair.second;if(comma)std::cout<<',';comma=true;std::cout<<"{\"H\":"<<a.H<<",\"target\":"<<a.target<<",\"mask\":"<<a.x<<",\"port_count\":"<<a.y<<",\"rows\":[";for(unsigned i=0;i<a.rows.size();i++){if(i)std::cout<<',';std::cout<<a.rows[i];}std::cout<<"]}";}
 std::cout<<"]}\n";return 0;
}
