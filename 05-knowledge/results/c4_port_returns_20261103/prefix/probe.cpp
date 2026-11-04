#include <algorithm>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using T=std::vector<unsigned>;
using I=std::int64_t;

unsigned canonical(const T& t){
 int n=t.size();std::vector<int> v(n);for(int i=0;i<n;i++)v[i]=i;
 std::sort(v.begin(),v.end(),[&](int a,int b){return std::popcount(t[a])<std::popcount(t[b])||(std::popcount(t[a])==std::popcount(t[b])&&a<b);});
 std::vector<std::pair<int,int>> groups;
 for(int i=0;i<n;){int j=i+1;while(j<n&&std::popcount(t[v[j]])==std::popcount(t[v[i]]))++j;groups.push_back({i,j});i=j;}
 unsigned best=~0u;
 std::function<void(int)> visit=[&](int at){
  if(at==(int)groups.size()){unsigned code=0;int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(t[v[i]]>>v[j]&1)code|=1u<<k;best=std::min(best,code);return;}
  auto[a,b]=groups[at];std::sort(v.begin()+a,v.begin()+b);
  do{visit(at+1);}while(std::next_permutation(v.begin()+a,v.begin()+b));
 };visit(0);return best;
}
bool strong(const T&t){
 int n=t.size();for(int start=0;start<n;start++){unsigned seen=1u<<start,last=0;while(seen!=last){last=seen;for(int u=0;u<n;u++)if(seen>>u&1)seen|=t[u];}if(seen!=(1u<<n)-1)return false;}return true;
}
std::pair<I,I> counts(const T&t){
 int n=t.size(),N=1<<n;std::vector<I>dp(N*n),h(N);h[0]=1;
 for(int u=0;u<n;u++)dp[(1<<u)*n+u]=1;
 for(int mask=1;mask<N;mask++)for(int u=0;u<n;u++){
  I c=dp[mask*n+u];h[mask]+=c;if(!c)continue;
  for(unsigned avail=t[u]&~mask;avail;avail&=avail-1){unsigned bit=avail&-avail;dp[(mask|bit)*n+std::countr_zero(bit)]+=c;}
 }
 I p=0;for(int m=0;m<N;m++)p+=h[m]*h[(N-1)^m];return{h.back(),p};
}
T birth(T t,unsigned cut){int n=t.size();for(int i=0;i<n;i++)if(cut>>i&1)t[i]|=1u<<n;t.push_back(((1u<<n)-1)^cut);return t;}
int main(){
 std::map<unsigned,T> bank{{0,T{0}}};std::cout<<"{\"orders\":[";bool comma=false;
 for(int n=1;n<=7;n++){
  if(comma)std::cout<<',';comma=true;std::cout<<"{\"n\":"<<n<<",\"classes\":"<<bank.size()<<",\"selected\":[";bool hit=false;
  for(auto&[key,t]:bank){if(!strong(t))continue;auto[h,p]=counts(t);if(h!=5&&h!=29&&h!=35)continue;
   if(hit)std::cout<<',';hit=true;std::cout<<"{\"h\":"<<h<<",\"p\":"<<p<<",\"rows\":[";
   for(int i=0;i<n;i++){if(i)std::cout<<',';std::cout<<t[i];}std::cout<<"]}";
  }std::cout<<"]}";
  if(n<7){std::map<unsigned,T> next;for(auto&[key,t]:bank)for(unsigned cut=0;cut<(1u<<n);cut++){auto u=birth(t,cut);next.emplace(canonical(u),u);}bank=std::move(next);}
 }std::cout<<"]}\n";
}
