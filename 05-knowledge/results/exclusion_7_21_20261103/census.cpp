// Independent complete Hamilton-frame census, with literal path recursion.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using Graph=vector<unsigned>;
void need(bool x,const char* why){if(!x)throw runtime_error(why);}
unsigned code(const Graph& a,const vector<int>& p){
 unsigned c=0;int b=0,n=a.size();
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++)if((a[p[i]]>>p[j])&1)c|=1u<<b;
 return c;
}
Graph decode(int n,unsigned c){Graph a(n);int b=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++)
 if((c>>b)&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
int hp(const Graph& a){
 int n=a.size(),N=1<<n; if(!n)return 1;
 vector<array<int,7>> d(N);for(int i=0;i<n;i++)d[1<<i][i]=1;
 for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(d[s][i]){
  unsigned options=a[i]&(~s)&(N-1);while(options){int j=__builtin_ctz(options);options&=options-1;d[s|(1<<j)][j]+=d[s][i];}}
 return accumulate(d.back().begin(),d.back().end(),0);
}
int walk(const Graph& a,unsigned used,int v){
 if(used==((1u<<a.size())-1))return 1;
 int total=0;for(unsigned q=a[v]&~used;q;q&=q-1){int j=__builtin_ctz(q);total+=walk(a,used|(1u<<j),j);}return total;
}
bool strong(const Graph& a){
 int n=a.size();if(!n)return false;Graph r=a;for(int i=0;i<n;i++)r[i]|=1u<<i;
 for(int k=0;k<n;k++)for(int i=0;i<n;i++)if((r[i]>>k)&1)r[i]|=r[k];
 return all_of(r.begin(),r.end(),[n](unsigned x){return x==(1u<<n)-1;});
}
bool no_cut(const Graph& a){
 unsigned N=1u<<a.size();for(unsigned s=1;s<N-1;s++){
  unsigned outside=(N-1)^s;bool dominates=true;
  for(unsigned q=s;q;q&=q-1)if((a[__builtin_ctz(q)]&outside)!=outside){dominates=false;break;}
  if(dominates)return false;
 }return true;
}
Graph induced(const Graph& a,unsigned mask){
 vector<int> v;for(int i=0;i<(int)a.size();i++)if((mask>>i)&1)v.push_back(i);
 Graph b(v.size());for(int i=0;i<(int)v.size();i++)for(int j=0;j<(int)v.size();j++)if((a[v[i]]>>v[j])&1)b[i]|=1u<<j;return b;
}
unsigned canonical(const Graph& a){
 int n=a.size();vector<int> p(n);iota(p.begin(),p.end(),0);vector<int>d(n);for(int i=0;i<n;i++)d[i]=__builtin_popcount(a[i]);
 unsigned best=~0u;
 do{bool sorted=true;for(int i=1;i<n;i++)if(d[p[i-1]]>d[p[i]]){sorted=false;break;}
 if(sorted)best=min(best,code(a,p));}while(next_permutation(p.begin(),p.end()));
 return best;
}
vector<vector<int>> autos(const Graph&a){
 int n=a.size();vector<int>p(n);iota(p.begin(),p.end(),0);vector<vector<int>> out;
 do{bool ok=true;for(int i=0;i<n&&ok;i++){if(__builtin_popcount(a[i])!=__builtin_popcount(a[p[i]])){ok=false;break;}
 for(int j=i+1;j<n;j++)if(((a[i]>>j)&1)!=((a[p[i]]>>p[j])&1)){ok=false;break;}}
 if(ok)out.push_back(p);}while(next_permutation(p.begin(),p.end()));return out;
}
int orbit_count(const vector<vector<int>>& group,int n,int arity){
 set<int>seen;int count=0;for(int i=0;i<n;i++)for(int j=0;j<(arity==1?1:n);j++){
 if(arity==2&&i==j)continue;int x=arity==1?i:i*n+j;if(seen.count(x))continue;count++;
 for(const auto&p:group)seen.insert(arity==1?p[i]:p[i]*n+p[j]);}return count;
}
int main(){
 cout<<"{\"status\":\"PASS\",\"universes\":[";
 for(int n=1;n<=7;n++){
  vector<pair<int,int>> free;for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)free.emplace_back(i,j);
  unsigned N=1u<<free.size();set<unsigned>classes;map<int,int> hist,stronghist;int floor=100000;
  for(unsigned mask=0;mask<N;mask++){
   Graph a(n);for(int i=0;i<n-1;i++)a[i]|=1u<<(i+1);
   for(int b=0;b<(int)free.size();b++){auto [i,j]=free[b];if((mask>>b)&1)a[i]|=1u<<j;else a[j]|=1u<<i;}
   int h=hp(a),ind=0;for(int i=0;i<n;i++)ind+=walk(a,1u<<i,i);need(h==ind,"path algorithms disagree");need(h%2==1,"oddness");
   bool s=strong(a);need(s==no_cut(a),"strong tests disagree");hist[h]++;if(s){stronghist[h]++;floor=min(floor,h);}
   need(h!=7&&h!=21,"forbidden finite value");classes.insert(canonical(a));
  }
  if(n>1)cout<<',';cout<<"{\"n\":"<<n<<",\"frames\":"<<N<<",\"class_count\":"<<classes.size()<<",\"strong_floor\":"<<(floor==100000?0:floor)<<",\"histogram\":{";
  bool first=true;for(auto [h,c]:hist){if(!first)cout<<',';first=false;cout<<'\"'<<h<<"\":"<<c;}cout<<"},\"strong_histogram\":{";
  first=true;for(auto [h,c]:stronghist){if(!first)cout<<',';first=false;cout<<'\"'<<h<<"\":"<<c;}cout<<"},\"classes\":[";
  first=true;long long labelled=0;
  for(unsigned c:classes){Graph a=decode(n,c);int h=hp(a);bool s=strong(a);auto group=autos(a);long long fac=1;for(int i=2;i<=n;i++)fac*=i;labelled+=fac/group.size();
   int maxdel=0;bool removable=n<4;for(int i=0;i<n;i++){auto b=induced(a,((1u<<n)-1)^(1u<<i));maxdel=max(maxdel,hp(b));removable|=strong(b);}
   if(s&&n>=4)need(removable,"strong deletion");
   if(!first)cout<<',';first=false;cout<<"{\"code\":"<<c<<",\"H\":"<<h<<",\"strong\":"<<(s?"true":"false")<<",\"aut\":"<<group.size()<<",\"roots\":"<<orbit_count(group,n,1)<<",\"ordered_distinct_ports\":"<<orbit_count(group,n,2)<<",\"max_delete_H\":"<<maxdel<<",\"minimal_above_21\":"<<((h>21&&maxdel<=21)?"true":"false")<<'}';
  }
  need(labelled==(1ll<<(n*(n-1)/2)),"labelled orbit mass");cout<<"],\"labelled_mass\":"<<labelled<<'}';
 }
 cout<<"]}\n";
}
