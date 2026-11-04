#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <array>
#include <limits>
using namespace std; using I=int64_t; using T=vector<unsigned>;
struct Profile {I h,f;vector<I> all_f;};
Profile profile(const T&t,int w){int n=t.size(),all=(1<<n)-1; vector<I>dp((all+1)*n),h(all+1); h[0]=1;
 for(int v=0;v<n;v++)dp[(1<<v)*n+v]=1;
 for(int s=1;s<=all;s++)for(unsigned q=s;q;q&=q-1){int v=countr_zero(q); I k=dp[s*n+v]; h[s]+=k; if(k)for(unsigned z=t[v]&~s;z;z&=z-1)dp[(s|(z&-z))*n+countr_zero(z)]+=k;}
 vector<I>fs(n);for(int s=1;s<=all;s++)for(unsigned q=s;q;q&=q-1){int v=countr_zero(q);fs[v]+=dp[s*n+v]*h[all^s];}return {h[all],fs[w],fs};}
bool strong(const T&t){unsigned all=(1<<t.size())-1;for(int rev=0;rev<2;rev++){unsigned reach=1,last=0;while(reach!=last){last=reach;for(int v=0;v<(int)t.size();v++)if(reach>>v&1)reach|=rev?(all^t[v]^(1<<v)):t[v];}if(reach!=all)return false;}return true;}
void rows(const T&t){cout<<'[';for(int v=0;v<(int)t.size();v++){if(v)cout<<',';cout<<t[v];}cout<<']';}
int main(){int n,w,depth,moving,min_prefix,max_prefix;I source_H;cin>>n>>w>>depth>>moving>>source_H>>min_prefix>>max_prefix; if(n<1||n>16||w<0||w>=n||depth<0||depth>3||min_prefix<1||max_prefix>60||min_prefix>max_prefix)return 2;T source(n);for(auto&v:source)cin>>v;auto base=profile(source,w);vector<pair<int,int>>arcs;for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)arcs.push_back({a,b});
 array<I,2>targets={3*source_H+2,3*source_H+4};array<I,2>near={numeric_limits<I>::max(),numeric_limits<I>::max()};int checked=0,strong_count=0,hits=0;cout<<"{\"source_H\":"<<source_H<<",\"source_h\":"<<base.h<<",\"source_F\":"<<base.f<<",\"maximum_flips\":"<<depth<<",\"minimum_prefix\":"<<min_prefix<<",\"maximum_prefix\":"<<max_prefix<<",\"hits\":[";
 auto test=[&](vector<int>edits){T t=source;for(int j:edits){auto[a,b]=arcs[j];t[a]^=1<<b;t[b]^=1<<a;}checked++;if(!strong(t))return;strong_count++;auto p=profile(t,w);
  for(int mark=0;mark<n;mark++){if(!moving&&mark!=w)continue;I f=p.all_f[mark];for(int s=min_prefix;s<=max_prefix;s++){__int128 wide=p.h+(__int128(1)<<s)*f;if(wide>numeric_limits<I>::max())break;I value=I(wide);for(int branch=0;branch<2;branch++){I delta=value-targets[branch];near[branch]=min(near[branch],delta<0?-delta:delta);if(value==targets[branch]){if(hits++)cout<<',';cout<<"{\"offset\":"<<2+2*branch<<",\"prefix_size\":"<<s<<",\"mark\":"<<mark<<",\"h\":"<<p.h<<",\"F\":"<<f<<",\"uniform_same_clock\":"<<((p.h==3*base.h+2+2*branch&&f==3*base.f)?"true":"false")<<",\"edges\":[";for(int k=0;k<(int)edits.size();k++){if(k)cout<<',';auto[a,b]=arcs[edits[k]];cout<<'['<<a<<','<<b<<']';}cout<<"],\"rows\":";rows(t);cout<<'}'<<flush;}}if(value>targets[1])break;}}
 };
 test({});if(depth>=1)for(int i=0;i<(int)arcs.size();i++)test({i});if(depth>=2)for(int i=0;i<(int)arcs.size();i++)for(int j=0;j<i;j++)test({j,i});if(depth>=3)for(int i=0;i<(int)arcs.size();i++)for(int j=0;j<i;j++)for(int k=0;k<j;k++)test({k,j,i});
 cout<<"],\"checked\":"<<checked<<",\"strong\":"<<strong_count<<",\"nearest\":["<<near[0]<<','<<near[1]<<"]}\n";
}
