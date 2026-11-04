#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <stdexcept>
using namespace std;
using U=uint64_t;
void require(bool x,const char*s){if(!x)throw runtime_error(s);}
struct Profiles {U h;vector<U> f,s;};
struct Counter {
 vector<int> r;int n,full;vector<U> h,e,s;U count=0,cap;
 Counter(const vector<int>& a):r(a),n(a.size()),full((1<<n)-1){}
 // Enumerate actual full directed words, globally stopping only above the required cutoff.
 void fullpaths(int used,int last){if(count>=cap)return;if(used==full){count++;return;}int next=r[last]&~used;while(next){int bit=next&-next;next-=bit;fullpaths(used|bit,__builtin_ctz((unsigned)bit));if(count>=cap)return;}}
 U bounded(U cutoff=202){count=0;cap=cutoff;for(int i=0;i<n&&count<cap;i++)fullpaths(1<<i,i);return count;}
 // Enumerate every actual simple directed path, not subset dynamic programming.
 void partial(int used,int first,int last){h[used]++;e[used*n+last]++;s[used*n+first]++;int next=r[last]&~used;while(next){int bit=next&-next;next-=bit;partial(used|bit,first,__builtin_ctz((unsigned)bit));}}
 Profiles profile(){h.assign(full+1,0);h[0]=1;e.assign((full+1)*n,0);s.assign((full+1)*n,0);for(int i=0;i<n;i++)partial(1<<i,i,i);Profiles p{h[full],vector<U>(n),vector<U>(n)};for(int v=0;v<n;v++){p.s[v]=s[full*n+v];for(int m=1;m<=full;m++)p.f[v]+=e[m*n+v]*h[full^m];}return p;}
};
bool strong(const vector<int>&r){int n=r.size();vector<vector<bool>> reach(n,vector<bool>(n));for(int i=0;i<n;i++)for(int j=0;j<n;j++)reach[i][j]=(i==j)||((r[i]>>j)&1);for(int k=0;k<n;k++)for(int i=0;i<n;i++)for(int j=0;j<n;j++)reach[i][j]=reach[i][j]||(reach[i][k]&&reach[k][j]);for(auto&row:reach)for(bool v:row)if(!v)return false;return true;}
void valid(const vector<int>&r){int n=r.size();for(int i=0;i<n;i++){require(r[i]>=0&&r[i]<(1<<n)&&!(r[i]>>i&1),"row bounds");for(int j=0;j<i;j++)require(((r[i]>>j)&1)+((r[j]>>i)&1)==1,"orientation");}}
vector<int> erase_vertex(const vector<int>&r,int v){vector<int> keep;for(int i=0;i<(int)r.size();i++)if(i!=v)keep.push_back(i);vector<int> a(keep.size());for(int i=0;i<(int)keep.size();i++)for(int j=0;j<(int)keep.size();j++)if(r[keep[i]]>>keep[j]&1)a[i]|=1<<j;return a;}
int main(int argc,char**argv){require(argc>=2,"mode required");string mode=argv[1];
 if(mode=="small"){set<pair<U,U>>profiles;for(int n=3;n<=6;n++){U total=0,st=0,unit=0;for(int mask=0;mask<(1<<(n*(n-1)/2));mask++){vector<int>r(n);int b=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++)if(mask>>b&1)r[i]|=1<<j;else r[j]|=1<<i;total++;if(!strong(r))continue;st++;auto p=Counter(r).profile();for(int v=0;v<n;v++)if(p.s[v]==1){unit++;profiles.insert({p.h,p.f[v]});}}cout<<"SMALL "<<n<<" "<<total<<" "<<st<<" "<<unit<<"\n";}U matches=0;for(auto [h,f]:profiles)matches+=profiles.count({h+2,f});require(!matches,"small twin");cout<<"PROFILES "<<profiles.size()<<" MATCHES "<<matches<<"\n";return 0;}
 int n,K;cin>>n>>K;require(n==8||n==9,"input order");U tested=0,st=0,low=0,found=0,max7=0;U tc[2]={0,0};set<U>fs[2];set<vector<int>> lows;
 auto inspect=[&](const vector<int>&r){valid(r);tested++;Counter counter(r);U h=counter.bounded();if(r.size()==7)max7=max(max7,h);if(!strong(r))return;st++;if(r.size()==9&&h<=201){low++;lows.insert(r);}if(h==191||h==201){int target=h==191?0:1;tc[target]++;auto p=counter.profile();require(p.h==h,"bounded Hamilton count");for(auto f:p.f){fs[target].insert(f);if(f==(target==0?1273:5303))found++;}}};
 for(int index=0;index<K;index++){vector<int>r(n);for(int&i:r)cin>>i;valid(r);if(mode=="bank"){inspect(r);for(int v=0;v<n;v++)inspect(erase_vertex(r,v));}else {require(mode=="extend","mode");for(int mask=0;mask<(1<<n);mask++){vector<int>a=r;a.push_back(mask);for(int i=0;i<n;i++)if(!(mask>>i&1))a[i]|=1<<n;inspect(a);}}}
 require(found==0,"target found");cout<<"SUMMARY "<<mode<<" "<<n<<" "<<K<<" "<<tested<<" "<<st<<" "<<low<<" "<<found<<" "<<max7<<"\n";for(int t=0;t<2;t++){cout<<"TARGET "<<(t?201:191)<<" "<<tc[t]<<" "<<fs[t].size();for(U f:fs[t])cout<<" "<<f;cout<<"\n";}for(auto&r:lows){cout<<"LOW";for(int x:r)cout<<" "<<x;cout<<"\n";}
}
