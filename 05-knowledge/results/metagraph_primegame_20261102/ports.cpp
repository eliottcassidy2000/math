#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using U=uint64_t;
using T=std::vector<unsigned>;
struct Profile { U H,E; std::vector<long long> I; };
Profile profile(const T&r){
 int n=r.size(),M=1<<n; std::vector<U> e(M*n),s(M*n),h(M);h[0]=1;
 for(int mask=1;mask<M;mask++)for(int v=0;v<n;v++)if(mask>>v&1){
  int sub=mask^(1<<v);U a=0,b=0;
  if(!sub)a=b=1;else for(int u=0;u<n;u++)if(sub>>u&1){if(r[u]>>v&1)a+=e[sub*n+u];if(r[v]>>u&1)b+=s[sub*n+u];}
  e[mask*n+v]=a;s[mask*n+v]=b;h[mask]+=a;
 }
 Profile P{h.back(),0,{}};for(int mask=0;mask<M;mask++)P.E+=h[mask]*h[(M-1)^mask];
 for(int x=0;x<n;x++){
  unsigned A=(M-1)^(1<<x),S=r[x]&A;U q=0,g=0;
  for(unsigned Q=A;Q;Q=(Q-1)&A){U a=0;for(int u=0;u<n;u++)if(S>>u&1)a+=s[Q*n+u];q+=h[A^Q]*a;}
  for(int u=0;u<n;u++)if(S>>u&1)g+=e[A*n+u];
  P.I.push_back((long long)(q+g)-(long long)P.H);
 }
 return P;
}
T birth(T r,int x){int n=r.size();for(int u=0;u<n;u++)if(u!=x)r[u]|=1<<n;r.push_back(1<<x);return r;}
T decode(const std::string& c,int n){T r(n);int k=0;for(int a=0;a<n;a++)for(int b=a+1;b<n;b++,k++)if(c[k]=='1')r[a]|=1<<b;else r[b]|=1<<a;return r;}
void rows(const T&r){std::cout<<'[';for(unsigned i=0;i<r.size();i++){if(i)std::cout<<',';std::cout<<r[i];}std::cout<<']';}
bool has(const Profile&p,int sign){return std::find(p.I.begin(),p.I.end(),sign)!=p.I.end();}
bool strong(const T&r){unsigned full=(1u<<r.size())-1;for(int rev=0;rev<2;rev++){unsigned a=1;while(true){unsigned b=a;for(unsigned u=0;u<r.size();u++)if(a>>u&1)for(unsigned v=0;v<r.size();v++)if((rev?r[v]>>u:r[u]>>v)&1)b|=1u<<v;if(a==b)break;a=b;}if(a!=full)return false;}return true;}
int main(int argc,char**argv){
 if(argc<2||argc>3)return 2;std::string path=argv[1];int maxn=argc==3?std::stoi(argv[2]):8;
 std::cout<<"{\"orders\":[";
 for(int n=3;n<=maxn;n++){
  std::ifstream in(path+"/quotients_n"+std::to_string(n)+".txt");if(!in)return 3;
  std::ifstream expected(path+"/expected_n"+std::to_string(n)+".txt");if(!expected)return 6;
  std::string c;int classes=0,both=0,eligible=0,closed=0;std::map<U,int> counts;std::vector<T> representatives,bi,low;
  while(in>>c){if(c.size()!=unsigned(n*(n-1)/2)||c.find_first_not_of("01")!=std::string::npos)return 7;T r=decode(c,n);Profile p=profile(r);U ref;if(!(expected>>ref)||p.H!=ref)return 5;classes++;if(p.H<125&&strong(r))low.push_back(r);if(p.H<5)continue;
   if(has(p,-1)||has(p,1))eligible++;
   if(!has(p,-1)||!has(p,1))continue;
   both++;counts[p.H]++;bi.push_back(r);
   bool ok[2]={false,false};
   for(int x=0;x<n;x++)if(p.I[x]==-1||p.I[x]==1){auto b=profile(birth(r,x));if(b.H!=2*p.H+p.I[x])return 4;if(has(b,-1)&&has(b,1))ok[(p.I[x]+1)/2]=true;}
   if(ok[0]&&ok[1]){closed++;representatives.push_back(r);}
  }
  U extra;if(expected>>extra)return 8;
  if(n!=3)std::cout<<',';
  std::cout<<"{\"n\":"<<n<<",\"classes\":"<<classes<<",\"eligible\":"<<eligible<<",\"both\":"<<both<<",\"both_children_both\":"<<closed<<",\"both_H\":[";
  bool comma=false;for(auto [h,k]:counts){if(comma)std::cout<<',';comma=true;std::cout<<'['<<h<<','<<k<<']';}
  std::cout<<"],\"recursive_seeds\":[";for(unsigned i=0;i<representatives.size();i++){if(i)std::cout<<',';rows(representatives[i]);}
  std::cout<<"],\"both_sources\":[";for(unsigned i=0;i<bi.size();i++){
   if(i)std::cout<<',';auto p=profile(bi[i]);std::cout<<"{\"H\":"<<p.H<<",\"rows\":";rows(bi[i]);std::cout<<",\"children\":[";bool sep=false;
   for(int x=0;x<n;x++)if(p.I[x]==-1||p.I[x]==1){if(sep)std::cout<<',';sep=true;auto b=profile(birth(bi[i],x));std::cout<<"{\"mark\":"<<x<<",\"sign\":"<<p.I[x]<<",\"I\":[";for(unsigned y=0;y<b.I.size();y++){if(y)std::cout<<',';std::cout<<b.I[y];}std::cout<<"]}";}
   std::cout<<"]}";
  }std::cout<<"],\"strong_below125\":[";for(unsigned i=0;i<low.size();i++){
   if(i)std::cout<<',';auto p=profile(low[i]);std::cout<<"{\"H\":"<<p.H<<",\"rows\":";rows(low[i]);std::cout<<",\"I\":[";for(unsigned y=0;y<p.I.size();y++){if(y)std::cout<<',';std::cout<<p.I[y];}std::cout<<"]}";
  }std::cout<<"]}";std::cout.flush();
 }
 std::cout<<"]}\n";
}
