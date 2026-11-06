#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
using U=uint64_t;
vector<unsigned> decode(U code,int n){vector<unsigned>a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(code>>k&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
U encode(const vector<unsigned>&a){U r=0;int k=0,n=a.size();for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(a[i]>>j&1)r|=U(1)<<k;return r;}
U canon(const vector<unsigned>&a){int n=a.size();vector<int>p(n),d(n);iota(p.begin(),p.end(),0);for(int i=0;i<n;i++)d[i]=__builtin_popcount(a[i]);sort(p.begin(),p.end(),[&](int i,int j){return d[i]<d[j]||(d[i]==d[j]&&i<j);});vector<pair<int,int>>g;for(int i=0,j;i<n;i=j){j=i+1;while(j<n&&d[p[j]]==d[p[i]])j++;g.push_back({i,j});}U best=~U(0);auto rec=[&](auto&&self,int t)->void{if(t==(int)g.size()){U c=0;int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(a[p[i]]>>p[j]&1)c|=U(1)<<k;best=min(best,c);return;}auto [l,r]=g[t];sort(p.begin()+l,p.begin()+r);do{self(self,t+1);}while(next_permutation(p.begin()+l,p.begin()+r));};rec(rec,0);return best;}
int hp(const vector<unsigned>&a){int n=a.size();int dp[256][8]={};for(int i=0;i<n;i++)dp[1<<i][i]=1;for(int s=1;s<(1<<n);s++)for(int i=0;i<n;i++){unsigned z=a[i]&~s;while(z){int j=__builtin_ctz(z);z&=z-1;dp[s|(1<<j)][j]+=dp[s][i];}}return accumulate(dp[(1<<n)-1],dp[(1<<n)-1]+n,0);}
bool strong(const vector<unsigned>&a){int n=a.size();unsigned full=(1u<<n)-1;for(int v=0;v<n;v++){unsigned s=1u<<v,old;do{old=s;for(int i=0;i<n;i++)if(s>>i&1)s|=a[i];}while(s!=old);if(s!=full)return false;}return true;}
int brute(const vector<unsigned>&a){int n=a.size(),c=0;vector<int>p(n);iota(p.begin(),p.end(),0);do{bool ok=1;for(int i=1;i<n;i++)ok&=bool(a[p[i-1]]>>p[i]&1);c+=ok;}while(next_permutation(p.begin(),p.end()));return c;}
void arrayprint(const set<int>&s){cout<<"[";bool first=1;for(int x:s){if(!first)cout<<",";cout<<x;first=0;}cout<<"]";}
int main(){vector<U>reps{0};int expected[]={0,1,1,2,4,12,56,456};cout<<"{\"convention\":\"bit k for lexicographic pair i<j is 1 iff i beats j\",\"orders\":[\n";for(int n=1;n<=8;n++){set<int>all,sc;map<int,U>witness;set<U>next;U candidates=0;auto visit=[&](const vector<unsigned>&a){int h=hp(a);assert(h%2==1);if(n<=5)assert(h==brute(a));all.insert(h);if(strong(a))sc.insert(h);witness.emplace(h,encode(a));candidates++;if(n<=7)next.insert(canon(a));};if(n==1)visit(vector<unsigned>{0});else for(U r:reps){auto a=decode(r,n-1);for(unsigned pat=0;pat<(1u<<(n-1));pat++){auto b=a;b.push_back(pat);for(int i=0;i<n-1;i++)if(!(pat>>i&1))b[i]|=1u<<(n-1);visit(b);}}if(n<=7){assert(next.size()==(unsigned)expected[n]);reps.assign(next.begin(),next.end());}set<int>gaps;for(int i=1;i<=*all.rbegin();i+=2)if(!all.count(i))gaps.insert(i);if(n>1)cout<<",\n";cout<<"{\"n\":"<<n<<",\"augmentation_candidates\":"<<candidates<<",\"classes\":"<<(n<=7?(int)next.size():-1)<<",\"spectrum\":";arrayprint(all);cout<<",\"strong_spectrum\":";arrayprint(sc);cout<<",\"odd_gaps\":";arrayprint(gaps);cout<<",\"witness_codes\":{";bool first=1;for(auto [h,c]:witness){if(!first)cout<<",";cout<<"\""<<h<<"\":"<<c;first=0;}cout<<"}}";cerr<<"n="<<n<<" count="<<all.size()<<" max="<<*all.rbegin()<<" candidates="<<candidates<<"\n";}cout<<"\n]}\n";}
