#include <algorithm>
#include <array>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>
using namespace std;
using Poly=array<long long,10>;
vector<int> decode(string code,int n){vector<int> r(n);int b=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++)if(code[b]=='1')r[i]|=1<<j;else r[j]|=1<<i;return r;}
bool strong(vector<int> r){int n=r.size();for(int pass=0;pass<2;pass++){int s=1,old;do{old=s;for(int i=0;i<n;i++)if(s>>i&1)s|=r[i];}while(s!=old);if(s!=(1<<n)-1)return false;vector<int>t(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(r[i]>>j&1)t[j]|=1<<i;r=t;}return true;}
Poly skeleton(const vector<int>& rows,int v){
 int n=rows.size(),m=n-1;vector<int>e;for(int u=0;u<n;u++)if(u!=v)e.push_back(u);e.push_back(v);
 vector<vector<Poly>> memo(1<<m,vector<Poly>(n));vector<vector<bool>>done(1<<m,vector<bool>(n));
 function<Poly(int,int)> f=[&](int S,int last){
  if(done[S][last])return memo[S][last];Poly ans={};if(!S)ans[0]=1;
  for(int j=0;j<m;j++)if((S>>j&1)&&(rows[e[last]]>>e[j]&1)){auto z=f(S^(1<<j),j);for(int k=0;k<=n;k++)ans[k]+=z[k];}
  if(last!=m&&(rows[e[last]]>>v&1)){auto z=f(S,m);for(int k=0;k<n;k++)ans[k+1]+=z[k];}
  done[S][last]=true;return memo[S][last]=ans;
 };
 int full=(1<<m)-1;Poly out={};for(int j=0;j<m;j++){auto z=f(full^(1<<j),j);for(int k=0;k<=n;k++)out[k]+=z[k];}
 auto z=f(full,m);for(int k=0;k<n;k++)out[k+1]+=z[k];return out;
}
long long binom(int n,int k){if(k<0||k>n)return 0;long long x=1;for(int i=1;i<=k;i++)x=x*(n-i+1)/i;return x;}
Poly covers(const vector<int>&r){int n=r.size();vector<int>w(n);iota(w.begin(),w.end(),0);Poly c={};do{int bad=0;for(int i=0;i<n-1;i++)if(!(r[w[i]]>>w[i+1]&1))bad++;for(int k=1;k<=n;k++)c[k]+=binom(n-1-bad,k-1-bad);}while(next_permutation(w.begin(),w.end()));return c;}
bool prime(long long n){if(n<2)return false;for(long long d=2;d*d<=n;d++)if(n%d==0)return false;return true;}
int main(int argc,char**argv){
 vector<pair<int,Poly>>banks;for(int p:{3,5,7}){set<Poly>seen;for(int b=0;b<(1<<((p-1)/2));b++){vector<int>r(p);for(int u=0;u<p;u++)for(int d=1;d<=p/2;d++){int step=b>>(d-1)&1?p-d:d;r[u]|=1<<((u+step)%p);}Poly c=covers(r);if(seen.insert(c).second)banks.push_back({p,c});}}
 cout<<"{\"cover_vectors\":[";for(size_t i=0;i<banks.size();i++){if(i)cout<<',';cout<<'[';for(int k=0;k<=banks[i].first;k++){if(k)cout<<',';cout<<banks[i].second[k];}cout<<']';}cout<<"],\"orders\":[";
 set<long long>hits;bool comma=false;
 for(int n=3;n<=8;n++){
  if(argc!=2){cerr<<"usage: module_audit INPUT_DIRECTORY\n";return 2;}
  string name=string(argv[1])+"/quotients_n"+to_string(n)+".txt";FILE*p=fopen(name.c_str(),"r");if(!p)return 2;
  char buf[1024];long long all=0,sc=0;vector<long long>mins(banks.size(),1LL<<60);
  while(fgets(buf,sizeof(buf),p)){string code(buf);auto r=decode(code,n);all++;if(!strong(r))continue;sc++;
   for(int v=0;v<n;v++){Poly b=skeleton(r,v);for(size_t i=0;i<banks.size();i++){long long h=0;for(int k=1;k<=n;k++)h+=b[k]*banks[i].second[k];mins[i]=min(mins[i],h);if(banks[i].first==7&&h<30625&&h%7==0&&prime(h/7))hits.insert(h);}}
  }
  if(fclose(p)!=0)return 3;
  if(comma)cout<<',';comma=true;cout<<"{\"n\":"<<n<<",\"all_classes\":"<<all<<",\"strong_classes\":"<<sc<<",\"minima\":[";for(size_t i=0;i<mins.size();i++){if(i)cout<<',';cout<<mins[i];}cout<<"]}";
 }
 cout<<"],\"symmetric_7p_below_30625\":[";comma=false;for(auto h:hits){if(comma)cout<<',';comma=true;cout<<h;}cout<<"]}\n";
}
