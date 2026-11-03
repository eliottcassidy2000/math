// Exact rooted neighbourhood census. Root 0 beats a vertices and loses to b.
// Complete side-orbit banks + every a*b cross bit cover every rooted class.
#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using namespace std;using T=vector<unsigned>;using I=int64_t;
I hp(const T&t){int n=t.size(),f=(1<<n)-1;vector<I>d((f+1)*n);
 for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;
 for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(I c=d[s*n+i]){
  unsigned bits=t[i]&~s;while(bits){unsigned bit=bits&-bits;bits-=bit;d[(s|bit)*n+countr_zero(bit)]+=c;}}
 I h=0;for(int i=0;i<n;i++)h+=d[f*n+i];return h;}
T clone(const T&t,int m){int n=t.size();T r(n+m-1);for(int i=0;i<(int)r.size();i++)for(int j=0;j<(int)r.size();j++){
 int a=i<m?0:i-m+1,b=j<m?0:j-m+1;
 if((a==b&&i<j)||(a!=b&&(t[a]>>b&1)))r[i]|=1u<<j;}return r;}
bool strong(const T&t){int n=t.size();unsigned f=(1u<<n)-1;for(int side=0;side<2;side++){
 unsigned seen=1,old=0;while(seen!=old){old=seen;for(int i=0;i<n;i++)if(seen>>i&1)seen|=side?(f^t[i]^(1u<<i)):t[i];}if(seen!=f)return false;}return true;}
vector<T> read(const string&p,int n){ifstream f(p);int k;f>>k;if(!f)throw runtime_error("bank missing");vector<T>z(k,T(n));for(auto&t:z)for(auto&x:t)f>>x;return z;}
void print(const T&t){cout<<'[';for(int i=0;i<(int)t.size();i++){if(i)cout<<',';cout<<t[i];}cout<<']';}
int main(int argc,char**argv){if(argc<5)return 2;string dir=argv[1];int a=stoi(argv[2]),b=stoi(argv[3]),level=stoi(argv[4]);
 auto aa=read(dir+"/bank"+to_string(a)+".txt",a),bb=read(dir+"/bank"+to_string(b)+".txt",b);
 I multiplier=1;for(int i=0;i<level;i++)multiplier*=3;vector<int> targets=level==1?vector<int>{95,97}:vector<int>{287,289,293,295};
 I A=8*multiplier,B=6*multiplier;vector<I> hits(targets.size()),scalar(targets.size());I examined=0;
 vector<unsigned> crosses;
 if(argc>5 && string(argv[5])=="cover2"){
  set<unsigned> masks;
  for(int u=0;u<a+b;u++)for(int v=u+1;v<a+b;v++){
   unsigned allowed=0;
   for(int i=0;i<a;i++)for(int j=0;j<b;j++)if(i==u||i==v||a+j==u||a+j==v)allowed|=1u<<(i*b+j);
   unsigned s=allowed;while(true){masks.insert(s);if(!s)break;s=(s-1)&allowed;}
  }
  crosses.assign(masks.begin(),masks.end());
 }else{for(unsigned c=0;c<(1u<<(a*b));c++)crosses.push_back(c);}
 for(auto&x:aa)for(auto&y:bb)for(unsigned cross:crosses){
  int n=1+a+b;T t(n);for(int i=0;i<a;i++){t[0]|=1u<<(1+i);t[1+i]|=x[i]<<1;}
  for(int j=0;j<b;j++){t[1+a+j]|=1;t[1+a+j]|=y[j]<<(1+a);}
  for(int i=0;i<a;i++)for(int j=0;j<b;j++){if(cross>>(i*b+j)&1)t[1+i]|=1u<<(1+a+j);else t[1+a+j]|=1u<<(1+i);}
  examined++;I h=hp(t);auto it=find(targets.begin(),targets.end(),h);if(it==targets.end())continue;int idx=it-targets.begin();scalar[idx]++;
  if(!strong(t))continue;I C=h-3*A-2*B;bool ok=true;int maxruns=min(a,b)+1;
  for(int m=2;m<=maxruns;m++){I p3=1,p2=1;for(int j=0;j<m;j++){p3*=3;p2*=2;}if(hp(clone(t,m))!=A*p3+B*p2+C){ok=false;break;}}
  if(ok){hits[idx]++;cout<<"{\"H\":"<<h<<",\"root\":0,\"a\":"<<a<<",\"b\":"<<b<<",\"level\":"<<level<<",\"rows\":";print(t);cout<<"}\n";}
 }
 cerr<<"examined "<<examined<<" scalar ";for(auto k:scalar)cerr<<k<<' ';cerr<<" uniform ";for(auto k:hits)cerr<<k<<' ';cerr<<'\n';
}
