// Independent path-word/permutation checker; no subset DP or clone expansion.
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using Matrix=vector<vector<int>>;
vector<vector<int>> bank(string path,int n){ifstream f(path);int k;f>>k;if(!f)throw 1;
 vector<vector<int>>b(k,vector<int>(n));for(auto&t:b)for(auto&x:t)f>>x;return b;}
int count_paths(const Matrix&t,int limit){int n=t.size(),count=0;
 auto go=[&](auto&&self,int last,int mask)->void{
  if(count>limit)return;if(mask==(1<<n)-1){count++;return;}
  for(int v=0;v<n;v++)if(!(mask>>v&1)&&(last<0||t[last][v]))self(self,v,mask|(1<<v));
 };go(go,-1,0);return count;}
bool connected(const Matrix&t){int n=t.size();for(int dir=0;dir<2;dir++){
 vector<int>seen(n);seen[0]=1;for(int pass=0;pass<n;pass++)for(int u=0;u<n;u++)if(seen[u])for(int v=0;v<n;v++)if(dir?t[v][u]:t[u][v])seen[v]=1;
 if(accumulate(seen.begin(),seen.end(),0)!=n)return false;}return true;}
long long choose(int n,int k){if(k<0||k>n)return 0;long long a=1;for(int j=1;j<=k;j++)a=a*(n-j+1)/j;return a;}
vector<long long> skeleton(const Matrix&t){int n=t.size();vector<int>p(n-1);iota(p.begin(),p.end(),1);vector<long long>c(n+1);
 do{int mandatory=0,optional=t[0][p.front()]+t[p.back()][0];bool ok=true;
  for(int j=0;j<n-2;j++){int a=p[j],b=p[j+1];bool no=t[a][b],yes=t[a][0]&&t[0][b];
   if(!no&&!yes){ok=false;break;}if(!no)mandatory++;else if(yes)optional++;}
  if(ok)for(int r=0;r<=n;r++)c[r]+=choose(optional,r-mandatory);
 }while(next_permutation(p.begin(),p.end()));return c;}
int main(int argc,char**argv){if(argc!=5)return 2;string dir=argv[1];int a=stoi(argv[2]),b=stoi(argv[3]),level=stoi(argv[4]),n=1+a+b;
 auto out=bank(dir+"/bank"+to_string(a)+".txt",a),in=bank(dir+"/bank"+to_string(b)+".txt",b);
 vector<int>targets=level==1?vector<int>{95,97}:vector<int>{287,289,293,295};vector<long long>scalar(targets.size()),hits(targets.size());
 long long cases=0;int factor=level==1?3:9;
 for(auto&o:out)for(auto&i:in)for(int mask=0;mask<(1<<(a*b));mask++){
  Matrix t(n,vector<int>(n));for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
   int bit;if(u==0)bit=v<=a;else if(v<=a)bit=o[u-1]>>(v-1)&1;
   else if(u>a)bit=i[u-a-1]>>(v-a-1)&1;else bit=mask>>((u-1)*b+v-a-1)&1;
   t[u][v]=bit;t[v][u]=1-bit;
  }
  cases++;int h=count_paths(t,targets.back());auto found=find(targets.begin(),targets.end(),h);if(found==targets.end())continue;
  int index=found-targets.begin();scalar[index]++;if(!connected(t))continue;
  auto c=skeleton(t);if(c[1]!=h)throw 2;bool ok=c[2]==30*factor&&c[3]==8*factor;
  for(int r=4;r<=n;r++)if(c[r])ok=false;
  if(ok)hits[index]++;
 }
 cout<<"examined "<<cases<<" scalar ";for(auto x:scalar)cout<<x<<' ';cout<<" uniform ";for(auto x:hits)cout<<x<<' ';cout<<'\n';
}
