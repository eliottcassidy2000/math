#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using I=long long;
struct Hit {vector<int> rows;int mask;I F,q,g,h;};
map<I,Hit> hits;
uint64_t parents=0,cuts=0;
void test(vector<int> rows){
 int n=rows.size(),full=(1<<n)-1;parents++;
 vector<vector<I>> end(1<<n,vector<I>(n)),start=end;
 vector<I> H(1<<n);H[0]=1;
 for(int v=0;v<n;v++)end[1<<v][v]=start[1<<v][v]=1;
 for(int m=1;m<=full;m++){
  for(int v=0;v<n;v++)if(m>>v&1){
   int a=m^(1<<v);
   if(a)for(int u=0;u<n;u++)if(a>>u&1){
    if(rows[u]>>v&1)end[m][v]+=end[a][u];
    if(rows[v]>>u&1)start[m][v]+=start[a][u];
   }
   H[m]+=end[m][v];
  }
 }
 vector<I> qv(n),gv=end[full];vector<vector<I>> W(n,vector<I>(n)),b=W;
 for(int v=0;v<n;v++)for(int m=1;m<=full;m++)if(m>>v&1)qv[v]+=start[m][v]*H[full^m];
 for(int m=1;m<full;m++)for(int u=0;u<n;u++)if(m>>u&1)for(int v=0;v<n;v++)if(!(m>>v&1))W[u][v]+=end[m][u]*start[full^m][v];
 for(int u=0;u<n;u++)for(int v=0;v<n;v++){
  if((W[u][v]+W[v][u])%2)throw 1;
  b[u][v]=(W[u][v]+W[v][u])/2;
 }
 for(int mask=1;mask<=full;mask++){
  cuts++;I q=0,g=0,L=0,C=0;
  for(int v=0;v<n;v++)if(mask>>v&1){q+=qv[v];g+=gv[v];for(int u=0;u<n;u++){
    if(mask>>u&1){if(u<v)L+=b[u][v];}
    else if(rows[v]>>u&1)C+=b[u][v];
  }}
  if(q!=2*C+g+2*L)throw 2;
  if(q && !(q&(q-1))){
   I c=2*L-H[full];
   if(c>0 && !hits.count(c))hits[c]={rows,mask,H[full]+2*C,q,g,H[full]};
  }
 }
}
uint64_t rnd=123456789;
uint64_t next(){rnd^=rnd<<13;rnd^=rnd>>7;rnd^=rnd<<17;return rnd;}
int main(int argc,char**argv){
 int maxn=argc>1?stoi(argv[1]):6,samples=argc>2?stoi(argv[2]):0;
 bool fixed=argc>3 && string(argv[3])=="fixed";
 for(int n=1;n<=maxn;n++){
  int e=n*(n-1)/2-(fixed?n-1:0);
  for(uint64_t a=0;a<(1ULL<<e);a++){
   vector<int> rows(n);int j=0;
   for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
    if(fixed && v==u+1){rows[u]|=1<<v;continue;}
    if(a>>j&1)rows[u]|=1<<v;else rows[v]|=1<<u;j++;
   }
   test(rows);
  }
  cerr<<"complete n="<<n<<" deficits="<<hits.size()<<" parents="<<parents<<"\n";
 }
 for(int n=maxn+1;n<=10 && samples;n++)for(int s=0;s<samples;s++){
  vector<int> rows(n);for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)if(next()&1)rows[u]|=1<<v;else rows[v]|=1<<u;
  test(rows);
 }
 cout<<"{\"parent_gauge\":\""<<(fixed?"fixed_hamilton_path":"all_labelled")<<"\",\"complete_parent_order\":"<<maxn<<",\"samples_each_higher_order_to10\":"<<samples<<",\"parents\":"<<parents<<",\"cuts\":"<<cuts<<",\"hits\":[";bool first=true;
 for(auto [c,h]:hits){if(!first)cout<<",";first=false;cout<<"{\"deficit\":"<<c<<",\"mask\":"<<h.mask<<",\"state\":["<<h.F<<","<<h.q<<","<<h.g<<","<<h.h<<"],\"rows\":[";
  for(int i=0;i<h.rows.size();i++){if(i)cout<<",";cout<<h.rows[i];}cout<<"]}";
 }
 cout<<"]}\n";
 return 0;
}
