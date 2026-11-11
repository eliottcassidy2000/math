// Local one/two-arc patch search on the literal 6853 source.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using I=long long;
const int n=15,full=(1<<n)-1;
vector<I> ends(const vector<int>& rows){
 vector<I> d((1<<n)*n);
 for(int v=0;v<n;v++)d[((1<<v)*n)+v]=1;
 for(int m=1;m<=full;m++)for(int v=0;v<n;v++)if(I x=d[m*n+v]){
  int a=rows[v]&~m&full;
  while(a){int b=a&-a;a-=b;int w=__builtin_ctz((unsigned)b);d[(m|b)*n+w]+=x;}
 }
 return d;
}
vector<int> converse(const vector<int>& rows){
 vector<int> out(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(rows[j]>>i&1)out[i]|=1<<j;
 return out;
}
I total(const vector<I>& d){I h=0;for(int v=0;v<n;v++)h+=d[full*n+v];return h;}
I passage(const vector<I>& d,const vector<I>& r,int u,int v){
 int rem=full^(1<<u)^(1<<v);I sum=0;
 for(int s=rem;;s=(s-1)&rem){int left=s|(1<<u);sum+=d[left*n+u]*r[(full^left)*n+v];if(!s)break;}
 return sum;
}
void flip(vector<int>& rows,int u,int v){rows[u]^=1<<v;rows[v]^=1<<u;}
vector<I> runs(const vector<int>& rows,int star){
 const int z=n-1,N=1<<z,D=n+1;vector<int> verts;
 for(int v=0;v<n;v++)if(v!=star)verts.push_back(v);
 vector<I> dp((size_t)N*z*D);
 auto at=[&](int mask,int v,int j)->I&{return dp[((size_t)mask*z+v)*D+j];};
 for(int v=0;v<z;v++){at(1<<v,v,0)=1;if(rows[star]>>verts[v]&1)at(1<<v,v,1)=1;}
 for(int mask=1;mask<N;mask++)for(int u=0;u<z;u++)if(mask>>u&1){
  int rest=(N-1)^mask;
  while(rest){int b=rest&-rest;rest-=b;int v=__builtin_ctz((unsigned)b);
   bool direct=rows[verts[u]]>>verts[v]&1;
   bool via=(rows[verts[u]]>>star&1)&&(rows[star]>>verts[v]&1);
   if(!direct&&!via)continue;
   for(int j=0;j<D;j++)if(I x=at(mask,u,j)){
    if(direct)at(mask|b,v,j)+=x;
    if(via&&j+1<D)at(mask|b,v,j+1)+=x;
   }
  }
 }
 vector<I> ans(D);
 for(int v=0;v<z;v++)for(int j=0;j<D;j++)if(I x=at(N-1,v,j)){
  ans[j]+=x;if((rows[verts[v]]>>star&1)&&j+1<D)ans[j+1]+=x;
 }
 return ans;
}
// Complete ordered-run coefficients for H,E,V,d,K at a marked nonstar root.
vector<vector<I>> port_runs(const vector<int>& rows,int star,int root){
 const int z=n-1,N=1<<z,D=n+2;vector<int> verts;
 for(int v=0;v<n;v++)if(v!=star)verts.push_back(v);
 int rt=find(verts.begin(),verts.end(),root)-verts.begin();
 vector<I> dp((size_t)N*z*D),P((size_t)N*D),E((size_t)N*D);
 auto at=[&](int mask,int v,int j)->I&{return dp[((size_t)mask*z+v)*D+j];};
 for(int v=0;v<z;v++){at(1<<v,v,0)=1;if(rows[star]>>verts[v]&1)at(1<<v,v,1)=1;}
 for(int mask=1;mask<N;mask++)for(int u=0;u<z;u++)if(mask>>u&1){
  int rest=(N-1)^mask;
  while(rest){int b=rest&-rest;rest-=b;int v=__builtin_ctz((unsigned)b);
   bool direct=rows[verts[u]]>>verts[v]&1,via=(rows[verts[u]]>>star&1)&&(rows[star]>>verts[v]&1);
   if(!direct&&!via)continue;
   for(int j=0;j<D;j++)if(I x=at(mask,u,j)){
    if(direct)at(mask|b,v,j)+=x;if(via&&j+1<D)at(mask|b,v,j+1)+=x;
   }
  }
 }
 P[0]=P[1]=1; // Empty outside set: empty path, or one placeholder run.
 for(int mask=1;mask<N;mask++)for(int v=0;v<z;v++)for(int j=0;j<D;j++)if(I x=at(mask,v,j)){
  P[(size_t)mask*D+j]+=x;
  if((rows[verts[v]]>>star&1)&&j+1<D)P[(size_t)mask*D+j+1]+=x;
  if(v==rt)E[(size_t)mask*D+j]+=x;
 }
 vector<vector<I>> out(5,vector<I>(D));int all=N-1,deleted=all^(1<<rt);
 for(int j=0;j<D;j++){out[0][j]=P[(size_t)all*D+j];out[1][j]=E[(size_t)all*D+j];out[3][j]=P[(size_t)deleted*D+j];}
 for(int mask=0;mask<N;mask++)if(mask>>rt&1){int other=all^mask;
  for(int a=0;a<D;a++)if(I x=E[(size_t)mask*D+a])for(int b=0;a+b<D;b++)out[2][a+b]+=x*P[(size_t)other*D+b];
 }
 for(int mask=deleted;;mask=(mask-1)&deleted){int other=deleted^mask;
  for(int a=0;a<D;a++)if(I x=P[(size_t)mask*D+a])for(int b=0;a+b<D;b++)out[4][a+b]+=x*P[(size_t)other*D+b];
  if(!mask)break;
 }
 return out;
}
void vec(const vector<int>& a){cout<<'[';for(size_t i=0;i<a.size();i++){if(i)cout<<',';cout<<a[i];}cout<<']';}
void vec(const vector<I>& a){cout<<'[';for(size_t i=0;i<a.size();i++){if(i)cout<<',';cout<<a[i];}cout<<']';}
vector<I> marked(const vector<I>& d,const vector<I>& r,int root){
 vector<I> h(1<<n);h[0]=1;
 for(int mask=1;mask<=full;mask++)for(int v=0;v<n;v++)h[mask]+=d[mask*n+v];
 I U=0,V=0,K=0;int deletion=full^(1<<root);
 for(int mask=1;mask<=full;mask++)if(mask>>root&1){U+=r[mask*n+root]*h[full^mask];V+=d[mask*n+root]*h[full^mask];}
 for(int s=deletion;;s=(s-1)&deletion){K+=h[s]*h[deletion^s];if(!s)break;}
 return {h[full],d[full*n+root],r[full*n+root],U,V,h[deletion],K};
}
struct Near {vector<int> rows,edges;I H;};
int main(){
 vector<int> rows={162,1020,1017,1009,993,960,897,768,513,1};
 int root=0;
 for(int birth=0;birth<5;birth++){
  int v=rows.size(),mask=(1<<v)-1;rows[root]|=1<<v;rows.push_back(mask^(1<<root));root=v;
 }
 auto base=ends(rows);I H=total(base);
 cout<<"{\"source_rows\":";vec(rows);cout<<",\"source_H\":"<<H<<",\"placeholder\":7,\"root\":14,\"target\":5993,\"hits\":[";
 vector<pair<int,int>> edges;for(int u=0;u<n-1;u++)for(int v=u+1;v<n-1;v++)edges.emplace_back(u,v);
 int hits=0,singles=0,pairs=0;I nearest=1LL<<60,bestvalue=0;vector<int> best,bestrows;
 vector<Near> near;
 auto report=[&](const vector<int>& graph,vector<int> changed,I value){
  I distance=llabs(value-5993);if(distance<nearest){nearest=distance;best=changed;bestvalue=value;bestrows=graph;}
  if(distance==2)near.push_back({graph,changed,value});
  if(value!=5993)return;
  if(hits++)cout<<',';
  cout<<"{\"changed_endpoints\":";vec(changed);cout<<",\"rows\":";vec(graph);
  cout<<",\"placeholder_run_coefficients\":";vec(runs(graph,7));cout<<'}';
 };
 for(size_t i=0;i<edges.size();i++){
  auto [u,v]=edges[i];auto first=rows;flip(first,u,v);
  auto d=ends(first),r=ends(converse(first));I h=total(d);singles++;
  report(first,{u,v},h);
  for(size_t j=i+1;j<edges.size();j++){
   auto [a,b]=edges[j];int tail=(first[a]>>b&1)?a:b,head=tail==a?b:a;
   I after=h-passage(d,r,tail,head)+passage(d,r,head,tail);pairs++;
   if(after==5993){auto second=first;flip(second,a,b);report(second,{u,v,a,b},after);}
   else {
    if(llabs(after-5993)==2){auto second=first;flip(second,a,b);near.push_back({second,{u,v,a,b},after});}
    if(llabs(after-5993)<nearest){nearest=llabs(after-5993);best={u,v,a,b};bestvalue=after;bestrows=first;flip(bestrows,a,b);}
   }
  }
 }
 // Move the sole in-neighbor while retaining indegree one at the marked root.
 int relocations=0;
 for(int pred=0;pred<n-2;pred++){
  auto graph=rows;flip(graph,n-2,n-1);flip(graph,pred,n-1);relocations++;
  report(graph,{n-2,n-1,pred,n-1},total(ends(graph)));
 }
 cout<<"],\"single_internal_reversals\":"<<singles<<",\"double_internal_reversals\":"<<pairs;
 cout<<",\"indegree_one_relocations\":"<<relocations<<",\"hit_count\":"<<hits;
 cout<<",\"nearest_absolute_error\":"<<nearest<<",\"nearest_H\":"<<bestvalue<<",\"nearest_changed_endpoints\":";vec(best);
 cout<<",\"nearest_rows\":";vec(bestrows);cout<<",\"nearest_run_coefficients\":";vec(runs(bestrows,7));
 cout<<",\"nearest_marked_HESUVdK\":";
 auto bd=ends(bestrows),br=ends(converse(bestrows));
 vec(marked(bd,br,14));
 cout<<",\"all_distance_two_patches\":[";
 for(size_t ni=0;ni<near.size();ni++){
  if(ni)cout<<',';auto &nn=near[ni];auto nd=ends(nn.rows),nr=ends(converse(nn.rows));
  cout<<"{\"changed_endpoints\":";vec(nn.edges);cout<<",\"H\":"<<nn.H<<",\"marked_HESUVdK\":";vec(marked(nd,nr,14));
  auto pr=port_runs(nn.rows,7,14),cr=port_runs(converse(nn.rows),7,14);
  cout<<",\"rows\":";vec(nn.rows);cout<<",\"run_coefficients_HESUVdK\":[";
  vector<vector<I>> allpr={pr[0],pr[1],cr[1],cr[2],pr[2],pr[3],pr[4]};
  for(size_t p=0;p<allpr.size();p++){if(p)cout<<',';vec(allpr[p]);}cout<<']';
  cout<<",\"extra_edge_hits\":[";int extra=0;
  for(auto [u,v]:edges){
   int tail=(nn.rows[u]>>v&1)?u:v,head=tail==u?v:u;
   I value=nn.H-passage(nd,nr,tail,head)+passage(nd,nr,head,tail);
   if(value!=5993)continue;
   auto graph=nn.rows;flip(graph,u,v);if(extra++)cout<<',';
   cout<<"{\"extra_edge\":["<<u<<','<<v<<"],\"rows\":";vec(graph);
   cout<<",\"run_coefficients\":";vec(runs(graph,7));cout<<'}';
  }
  cout<<"]}";
 }
 cout<<"],\"distance_two_patch_count\":"<<near.size()<<",\"extra_edge_controls\":"<<edges.size()*near.size()<<"}\n";
}
