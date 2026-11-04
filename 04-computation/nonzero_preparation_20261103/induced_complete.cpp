// RESERVED pending independent audit. Fixed-source, exhaustive induced extension certificate.
// Compile: c++ -std=c++17 -O3 induced_complete.cpp -o /tmp/induced_complete
// Run: /tmp/induced_complete > induced_complete.json
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
using Count=uint64_t;
using Graph=vector<unsigned>;
static Count gates=0;
static void need(bool value,const char* label){++gates;if(!value){cerr<<label<<'\n';exit(1);}}
static const Graph SOURCE={1790,1668,1672,1778,1670,1942,1718,1280,1631,128,512};
static const int ROOT=8, MAX_BIRTHS=5;
static const Count SOURCE_H=1429,SOURCE_F=1563,CEILING=4291;
static Count visited[6]={},retained[6]={},raw_hits[6][3]={},strong_hits[6][3]={};
static uint64_t digest=14695981039346656037ULL;
static void mix(uint64_t v){for(int j=0;j<8;j++){digest^=(v>>(8*j))&255;digest*=1099511628211ULL;}}
static bool tournament(const Graph&a){
 for(int i=0;i<(int)a.size();i++){
  if((a[i]>>i&1)||(a[i]>>a.size()))return false;
  for(int j=0;j<i;j++)if(((a[i]>>j)&1)+((a[j]>>i)&1)!=1)return false;
 }return true;
}
static Graph birth(Graph a,unsigned incoming){
 int n=a.size();need(incoming<(1u<<n),"invalid birth cut");
 for(int i=0;i<n;i++)if(incoming>>i&1)a[i]|=1u<<n;
 a.push_back(((1u<<n)-1)^incoming);return a;
}
static bool strong(const Graph&a){
 unsigned full=(1u<<a.size())-1;
 for(int i=0;i<(int)a.size();i++){
  unsigned seen=1u<<i,old=0;
  while(old!=seen){old=seen;for(int v=0;v<(int)a.size();v++)if(seen>>v&1)seen|=a[v];}
  if(seen!=full)return false;
 }return true;
}
struct Observers{Count h,f;};
// Independent direction: append the next vertex to an existing directed path.
static Observers forward(const Graph&a,int root){
 int n=a.size();need(n>0&&n<=16,"forward universe");unsigned size=1u<<n,full=size-1;
 vector<Count>E(size*n,0),P(size,0);P[0]=1;
 for(int i=0;i<n;i++)E[(1u<<i)*n+i]=1;
 for(unsigned mask=1;mask<size;mask++)for(unsigned vs=mask;vs;){
  unsigned bit=vs&-vs;vs-=bit;int v=__builtin_ctz(bit);Count value=E[mask*n+v];
  P[mask]+=value;unsigned next=a[v]&(full^mask);
  while(next){unsigned b=next&-next;next-=b;E[(mask|b)*n+__builtin_ctz(b)]+=value;}
 }
 Count f=0;if(root>=0)for(unsigned mask=1;mask<size;mask++)f+=E[mask*n+root]*P[full^mask];
 return {P[full],f};
}
struct Response{int n;Count h;vector<Count> e,s,K;};
// Endpoint deletion computes all subgraph path endpoints and starts.
static Response response(const Graph&a){
 int n=a.size();need(n>0&&n<=15,"response universe");unsigned size=1u<<n,full=size-1;
 vector<Count>E(size*n,0),S(size*n,0);vector<unsigned> in(n,0);
 for(int v=0;v<n;v++){
  E[(1u<<v)*n+v]=S[(1u<<v)*n+v]=1;
  for(int w=0;w<n;w++)if(a[w]>>v&1)in[v]|=1u<<w;
 }
 for(unsigned mask=1;mask<size;mask++)for(unsigned vs=mask;vs;){
  unsigned bit=vs&-vs;vs-=bit;int v=__builtin_ctz(bit);unsigned rest=mask^bit;if(!rest)continue;
  for(unsigned z=rest&in[v];z;){unsigned b=z&-z;z-=b;E[mask*n+v]+=E[rest*n+__builtin_ctz(b)];}
  for(unsigned z=rest&a[v];z;){unsigned b=z&-z;z-=b;S[mask*n+v]+=S[rest*n+__builtin_ctz(b)];}
 }
 Response q{n,0,vector<Count>(n),vector<Count>(n),vector<Count>(n*n,0)};
 for(int v=0;v<n;v++){q.e[v]=E[full*n+v];q.s[v]=S[full*n+v];q.h+=q.e[v];}
 need(q.h==accumulate(q.s.begin(),q.s.end(),Count(0)),"start/end Hamilton equality");
 for(unsigned mask=1;mask<full;mask++)for(unsigned vs=mask;vs;){
  unsigned bit=vs&-vs;vs-=bit;int i=__builtin_ctz(bit);Count e=E[mask*n+i];if(!e)continue;
  for(unsigned ws=full^mask;ws;){unsigned b=ws&-ws;ws-=b;int j=__builtin_ctz(b);q.K[i*n+j]+=e*S[(full^mask)*n+j];}
 }
 return q;
}
static Count direct_cut(const Response&q,unsigned mask){
 Count h=0;
 for(int i=0;i<q.n;i++)if(mask>>i&1){h+=q.e[i];for(int j=0;j<q.n;j++)if(!(mask>>j&1))h+=q.K[i*q.n+j];}else h+=q.s[i];
 return h;
}
static vector<Count> cuts(const Response&q){
 int n=q.n;unsigned size=1u<<n;vector<int64_t>A(n,0),B(n*n,0);vector<Count> values(size,0);
 values[0]=q.h;
 for(int i=0;i<n;i++){
  A[i]=(int64_t)q.e[i]-(int64_t)q.s[i];
  for(int j=0;j<n;j++){A[i]+=q.K[i*n+j];B[i*n+j]=q.K[i*n+j]+q.K[j*n+i];}
 }
 for(unsigned mask=1;mask<size;mask++){
  unsigned bit=mask&-mask;int i=__builtin_ctz(bit);unsigned rest=mask^bit;
  int64_t value=(int64_t)values[rest]+A[i];
  for(unsigned z=rest;z;){unsigned b=z&-z;z-=b;value-=B[i*n+__builtin_ctz(b)];}
  need(value>=(int64_t)q.h,"induced Hamilton monotonicity");values[mask]=(Count)value;
 }
 need(values.back()==q.h,"universal sink boundary");return values;
}
struct Hit{int depth,digit;Count h,f;vector<unsigned> sequence;Graph graph;};
static vector<Hit> hits;
static bool clone_target(const Graph&a){
 if(a.size()!=13)return false;
 unsigned block=(1u<<10)|(1u<<11)|(1u<<12);
 for(int v:{10,11,12})for(int w=0;w<13;w++)if(!(block>>w&1)){
  bool expected=(w==9);if(bool(a[v]>>w&1)!=expected)return false;
 }
 // A tournament on three vertices is transitive exactly when its outdegrees differ.
 vector<int> degrees;for(int v:{10,11,12})degrees.push_back(__builtin_popcount(a[v]&block));
 sort(degrees.begin(),degrees.end());return degrees==vector<int>({0,1,2});
}
static void search(const Graph&a,const vector<unsigned>&sequence,int depth,Count known_h){
 if(depth==MAX_BIRTHS)return;
 Response q=response(a);need(q.h==known_h,"retained parent count");auto values=cuts(q);
 for(unsigned incoming=0;incoming<values.size();incoming++){
  ++visited[depth+1];Count h=values[incoming];if(h>CEILING)continue;
  ++retained[depth+1];Graph child=birth(a,incoming);auto word=sequence;word.push_back(incoming);
  mix(depth+1);mix(h);for(unsigned row:child)mix(row);
  if(h==3*SOURCE_H||h==3*SOURCE_H+2||h==3*SOURCE_H+4){
   int digit=(int)(h-3*SOURCE_H);++raw_hits[depth+1][digit/2];
   if(strong(child)){
    ++strong_hits[depth+1][digit/2];auto o=forward(child,ROOT);need(o.h==h,"independent hit Hamilton count");
    need(tournament(child),"actual hit tournament");
    for(int v=0;v<11;v++)need((child[v]&2047)==SOURCE[v],"actual old induced source");
    if(digit==4&&depth+1==1)need(incoming==289&&o.f==4677,"one-birth clock coincidence");
    if(digit==4&&depth+1==2)need(clone_target(child)&&o.f==3*SOURCE_F,"exact clone-only pair hit");
    hits.push_back({depth+1,digit,h,o.f,word,child});
   }
  }
  search(child,word,depth+1,h);
 }
}
static void controls(){
 need(tournament(SOURCE)&&strong(SOURCE),"literal source shape");auto obs=forward(SOURCE,ROOT);
 need(obs.h==SOURCE_H&&obs.f==SOURCE_F,"literal source observers");
 need(9*625==5625&&5625>CEILING,"cited strong order17 cutoff arithmetic");
 Count cases=0;
 for(int n=1;n<=4;n++){
  int edges=n*(n-1)/2;
  for(unsigned code=0;code<(1u<<edges);code++){
   Graph a(n,0);int e=0;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,e++)if(code>>e&1)a[i]|=1u<<j;else a[j]|=1u<<i;
   auto q=response(a);auto v=cuts(q);
   for(unsigned incoming=0;incoming<(1u<<n);incoming++){
    need(v[incoming]==direct_cut(q,incoming),"quadratic versus slot cut law");
    need(v[incoming]==forward(birth(a,incoming),-1).h,"all-small literal birth count");++cases;
   }
  }
 }
 need(cases==1098,"complete small-control universe");
}
static void print_vector(const vector<unsigned>&v){cout<<'[';for(size_t i=0;i<v.size();i++)cout<<(i?",":"")<<v[i];cout<<']';}
int main(){
 controls();search(SOURCE,{},0,SOURCE_H);
 const Count expected_visits[]={0,2048,151552,1343488,12550144,136839168};
 const Count expected_retained[]={0,37,164,766,4176,26520};
 for(int d=1;d<=5;d++){
  need(visited[d]==expected_visits[d]&&retained[d]==expected_retained[d],"declared level universe");
  need(raw_hits[d][1]==0,"no plus2 even in weak range");
  need(strong_hits[d][2]==(d==1?1:d==2?6:0),"complete strong plus4 classification");
  need(strong_hits[d][0]==(d==2?4:0),"zero positive-control classification");
 }
 cout<<"{\n  \"status\":\"FINITE-EXACT; analytic all-order consumer requires cited strong minimum\",\n  \"source\":";print_vector(SOURCE);
 cout<<",\n  \"root\":"<<ROOT<<",\n  \"H\":"<<SOURCE_H<<",\n  \"F\":"<<SOURCE_F<<",\n  \"ceiling\":"<<CEILING<<",\n  \"small_birth_controls\":1098,\n  \"levels\":[\n";
 for(int d=1;d<=5;d++){
  cout<<(d>1?",\n":"")<<"    {\"births\":"<<d<<",\"raw_universe_power_of_two\":"<<11*d+d*(d-1)/2<<",\"tested_cuts\":"<<visited[d]<<",\"retained\":"<<retained[d]<<",\"all_target_hits_0_2_4\":[";
  for(int j=0;j<3;j++)cout<<(j?",":"")<<raw_hits[d][j];cout<<"],\"strong_target_hits_0_2_4\":[";
  for(int j=0;j<3;j++)cout<<(j?",":"")<<strong_hits[d][j];cout<<"]}";
 }
 cout<<"\n  ],\n  \"retained_stream_fnv1a64\":"<<digest<<",\n  \"hits\":[\n";
 for(size_t j=0;j<hits.size();j++){
  auto &v=hits[j];cout<<(j?",\n":"")<<"    {\"births\":"<<v.depth<<",\"digit\":"<<v.digit<<",\"H\":"<<v.h<<",\"F\":"<<v.f<<",\"cuts\":";print_vector(v.sequence);cout<<",\"rows\":";print_vector(v.graph);cout<<'}';
 }
 cout<<"\n  ],\n  \"gates\":"<<gates<<"\n}\n";
}
