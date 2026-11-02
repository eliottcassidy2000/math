// Exhaustive endpoint certificates over an explicitly supplied graph6 universe.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct Graph {int n,a[11]{},d[11]{},m=0,w=0,p21=0;};
string quote(const string&s){string r="\"";for(char c:s){if(c=='\\'||c=='\"')r+='\\';r+=c;}return r+'\"';}
Graph decode(const string&s){
 Graph g;g.n=(unsigned char)s[0]-63;if(g.n<1||g.n>11)throw runtime_error("graph order");
 int b=0;for(int j=1;j<g.n;j++)for(int i=0;i<j;i++,b++)if(((unsigned char)s[1+b/6]-63)>>(5-b%6)&1){g.a[i]|=1<<j;g.a[j]|=1<<i;g.m++;}
 for(int i=0;i<g.n;i++){g.d[i]=__builtin_popcount((unsigned)g.a[i]);g.w+=g.d[i]*(g.d[i]-1)/2;}
 for(int i=0;i<g.n;i++)for(int j=i+1;j<g.n;j++)if(g.a[i]>>j&1)
  for(int v=0;v<g.n;v++)if(v!=i&&v!=j){int d=g.d[v]-((g.a[v]>>i)&1)-((g.a[v]>>j)&1);g.p21+=d*(d-1)/2;}
 return g;
}
array<int,6> matches(const Graph&g){
 array<array<int,6>,2048>dp{};dp[0][0]=1;
 for(int s=1;s<(1<<g.n);s++){int i=__builtin_ctz((unsigned)s),t=s^(1<<i);dp[s]=dp[t];
  for(int z=g.a[i]&t;z;z&=z-1){int j=__builtin_ctz((unsigned)z);for(int k=1;k<=5;k++)dp[s][k]+=dp[t^(1<<j)][k-1];}}
 return dp[(1<<g.n)-1];
}
Graph augment(const Graph&g,int s){
 Graph h=g;h.n++;for(int i=0;i<g.n;i++)if(s>>i&1){h.a[i]|=1<<g.n;h.a[g.n]|=1<<i;}
 h.m=h.w=h.p21=0;for(int i=0;i<h.n;i++){h.d[i]=__builtin_popcount((unsigned)h.a[i]);h.m+=h.d[i];h.w+=h.d[i]*(h.d[i]-1)/2;}h.m/=2;
 for(int i=0;i<h.n;i++)for(int j=i+1;j<h.n;j++)if(h.a[i]>>j&1)for(int v=0;v<h.n;v++)if(v!=i&&v!=j){int d=h.d[v]-((h.a[v]>>i)&1)-((h.a[v]>>j)&1);h.p21+=d*(d-1)/2;}return h;
}
int direct11(const Graph&g){auto m=matches(g);return 499*m[1]-320*m[2]+80*m[3]-20*m[4]+5*m[5]+80*g.p21+56*g.w;}
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("usage: graph_certificate triangle8|extension11 expected_count");
 string mode=argv[1],s;uint64_t expected=stoull(argv[2]),count=0,automatic=0,subsets=0,checks=0,zeros=0;int minimum=2147483647,maxneg=0;string witness;int witness_mask=0;
 map<int,uint64_t>neg_hist;vector<string>exceptional;
 while(getline(cin,s)){
  if(s.empty()||s[0]=='>')continue;Graph g=decode(s);count++;
  if(mode=="triangle8"){
   if(g.n!=8)throw runtime_error("expected n8");auto c=matches(g);int f=512+33*g.w-32*c[2]+7*c[3];
   if(f<minimum){minimum=f;witness=s;}if(f<0)throw runtime_error("negative triangle bracket");
  }else if(mode=="extension11"){
   if(g.n!=10)throw runtime_error("expected n10");int Z[1024];Z[0]=1;
   for(int S=1;S<1024;S++){int i=__builtin_ctz((unsigned)S),t=S^(1<<i);int z=2*Z[t];for(int b=g.a[i]&t;b;b&=b-1)z-=Z[t^(1<<__builtin_ctz((unsigned)b))];Z[S]=z;}
   int base=5120-5*Z[1023]-781*g.m+80*g.p21+56*g.w;
   int L[10],Q[10][10]{},negative[10],k=0,lower=base;
   for(int i=0;i<10;i++){
    int nd=0;for(int j=0;j<10;j++)if(g.a[i]>>j&1)nd+=g.d[j];
    int wd=g.w-g.d[i]*(g.d[i]-1)/2-nd+g.d[i];
    int ell=wd+g.d[i]*(g.m-g.d[i]+1)-nd;
    if(Z[1023^(1<<i)]%2)throw runtime_error("odd scaled nine polynomial");
    L[i]=5*(Z[1023^(1<<i)]/2)-781+80*ell+56*g.d[i];
    if(L[i]<0){negative[k++]=i;lower+=L[i];}
    for(int j=0;j<i;j++){Q[i][j]=Q[j][i]=80*(g.m-g.d[i]-g.d[j]+((g.a[i]>>j)&1))+56;if(Q[i][j]<56)throw runtime_error("outside-edge count");}
   }
   maxneg=max(maxneg,k);neg_hist[k]++;
   // Independent direct pattern kernel, all1024 neighborhoods on periodic bases.
   if(count<=20||count%100003==0){
    for(int mask=0;mask<1024;mask++){
     int f=base;for(int i=0;i<10;i++)if(mask>>i&1){f+=L[i];for(int j=0;j<i;j++)if(mask>>j&1)f+=Q[i][j];}
     if(f!=direct11(augment(g,mask)))throw runtime_error("extension identity mismatch");checks++;
    }
   }
   if(lower>0){automatic++;continue;}
   exceptional.push_back(s);
   int dp[1024];dp[0]=base;
   for(int mask=0;mask<(1<<k);mask++){
    int f;if(mask==0)f=base;else{int b=__builtin_ctz((unsigned)mask),t=mask^(1<<b),i=negative[b];f=dp[t]+L[i];for(int z=t;z;z&=z-1)f+=Q[i][negative[__builtin_ctz((unsigned)z)]];dp[mask]=f;}
    subsets++;if(f<minimum){minimum=f;witness=s;witness_mask=0;for(int j=0;j<k;j++)if(mask>>j&1)witness_mask|=1<<negative[j];}
    if(f==0)zeros++;if(f<0){cerr<<s<<" "<<witness_mask<<" "<<f<<"\n";throw runtime_error("negative extension");}
   }
  }else throw runtime_error("unknown mode");
  if(count%1000000==0)cerr<<"bases="<<count<<" automatic="<<automatic<<" min="<<minimum<<"\n";
 }
 if(count!=expected)throw runtime_error("incomplete universe count");
 cout<<"{\"status\":\"passed\",\"mode\":\""<<mode<<"\",\"base_graphs\":"<<count<<",\"minimum_checked\":"<<minimum<<",\"witness_graph6\":"<<quote(witness)<<",\"witness_neighborhood\":"<<witness_mask<<",\"strict_positive_lower_bound_bases\":"<<automatic<<",\"minimizer_subsets_checked\":"<<subsets<<",\"direct_extension_checks\":"<<checks<<",\"zero_subsets\":"<<zeros<<",\"max_negative_coordinates\":"<<maxneg<<",\"negative_coordinate_histogram\":{";
 bool first=true;for(auto[k,c]:neg_hist){if(!first)cout<<",";first=false;cout<<"\""<<k<<"\":"<<c;}cout<<"},\"exceptional_bases\":[";
 first=true;for(auto&g:exceptional){if(!first)cout<<",";first=false;cout<<quote(g);}cout<<"]}\n";
}catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
