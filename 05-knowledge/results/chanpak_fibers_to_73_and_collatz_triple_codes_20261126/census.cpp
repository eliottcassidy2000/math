#include <algorithm>
#include <array>
#include <cassert>
#include <bit>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
int n,U; vector<pair<int,int>> E;
using G=array<int,8>;
using P=array<int,5>; // triangles, 5-cycles, 7-cycles, triangle pairs, (3,5) pairs
int paths(const G&a){
 int d[256][8]={};
 for(int s=1;s<=U;s++)for(int b=s;b;b&=b-1){
   int i=countr_zero((unsigned)b),r=s^(1<<i),x=!r;
   for(int c=r&~a[i];c;c&=c-1)x+=d[r][countr_zero((unsigned)c)];
   d[s][i]=x;
 }
 int h=0;for(int i=0;i<n;i++)h+=d[U][i];return h;
}
pair<int,int> canon(const G&a){
 array<int,8> ord{},deg{};for(int i=0;i<n;i++)ord[i]=i,deg[i]=popcount((unsigned)a[i]);
 sort(ord.begin(),ord.begin()+n,[&](int i,int j){return pair{deg[i],i}<pair{deg[j],j};});
 vector<pair<int,int>> ranges;
 for(int i=0,j;i<n;i=j){j=i+1;while(j<n&&deg[ord[i]]==deg[ord[j]])j++;if(j>i+1)ranges.emplace_back(i,j);}
 int b=1<<28,aut=0;
 auto rec=[&](auto&&f,int l)->void{
  if(l==(int)ranges.size()){
   int x=0;for(int k=0;k<(int)E.size();k++){auto[i,j]=E[k];x|=((a[ord[i]]>>ord[j])&1)<<k;}
   if(x<b)b=x,aut=0;if(x==b)aut++;return;
  }
  auto [i,j]=ranges[l];do{f(f,l+1);}while(next_permutation(ord.begin()+i,ord.begin()+j));
 };
 rec(rec,0);return {b,aut};
}
P profile(const G&a){
 int cy[256]={};
 auto r=[&](auto&&f,int x,int v,int s,int l)->void{
  if((l==3||l==5||l==7)&&(a[v]>>x&1))cy[s]++;
  if(l<7)for(int b=a[v]&U&~s&~((1<<x)-1);b;b&=b-1){int j=countr_zero((unsigned)b);f(f,x,j,s|(1<<j),l+1);}
 };
 for(int i=0;i<n;i++)r(r,i,i,1<<i,1);
 P p{};vector<int> t,f;
 for(int s=1;s<=U;s++)if(cy[s]){
  int k=popcount((unsigned)s);p[(k-3)/2]+=cy[s];
  if(k==3)t.push_back(s);if(k==5)f.push_back(s);
 }
 for(int i=0;i<(int)t.size();i++){
  for(int j=i+1;j<(int)t.size();j++)if(!(t[i]&t[j]))p[3]+=cy[t[i]]*cy[t[j]];
  for(int j:f)if(!(t[i]&j))p[4]+=cy[t[i]]*cy[j];
 }
 return p;
}
int main(){
 for(n=3;n<=8;n++){
  U=(1<<n)-1;E.clear();G a{};vector<pair<int,int>> chords;map<int,set<int>> seen;
  for(int i=0;i<n;i++)a[i]=1<<((i+1)%n);
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   E.emplace_back(i,j);
   if(j!=i+1&&!(i==0&&j==n-1))chords.emplace_back(i,j),a[i]|=1<<j;
  }
  int last=0;
  for(int it=0;it<(1<<(int)chords.size());it++){
   const int m=it^(it>>1);
   if(it){auto [i,j]=chords[countr_zero((unsigned)(m^last))];a[i]^=1<<j;a[j]^=1<<i;}last=m;
   int h=paths(a);if(h>73)continue;
   auto [k,aut]=canon(a);if(!seen[h].insert(k).second)continue;
   G b{};for(int i=0;i<n;i++)b[i]=U^a[i]^(1<<i);
   bool self=(canon(b).first==k);
   const auto p=profile(a);assert(h==1+2*(p[0]+p[1]+p[2])+4*(p[3]+p[4]));
   cout<<n<<' '<<h<<' '<<m<<' '<<aut<<' '<<self;
   for(int x:p)cout<<' '<<x; cout<<"\n";
  }
 }
}
