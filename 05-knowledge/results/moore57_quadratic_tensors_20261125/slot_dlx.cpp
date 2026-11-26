#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using namespace std;
struct Row {int h,y; array<int,3> slots;};
struct DLX {
 struct Node {int l,r,u,d,c,id;};
 vector<Node> a;
 vector<int>s,path;
 uint64_t visits=0,limit;
 DLX(int n,uint64_t b):a(n+1),s(n+1),limit(b) {
  for(int i=0;i<=n;i++)a[i]={(i+n)%(n+1),(i+1)%(n+1),i,i,i,-1};
 }
 void add(int id,const vector<int>& cols) {
  int first=a.size();
  for(int i=0;i<(int)cols.size();i++){
   int c=cols[i]+1,p=a.size(),u=a[c].u;
   a.push_back({i?p-1:first+(int)cols.size()-1,i+1==(int)cols.size()?first:p+1,u,c,c,id});
   a[u].d=p;a[c].u=p;s[c]++;
  }
 }
 void cover(int c){
  a[a[c].l].r=a[c].r;a[a[c].r].l=a[c].l;
  for(int i=a[c].d;i!=c;i=a[i].d)for(int j=a[i].r;j!=i;j=a[j].r){
   a[a[j].u].d=a[j].d;a[a[j].d].u=a[j].u;--s[a[j].c];
  }
 }
 void uncover(int c){
  for(int i=a[c].u;i!=c;i=a[i].u)for(int j=a[i].l;j!=i;j=a[j].l){
   a[a[j].u].d=j;a[a[j].d].u=j;++s[a[j].c];
  }
  a[a[c].l].r=c;a[a[c].r].l=c;
 }
 bool go() {
  if(++visits>limit)return false;
  if(!a[0].r)return true;
  int c=a[0].r;
  for(int i=a[c].r;i;i=a[i].r)if(s[i]<s[c])c=i;
  if(!s[c])return false;
  cover(c);
  for(int i=a[c].d;i!=c && visits<=limit;i=a[i].d){
   path.push_back(a[i].id);
   for(int j=a[i].r;j!=i;j=a[j].r)cover(a[j].c);
   if(go())return true;
   for(int j=a[i].l;j!=i;j=a[j].l)uncover(a[j].c);
   path.pop_back();
  }
  uncover(c);return false;
 }
};
int main(int argc,char**argv) {
 string mode=argc>1?argv[1]:"full";
 int runs=argc>2?stoi(argv[2]):1000;
 mt19937 gen(argc>3?stoi(argv[3]):1);
 uint64_t limit=argc>4?stoull(argv[4]):30000;
 int n;cin>>n; if(n!=53)return 2;
 vector<Row>C; Row mark;vector<int>ref(n);
 for(int i=0;i<n;i++){
  int k;cin>>ref[i]>>k;
  for(int z=0;z<k;z++){
   Row e;e.h=i;cin>>e.y;for(int &v:e.slots)cin>>v;
   if(i==0)mark=e;
   else if(e.slots[0]!=mark.slots[0]&&e.slots[1]!=mark.slots[1]&&e.slots[2]!=mark.slots[2])C.push_back(e);
  }
 }
 if(C.size()!=1741)return 3;
 array<vector<int>,3> labels;
 for(int k=0;k<3;k++){for(auto e:C)labels[k].push_back(e.slots[k]);auto&L=labels[k];sort(L.begin(),L.end());L.erase(unique(L.begin(),L.end()),L.end());}
 auto solve = [&](int drop,int omit,int forced,bool prefer) {
  int ids[4][56],m=0;
  for(int h=1;h<n;h++)ids[0][h]=(h==omit)?-1:m++;
  for(int k=0;k<3;k++)if(k!=drop)for(int y:labels[k])ids[k+1][y]=m++;
  int slack=m;if(omit>0)m+=3;
  vector<int> perm(m);iota(perm.begin(),perm.end(),0);shuffle(perm.begin(),perm.end(),gen);
  DLX D(m,limit);vector<int> order(C.size());iota(order.begin(),order.end(),0);shuffle(order.begin(),order.end(),gen);
  if(prefer)stable_partition(order.begin(),order.end(),[&](int e){return C[e].y==ref[C[e].h];});
  auto row=[&](int e){vector<int>v{perm[ids[0][C[e].h]]};for(int k=0;k<3;k++)if(k!=drop)v.push_back(perm[ids[k+1][C[e].slots[k]]]);return v;};
  if(omit>0)for(int k=0;k<3;k++){auto L=labels[k];shuffle(L.begin(),L.end(),gen);for(int y:L)D.add(-1,{perm[slack+k],perm[ids[k+1][y]]});}
  for(int e:order)if(C[e].h!=omit){D.add(e,row(e));}
  if(forced>=0){D.path.push_back(forced);for(int c:row(forced))D.cover(c+1);}
  if(!D.go())return vector<int>{};
  vector<int>ans(n,-1);ans[0]=mark.y;for(int e:D.path)if(e>=0)ans[C[e].h]=C[e].y;
  return ans;
 };
 if(mode=="project"){
  for(int d=0;d<3;d++)for(int k=0;k<(int)C.size();k++){
   vector<int> v;for(int t=0;t<runs&&v.empty();t++)v=solve(d,-1,k,true);
   if(v.empty()){cerr<<"unfinished "<<d<<" "<<k<<endl;return 4;}
   cout<<d<<" "<<k;for(int y:v)cout<<" "<<y;cout<<"\n";
  }
  return 0;
 }
 if(mode=="omit-all"){
  for(int h=1;h<53;h++){
   vector<int> v;for(int t=0;t<runs&&v.empty();t++)v=solve(-1,h,-1,t%4!=3);
   if(v.empty()){cerr<<"unfinished omit "<<h<<endl;continue;}
   cout<<h;for(int y:v)cout<<" "<<y;cout<<endl;
  }
  return 0;
 }
 if(mode=="near"){
  for(int t=0;t<runs;t++){auto v=solve(-1,t%2==0?11:30,-1,t%4!=3);if(!v.empty()){for(int y:v)cout<<y<<" ";cout<<endl;return 0;}}
 }else for(int t=0;t<runs;t++){auto v=solve(-1,-1,-1,t%4!=3);if(!v.empty()){for(int y:v)cout<<y<<" ";cout<<endl;return 0;}}
 return 1;
}
