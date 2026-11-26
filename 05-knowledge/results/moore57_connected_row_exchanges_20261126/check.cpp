#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
using U=uint64_t;
const U ALL=(1ull<<54)-1;
const string alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
using Row=array<int,6>;
using Set=array<int,5>;
struct Search{
 U allow[6][6][54]{}, freec[6][6][54]{};
 int owner[6][6][54][54];
 int cap;
 vector<Row> rows, candidates;
 array<uint64_t,6> inventory{},keyCount{},tested{};
 unordered_map<U,vector<int>> bank;
 vector<vector<U>> inverted;
 unordered_set<U> seen;
 Row work{};
 Search(string board,string witness,int cap):cap(cap){
  ifstream b(board),w(witness); for(auto &a:owner)for(auto &b:a)for(auto &c:b)fill_n(c,54,-1);
  for(int i=0;i<6;i++)for(int j=i+1;j<6;j++)for(int x=0;x<54;x++){
   string s;b>>s;allow[i][j][x]=freec[i][j][x]=stoull(s,nullptr,16);
  }
  string s;w>>s;
  if(s.size()%6)throw runtime_error("bad witness length");
  for(size_t t=0;t<s.size();t+=6){
   Row r;
   for(int i=0;i<6;i++){auto n=alphabet.find(s[t+i]);if(n>=54)throw runtime_error("bad symbol");r[i]=int(n);}
   for(int i=0;i<6;i++)for(int j=i+1;j<6;j++){
    if(!((freec[i][j][r[i]]>>r[j])&1))throw runtime_error("duplicate or forbidden");
    freec[i][j][r[i]]^=1ull<<r[j];owner[i][j][r[i]][r[j]]=rows.size();
   }
   rows.push_back(r);
  }
  if(rows.size()>=4095)throw runtime_error("key overflow");
 }
 U key(const Set& s,int k){
  U z=0;
  for(int i=k-1;i>=0;i--)z=(z<<12)|(s[i]+1);
  return z;
 }
 Set decode(U z,int& k){
  Set a{};k=0;
  while(z){a[k++]=int(z&4095)-1;z>>=12;}
  return a;
 }
 bool push(Set& s,int& k,int x){
  if(x<0)return true;
  int i=0;while(i<k&&s[i]<x)i++;
  if(i<k&&s[i]==x)return true;
  if(k==cap)return false;
  for(int j=k;j>i;j--)s[j]=s[j-1];
  s[i]=x;k++;return true;
 }
 void enumerate(int j,const Set& set,int size){
  if(j==6){
   if(size==1&&rows[set[0]]==work)return;
   inventory[size]++;
   bank[key(set,size)].push_back(candidates.size());
   candidates.push_back(work);return;
  }
  U opts=ALL;
  for(int i=0;i<j;i++)opts&=allow[i][j][work[i]];
  if(size==cap){
   U existing=0;
   for(int i=0;i<size;i++)existing|=1ull<<rows[set[i]][j];
   for(int i=0;i<j;i++)opts&=freec[i][j][work[i]]|existing;
  }
  while(opts){
   int x=__builtin_ctzll(opts);opts&=opts-1;work[j]=x;
   Set s=set;int k=size;bool ok=true;
   for(int i=0;i<j;i++)if(!push(s,k,owner[i][j][work[i]][x])){ok=false;break;}
   if(ok)enumerate(j+1,s,k);
  }
 }
 bool compatible(int x,int y){
  int eq=0;
  for(int j=0;j<6;j++)if(candidates[x][j]==candidates[y][j]&&++eq==2)return false;
  return true;
 }
 bool clique(const vector<int>& p,int n,vector<int>& chosen){
  if(!n)return true;
  for(size_t i=0;i+n<=p.size();i++){
   chosen.push_back(p[i]);
   if(n==1)return true;
   vector<int> q;
   for(size_t j=i+1;j<p.size();j++)if(compatible(p[i],p[j]))q.push_back(p[j]);
   if(clique(q,n-1,chosen))return true;
   chosen.pop_back();
  }
  return false;
 }
 void output(const Row& r){for(int x:r)cout<<alphabet[x];cout<<' ';}
 bool check(const Set& a,int n){
  vector<int> pool,chosen;
  for(int mask=0;mask<(1<<n);mask++){
   Set q{};int k=0;
   for(int i=0;i<n;i++)if(mask>>i&1)q[k++]=a[i];
   auto it=bank.find(key(q,k));
   if(it!=bank.end())pool.insert(pool.end(),it->second.begin(),it->second.end());
  }
  if(clique(pool,n+1,chosen)){
   cout<<"GAIN ";
   for(int i=0;i<n;i++)output(rows[a[i]]);
   cout<<"> ";
   for(int i:chosen)output(candidates[i]);
   cout<<'\n';return true;
  }
  return false;
 }
 bool visit(U z){
  if(!seen.insert(z).second)return false;
  int n;Set a=decode(z,n);tested[n]++;
  if(check(a,n))return true;
  if(n==cap)return false;
  for(int j=0;j<n;j++)for(U w:inverted[a[j]]){
   int m;Set b=decode(w,m);
   // A full-size support cannot grow any proper support of an earlier test.
   if(m==cap)continue;
   Set u=a;int k=n;bool ok=true;
   for(int i=0;i<m;i++)if(!push(u,k,b[i])){ok=false;break;}
   if(ok&&k>n&&visit(key(u,k)))return true;
  }
  return false;
 }
 bool run(){
  enumerate(0,Set{},0);
  vector<U> keys;
  inverted.resize(rows.size());
  for(auto& [z,p]:bank){
   int n;Set a=decode(z,n);keyCount[n]++;
   keys.push_back(z);
   if(n<cap)for(int i=0;i<n;i++)inverted[a[i]].push_back(z);
  }
  cout<<"cap "<<cap<<" rows "<<rows.size()<<" candidates";for(int i=0;i<=cap;i++)cout<<' '<<inventory[i];
  cout<<" indexed supports";for(int i=0;i<=cap;i++)cout<<' '<<keyCount[i];cout<<endl;
  sort(keys.begin(),keys.end());
  for(auto& v:inverted)sort(v.begin(),v.end());
  bool yes=false;
  for(U z:keys)if(visit(z)){yes=true;break;}
  cout<<"tested supports";for(int i=0;i<=cap;i++)cout<<' '<<tested[i];cout<<" gain "<<yes<<endl;
  return yes;
 }
};
int main(int argc,char**argv){
 if(argc!=4)return 1;
 int n=stoi(argv[3]);if(n<1||n>5)return 1;
 Search x(argv[1],argv[2],n); x.run();
}
