#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
using Word=array<int,7>;
long long gates=0,raw=0,accepted=0;
void need(bool ok,const string&why){++gates;if(!ok){cerr<<why<<'\n';exit(1);}}
array<unordered_set<unsigned long long>,7> profiles;
int gd[91][91];vector<int> masks;
bool valid(const Word&w){
 int whole=0;for(int a:w)whole=gd[whole][a];if(whole!=1)return false;
 int cache[128]={};
 auto divisor=[&](auto&&self,int m)->int{
  if(cache[m])return cache[m];int i=__builtin_ctz((unsigned)m),r=m&(m-1);
  return cache[m]=r?gd[w[i]][self(self,r)]:w[i];
 };
 for(int m:masks){
  int c=divisor(divisor,m),k=7-__builtin_popcount((unsigned)m),j=0;array<int,6>a{};
  for(int i=0;i<7;i++)if(!(m&(1<<i)))a[j++]=gd[c][w[i]];
  sort(a.begin(),a.begin()+k);unsigned long long key=c;
  for(int i=0;i<k;i++)key=(key<<7)|a[i];
  if(!profiles[k].count(key))return false;
 }return true;
}
void visit(const vector<int>&D,Word&w,int n,int first,vector<Word>&out){
 if(n==7){++raw;if(valid(w)){out.push_back(w);++accepted;}return;}
 for(int j=first;j<(int)D.size();j++){w[n]=D[j];visit(D,w,n+1,j,out);}
}
int main(int argc,char**argv){
 need(argc==3,"input and private output directory required");ifstream in(argv[1]);
 for(int i=0;i<=90;i++)for(int j=0;j<=90;j++)gd[i][j]=gcd(i,j);
 int n;in>>n;for(int i=0;i<n;i++){
  int k,c;in>>k>>c;need(1<=k&&k<=6&&1<=c&&c<=90,"profile scalar typing");
  unsigned long long key=c;for(int j=0;j<k;j++){int x;in>>x;need(1<=x&&x<=90,"profile margin typing");key=(key<<7)|x;}
  profiles[k].insert(key);
 }
 for(int k=6;k>=1;k--)for(int m=1;m<127;m++)if(__builtin_popcount((unsigned)m)==k)masks.push_back(m);
 need(masks.size()==126,"all proper nonempty subsets, no partial-word pruning");
 in>>n;ofstream counts(string(argv[2])+"/counts.txt");
 for(int i=0;i<n;i++){
  int id,size;in>>id>>size;vector<int>D(size);for(int&d:D){in>>d;need(d>=1&&d<=90,"complete alphabet bounded");}
  need(is_sorted(D.begin(),D.end()),"alphabet sorted");
  Word w{};vector<Word> words;long long before=raw;visit(D,w,0,0,words);
  counts<<id<<' '<<raw-before<<' '<<words.size()<<'\n';
  ofstream out(string(argv[2])+"/words_"+to_string(id)+".json");out<<'[';
  for(size_t j=0;j<words.size();j++){if(j)out<<',';out<<'[';for(int k=0;k<7;k++){if(k)out<<',';out<<words[j][k];}out<<']';}out<<']';
 }
 need(bool(in),"complete input consumed without a parse failure");
 cout<<"DOMAINS "<<n<<" RAW_MULTISETS "<<raw<<" ACCEPTED_WORDS "<<accepted<<" GATES "<<gates<<'\n';
}
