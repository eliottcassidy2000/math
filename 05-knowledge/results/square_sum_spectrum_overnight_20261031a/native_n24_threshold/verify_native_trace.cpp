// Independent full-vector replay: checks EVERY complete intermediate cycle.
// No search, treap, bad-edge induction, or external graph solver is used.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
uint64_t checks=0;
void need(bool b,const char* m){++checks;if(!b)throw runtime_error(m);}
uint32_t readle(istream& s,int n){uint32_t z=0;for(int i=0;i<n;++i){int c=s.get();need(c!=EOF,"truncated trace");z|=uint32_t(c)<<(8*i);}return z;}
int main(int argc,char**argv){try{
 need(argc==3,"usage: verifier uncompressed-trace uncompressed-final");
 ifstream f(argv[1],ios::binary),initial("initial_cycle.txt"),final(argv[2]);
 char magic[8];f.read(magic,8);const char expected[8]={'S','S','I','N','C','1',0,0};
 need(equal(magic,magic+8,expected),"magic");
 int first=readle(f,4),last=readle(f,4);need(first==128&&last==128452,"range");
 vector<int>a;int z;while(initial>>z)a.push_back(z);need(a.size()==128,"initial length");
 vector<bool>sq(2*last+1,false);for(int r=1;r*r<=2*last;++r)sq[r*r]=true;
 vector<int>seen(last+1,0),pos(last+1,0);
 const int path[]={1,8,17,19,6,10,15,21,4,12,24};
 uint64_t edges=0,moves=0;
 auto validate=[&](int N){
  need(a.size()==size_t(N),"cycle length");
  for(int i=0;i<N;++i){int x=a[i];need(x>=1&&x<=N,"vertex range");need(seen[x]!=N,"duplicate vertex");seen[x]=N;pos[x]=i;}
  for(int i=0;i<N;++i){need(sq[a[i]+a[(i+1)%N]],"nonsquare cycle edge");++edges;}
  for(int i=0;i<10;++i){int d=abs(pos[path[i]]-pos[path[i+1]]);need(d==1||d==N-1,"native path lost");}
 };
 validate(first);
 for(int N=first+1;N<=last;++N){
  int k=readle(f,4),c=readle(f,2);need(k>=0&&k<N-1,"insert index");a.insert(a.begin()+k+1,N);
  for(int j=0;j<c;++j){int l=readle(f,4),r=readle(f,4);need(0<=l&&l<r&&r<N,"reversal range");reverse(a.begin()+l+1,a.begin()+r+1);++moves;}
  validate(N);
  if(N%10000==0)cerr<<"Full cycles verified through "<<N<<'\n';
 }
 need(f.get()==EOF,"trailing trace bytes");
 for(int x:a){need(bool(final>>z),"truncated final");need(x==z,"final disagreement");}
 need(!(final>>z),"trailing final vertices");
 cout<<"{\"status\":\"passed\",\"first\":"<<first<<",\"last\":"<<last
     <<",\"full_cycles\":"<<last-first+1<<",\"square_edges_checked\":"<<edges
     <<",\"checks\":"<<checks<<",\"moves\":"<<moves<<"}\n";
 return 0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
