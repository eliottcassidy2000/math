// Independent vector-based replay: no implicit sequence tree or discovery code.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
uint64_t checks=0;
void check(bool x,const char*m){checks++;if(!x)throw runtime_error(m);}
uint32_t readu(ifstream&f,int bytes){uint32_t v=0;for(int k=0;k<bytes;k++){int c=f.get();check(c!=EOF,"truncated trace");v|=uint32_t(c)<<(8*k);}return v;}
bool square(int v){int r=sqrt(v);while(r*r>v)--r;while((r+1)*(r+1)<=v)++r;return r*r==v;}
bool protect(int x,int y){return(x==3362&&(y==482||y==6242))||(y==3362&&(x==482||x==6242));}
pair<int,int> edge(int x,int y){return minmax(x,y);}
int main(){
 ifstream initial("initial_cycle.txt");vector<int>a;int v;while(initial>>v)a.push_back(v);
 auto validate=[&](){auto b=a;sort(b.begin(),b.end());for(int k=0;k<(int)b.size();k++)check(b[k]==k+1,"coverage");int seed=0;for(int k=0;k<(int)a.size();k++){int x=a[k],y=a[(k+1)%a.size()];check(square(x+y),"square edge");seed+=protect(x,y);}check(seed==2,"seed");};
 validate();ifstream f("full_incremental_trace.bin",ios::binary);char magic[8];f.read(magic,8);check(string(magic,8)==string("SSINC1\0\0",8),"header");int first=readu(f,4),last=readu(f,4);check(first==6242&&int(a.size())==first&&last==128452,"range");set<pair<int,int>>bad;uint64_t moves=0;
 auto remove=[&](int x,int y){check(!protect(x,y),"protected removal");auto e=edge(x,y);if(!square(x+y)){check(bad.count(e),"bad deletion");bad.erase(e);}else check(!bad.count(e),"good deletion");};
 auto add=[&](int x,int y){if(!square(x+y))bad.insert(edge(x,y));};
 for(int N=first+1;N<=last;N++){
  int k=readu(f,4),count=readu(f,2);check(k>=0&&k<N-1,"insertion index");int x=a[k],y=a[(k+1)%(N-1)];remove(x,y);add(x,N);add(N,y);a.insert(a.begin()+k+1,N);check(bad.size()<=1,"one defect insertion");
  for(int t=0;t<count;t++){
   int i=readu(f,4),j=readu(f,4);check(0<=i&&i<j&&j<N,"reversal indices");int x=a[i],y=a[(i+1)%N],z=a[j],w=a[(j+1)%N];remove(x,y);remove(z,w);add(x,z);add(y,w);reverse(a.begin()+i+1,a.begin()+j+1);check(bad.size()<=1,"one defect repair");moves++;
  }
  check(bad.empty(),"every length closes");
 }
 check(f.get()==EOF,"trailing bytes");validate();ifstream final("final_cycle.txt");for(int x:a){int y;check(bool(final>>y)&&x==y,"final agreement");}int extra;check(!(final>>extra),"final trailing token");
 cout<<"{\"status\":\"passed\",\"first\":"<<first<<",\"last\":"<<last<<",\"cycles\":"<<last-first+1<<",\"moves\":"<<moves<<",\"checks\":"<<checks<<",\"representation\":\"independent explicit vector\"}"<<endl;
}
