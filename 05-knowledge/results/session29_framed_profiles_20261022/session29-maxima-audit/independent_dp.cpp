#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>
using namespace std;
long long paths(int m,unsigned q,unsigned neg){
 unsigned out[16]={};int e=0;
 for(int i=0;i<m;i++)for(int j=i+1;j<m;j++,e++){
  int c=(q>>e&1)?j:i,w=c==i?j:i;
  int minus=(neg>>e)&1;
  for(int sc=0;sc<2;sc++)for(int sw=0;sw<2;sw++){
   int a=2*c+sc,b=2*w+sw;
   // sc=0 is positive; controller goes first iff sign=-b.
   if(sc==1-minus)out[a]|=1u<<b;else out[b]|=1u<<a;
  }
 }
 long long d[256][16]={};for(int v=0;v<2*m;v++)d[1<<(v/2)][v]=1;
 for(int s=1;s<(1<<m);s++)for(int v=0;v<2*m;v++)if(d[s][v]){
  unsigned dest=out[v];for(int j=0;j<m;j++)if(s>>j&1)dest&=~(3u<<(2*j));
  while(dest){int w=__builtin_ctz(dest);dest&=dest-1;d[s|(1<<(w/2))][w]+=d[s][v];}
 }
 long long f=0;for(int v=0;v<2*m;v+=2)f+=d[(1<<m)-1][v];return f;
}
int main(int argc,char**argv){
 int m=atoi(argv[1]);ifstream in(argv[2]);string line;vector<unsigned> reps;vector<long long> expected;vector<unsigned> literal;
 while(getline(in,line)){istringstream s(line);int n;unsigned q;long long f;long long A;unsigned neg;s>>n>>q>>f>>A>>neg;reps.push_back(q);expected.push_back(f);literal.push_back(neg);}
 int E=m*(m-1)/2,ids[8][8];vector<pair<int,int>> es;for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){ids[i][j]=ids[j][i]=es.size();es.push_back({i,j});}
 vector<int> owner(1<<E,-1);unsigned total=0;unsigned long long configs=0;long long global=0;
 for(int ix=0;ix<(int)reps.size();ix++){
  unsigned q=reps[ix];if(paths(m,q,literal[ix])!=expected[ix]){cerr<<"LITERAL MISMATCH";return 5;}vector<int> p(m);iota(p.begin(),p.end(),0);int aut=0,orb=0;
  do {unsigned r=0;for(int k=0;k<E;k++){
   auto [i,j]=es[k];int src=(q>>k&1)?j:i,dst=(q>>k&1)?i:j;
   if(p[src]>p[dst])r|=1u<<ids[p[src]][p[dst]];
  }if(r==q)aut++;if(owner[r]!=-1 && owner[r]!=ix){cerr<<"OVERLAP";return 1;}if(owner[r]==-1){owner[r]=ix;orb++;}}
  while(next_permutation(p.begin(),p.end()));total+=orb;
  vector<unsigned> stars(m);vector<int> free;for(int k=0;k<E;k++){auto [i,j]=es[k];stars[(q>>k&1)?j:i]|=1u<<k;}
  for(auto star:stars)if(star){unsigned x=star&(star-1);while(x){free.push_back(__builtin_ctz(x));x&=x-1;}}
  configs+=1ull<<free.size();long long best=-1;unsigned witness=0;
  if(m==6){for(unsigned z=0;z<(1u<<free.size());z++){
   unsigned neg=0;for(int j=0;j<(int)free.size();j++)if(z>>j&1)neg|=1u<<free[j];
   long long f=paths(m,q,neg),opt=f;unsigned final=neg;
   for(auto star:stars)if(star){long long flipped=paths(m,q,neg^star);if((f-flipped)%2){cerr<<"ODD PHASE";return 2;}if(flipped>f){opt+=flipped-f;final^=star;}}
   if(opt>best){best=opt;witness=final;}
  }if(best!=expected[ix]||paths(m,q,witness)!=best){cerr<<"MISMATCH q="<<q;return 3;}global=max(global,best);}
  cout<<m<<" "<<q<<" orbit="<<orb<<" aut="<<aut<<" configs="<<(1ull<<free.size());if(m==6)cout<<" independent_max="<<best<<" witness="<<witness;cout<<"\n";
 }
 if(total!=(1u<<E)){cerr<<"INCOMPLETE";return 4;}
 cout<<"SUMMARY classes="<<reps.size()<<" covered="<<total<<" configs="<<configs;if(m==6)cout<<" independent_max="<<global;else cout<<" independent_max=not_run";cout<<"\n";
}
