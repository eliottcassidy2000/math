#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;
int ir(int v){int r=sqrt(v);while((r+1)*(r+1)<=v)r++;while(r*r>v)r--;return r;}
bool protect(int x,int y){return(x==3362&&(y==482||y==6242))||(y==3362&&(x==482||x==6242));}
void u32(ofstream&o,uint32_t x){for(int k=0;k<4;k++)o.put((x>>(8*k))&255);}
void u16(ofstream&o,uint16_t x){o.put(x&255);o.put(x>>8);}
int main(int argc,char**argv){
 int last=argc>1?stoi(argv[1]):128452,maxmoves=2000;double timebudget=1800;
 ifstream in("initial_cycle.txt");vector<int>a;int z;while(in>>z)a.push_back(z);if(a.size()!=6242)return 2;
 vector<int>pos(last+1,-1);for(int k=0;k<(int)a.size();k++)pos[a[k]]=k;
 vector<bool>sq(2*last+1,false);for(int r=1;r*r<=2*last;r++)sq[r*r]=true;
 mt19937 rng(20261019);auto start=chrono::steady_clock::now();
 ofstream trace("full_incremental_trace.bin",ios::binary);trace.write("SSINC1\0\0",8);u32(trace,6242);u32(trace,last);
 ofstream stats("full_step_stats.tsv");stats<<"N\tmoves\tdirect\n";
 uint64_t moves_total=0,direct=0;int worst=0,completed=6242;
 auto secs=[&](){return chrono::duration<double>(chrono::steady_clock::now()-start).count();};
 for(int N=6243;N<=last;N++){
  vector<int>options;int best=-1;
  for(int r=ir(N)+1;r<=ir(2*N-1);r++){int v=r*r-N;if(v<1||v>=N)continue;int k=pos[v],w=a[(k+1)%(N-1)];if(protect(v,w))continue;int q=sq[N+w];if(q>best){best=q;options.clear();}if(q==best)options.push_back(k);}
  if(options.empty()){cerr<<"No insertion at "<<N<<endl;break;}
  int insertion=options[rng()%options.size()],w=a[(insertion+1)%(N-1)];a.insert(a.begin()+insertion+1,N);for(int k=insertion+1;k<N;k++)pos[a[k]]=k;
  bool bad=!sq[N+w];int x=N,y=w;vector<pair<int,int>>moves;
  while(bad&&(int)moves.size()<maxmoves&&secs()<timebudget){
   int i=pos[x];if(a[(i+1)%N]!=y){swap(x,y);i=pos[x];}
   vector<int>candidates;
   for(int r=ir(x)+1;r<=ir(x+N);r++){int c=r*r-x;if(c>=1&&c<=N&&c!=x)candidates.push_back(pos[c]);}
   for(int r=ir(y)+1;r<=ir(y+N);r++){int c=r*r-y;if(c>=1&&c<=N&&c!=y)candidates.push_back((pos[c]+N-1)%N);}
   sort(candidates.begin(),candidates.end());candidates.erase(unique(candidates.begin(),candidates.end()),candidates.end());
   options.clear();best=-1;
   for(int j:candidates){if(j==i||j==(i+1)%N||j==(i+N-1)%N)continue;int c=a[j],d=a[(j+1)%N];if(protect(c,d))continue;int gain=int(sq[x+c])+int(sq[y+d])-1;if(gain<0)continue;if(gain>best){best=gain;options.clear();}if(gain==best)options.push_back(j);}
   if(options.empty())break;
   int j=options[rng()%options.size()],c=a[j],d=a[(j+1)%N];bool good1=sq[x+c],good2=sq[y+d];int l=min(i,j),r=max(i,j);moves.emplace_back(l,r);reverse(a.begin()+l+1,a.begin()+r+1);for(int k=l+1;k<=r;k++)pos[a[k]]=k;
   bad=!(good1&&good2);if(!good1)y=c;else if(!good2){x=y;y=d;}
  }
  if(bad){cerr<<"Incomplete at "<<N<<" moves "<<moves.size()<<" elapsed "<<secs()<<endl;break;}
  u32(trace,insertion);u16(trace,moves.size());for(auto[l,r]:moves){u32(trace,l);u32(trace,r);}
  stats<<N<<'\t'<<moves.size()<<'\t'<<(moves.empty()?1:0)<<'\n';moves_total+=moves.size();direct+=moves.empty();worst=max(worst,(int)moves.size());completed=N;
  if(N%1000==0){trace.flush();stats.flush();cout<<"{\"last\":"<<N<<",\"moves\":"<<moves_total<<",\"worst\":"<<worst<<",\"seconds\":"<<secs()<<"}"<<endl;}
  if(N%8192==0){ofstream checkpoint("checkpoint_"+to_string(N)+".txt");for(int v:a)checkpoint<<v<<'\n';ofstream meta("checkpoint_"+to_string(N)+".json");meta<<"{\"N\":"<<N<<",\"trace_offset\":"<<trace.tellp()<<"}";}
 }
 trace.close();stats.close();
 if(completed==last){ofstream out("final_cycle.txt");for(int v:a)out<<v<<'\n';}
 cout<<"{\"status\":\""<<(completed==last?"complete":"bounded search incomplete")<<"\",\"first\":6242,\"last\":"<<completed<<",\"requested_last\":"<<last<<",\"moves\":"<<moves_total<<",\"direct_insertions\":"<<direct<<",\"max_moves\":"<<worst<<",\"seconds\":"<<secs()<<"}"<<endl;
 return completed==last?0:1;
}
