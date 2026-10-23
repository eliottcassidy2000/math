#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
vector<vector<int>> neighbors;int N;uint64_t oriented, visits;
void extend(int v,uint64_t mask,int length){
 ++visits;
 if(length==N){++oriented;return;}
 for(int w:neighbors[v])if(!(mask&(uint64_t(1)<<w)))extend(w,mask|(uint64_t(1)<<w),length+1);
}
int main(){
 cout<<"[\n";bool first=true;
 for(N=33;N<=36;++N){
  neighbors.assign(N,{});
  for(int r=2;r*r<=2*N;++r)for(int a=1;a<=N;++a){int b=r*r-a;if(a<b&&b<=N){neighbors[a-1].push_back(b-1);neighbors[b-1].push_back(a-1);}}
  oriented=visits=0;for(int a=0;a<N;++a)extend(a,uint64_t(1)<<a,1);
  if(!first)cout<<",\n";first=false;
  cout<<"{\"N\":"<<N<<",\"oriented\":"<<oriented<<",\"unoriented\":"<<oriented/2<<",\"visits\":"<<visits<<"}"<<flush;
 }
 cout<<"\n]\n";
}
