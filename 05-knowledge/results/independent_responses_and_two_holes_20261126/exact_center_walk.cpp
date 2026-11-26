
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
constexpr int N=52;
struct Op{int y;array<int,3>s;};
uint64_t state;
uint64_t rnd(){state^=state>>12;state^=state<<25;state^=state>>27;return state*2685821657736338717ull;}
vector<Op> op[N];int reference[N], pick[N],sorter[N], own[56], use[3][56], k;
bool mark[N],slot[56];
bool augment(int i){
 for(int t=0,d=rnd()%op[i].size();t<(int)op[i].size();t++){
  int j=(t+d)%op[i].size(),z=op[i][j].s[k];
  if(slot[z])continue;slot[z]=true;
  if(own[z]<0||augment(own[z])){own[z]=i;pick[i]=j;return true;}
 }
 return false;
}
int add(int i,int e,int d){
 int a=0;
 for(int j=0;j<3;j++)if(j!=k){
  int &z=use[j][op[i][e].s[j]];
  a+=d>0?z:1-z;z+=d;
 }
 return a;
}
int main(int argc,char**argv){
 if(argc!=4)return 2;
 k=stoi(argv[1]);state=stoull(argv[2]);uint64_t maxit=stoull(argv[3]);
 int n;cin>>n;assert(n==53);
 int marker;array<int,3> ms;
 cin>>marker>>n>>marker>>ms[0]>>ms[1]>>ms[2];
 for(int i=0;i<N;i++){
  int q;cin>>reference[i]>>q;sorter[i]=i;
  while(q--){
   Op x;cin>>x.y>>x.s[0]>>x.s[1]>>x.s[2];
   if(x.s[0]!=ms[0]&&x.s[1]!=ms[1]&&x.s[2]!=ms[2])op[i].push_back(x);
  }
 }
 int best=1000;
 const int block=500000;
 for(uint64_t epoch=0;epoch<=maxit/block;epoch++){
  for(int j=N-1;j>0;j--)swap(sorter[j],sorter[rnd()%(j+1)]);
  fill(own,own+56,-1);
  bool seeded=(epoch%3!=2&&k==2);
  if(seeded){
   int holes=0;
   for(int i=0;i<N;i++){
    const int y=reference[i]>=0?reference[i]:(holes++==0?16:19);
    int j=0;while(j<(int)op[i].size()&&op[i][j].y!=y)j++;
    if(j==(int)op[i].size()){seeded=false;break;}
    const int s=op[i][j].s[k];if(own[s]>=0){seeded=false;break;}own[s]=i;pick[i]=j;
   }
   assert(seeded);
  }else for(int i:sorter){
   fill(slot,slot+56,false);
   if(!augment(i)){cout<<"projection matching absent"<<endl;return 0;}
  }
  for(auto &u:use)fill(u,u+56,0);
  int score=0;for(int i=0;i<N;i++)score+=add(i,pick[i],1);
  best=min(best,score);
  for(int it=0;it<block&&epoch*block+it<maxit;it++){
   const double t=0.92-0.72*(double(it)/block);
   int X[10],E[10], len=0,u=rnd()%N;uint64_t seen=0;
   bool ok=false;
   int lim=2+rnd()%((it%7==0)?8:4);
   while(len<lim){
    X[len]=u;seen|=1ull<<u;
    int e=rnd()%op[u].size();
    if(len==lim-1){
     int e2=-1,start=rnd()%op[u].size();
     for(int z=0;z<(int)op[u].size();z++){
      int j=(z+start)%op[u].size();
      int v=own[op[u][j].s[k]];
      if(v<0||v==X[0]){e2=j;break;}
     }
     if(e2<0)break;e=e2;
    }
    if(e==pick[u])break;
    int v=own[op[u][e].s[k]];
    E[len++]=e;
    if(v<0||v==X[0]){ok=true;break;}
    if(seen>>v&1)break;
    u=v;
   }
   if(!ok)continue;
   int diff=0;
   for(int j=0;j<len;j++)diff+=add(X[j],pick[X[j]],-1);
   for(int j=0;j<len;j++)diff+=add(X[j],E[j],1);
   const bool yes=diff<=0 || ((rnd()>>11)*0x1.0p-53 < exp(-diff/t));
   if(yes){
    score+=diff;
    for(int j=0;j<len;j++)own[op[X[j]][pick[X[j]]].s[k]]=-1;
    for(int j=0;j<len;j++){
     pick[X[j]]=E[j];int z=op[X[j]][E[j]].s[k];assert(own[z]<0);own[z]=X[j];
    }
    if(score<best){
     best=score;
     cerr<<"center "<<k<<" score "<<best<<" at "<<epoch*block+it<<endl;
     if(best<=2){
      array<int,3> q{};int edits=0;
      for(int j=0;j<3;j++)if(j!=k)for(int s:use[j])q[j]+=s*(s-1)/2;
      for(int i=0;i<N;i++)edits+=reference[i]>=0&&op[i][pick[i]].y!=reference[i];
      assert(q[0]+q[1]+q[2]==best);
      cout<<"center "<<k<<" score "<<best<<" profile "<<q[0]<<" "<<q[1]<<" "<<q[2]<<" edits "<<edits<<endl;
      for(int i=0;i<N;i++)cout<<op[i][pick[i]].y<<" ";cout<<endl;
      if(!best)return 0;
     }
    }
   }else{
    for(int j=0;j<len;j++)add(X[j],E[j],-1);
    for(int j=0;j<len;j++)add(X[j],pick[X[j]],1);
   }
  }
 }
 cerr<<"complete center "<<k<<" budget "<<maxit<<" best "<<best<<endl;
}
