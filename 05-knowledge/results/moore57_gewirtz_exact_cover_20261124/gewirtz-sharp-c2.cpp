#include <bitset>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
const int N=10152,W=56;
using B=bitset<N>;
static unsigned char a[N][W];
static B adj[N];

long long nodes=0;
void rec(B P,int d) {
 ++nodes;
 if(P.count()+d<27)return;
 if(d==27){cout<<"found"<<endl;exit(1);}
 vector<int> ord, bnd;
 for(int i=0;P.any();i++){
  B Q=P;
  while(Q.any()){
    int v=Q._Find_first();P.reset(v);Q.reset(v);Q&=~adj[v];
    ord.push_back(v);bnd.push_back(i+1);
  }
 }
 for(auto v:ord)P.set(v);
 for(int i=ord.size()-1;i>=0;--i){
   if(d+bnd[i]<27)return;
   int v=ord[i];P.reset(v);rec(P&adj[v],d+1);
 }
}
int main() {
 ifstream f("invol.dat",ios::binary);f.read((char*)a,sizeof(a)); if(!f)return 2;
 long long edges=0;
 for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){
  int k=0;for(;k<W&&a[i][k]!=a[j][k];k++);
  if(k==W){adj[i].set(j);adj[j].set(i);++edges;}
 }
 if(edges!=12533328)return 3;
 ifstream r("reps.dat");int i, nrep=0;
 while(r>>i){if(i<0||i>=N)return 4;rec(adj[i],1);++nrep;}
 if(nrep!=41||nodes!=90184)return 5;
 cout<<"No C2-stable sharply transitive set; "<<nrep<<" orbits and "<<nodes<<" search nodes.\n";
}
