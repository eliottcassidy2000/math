#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
using namespace std;
struct Perm { vector<int> ep; unsigned rev; };
int main(){
 cout<<"[\n";
 for(int n=0;n<=7;n++){
  vector<pair<int,int>> ed; int ids[7][7]={};
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ids[i][j]=ed.size();ed.emplace_back(i,j);}
  int N=ed.size(); vector<int> p(n);iota(p.begin(),p.end(),0);vector<Perm> pp;
  do{Perm pm;pm.rev=0;for(int e=0;e<N;e++){auto [i,j]=ed[e];int a=p[i],b=p[j];if(a>b){pm.rev|=1u<<e;swap(a,b);}pm.ep.push_back(ids[a][b]);}pp.push_back(pm);}while(next_permutation(p.begin(),p.end()));
  vector<bool> seen(1u<<N);int ext=0,roy=0,graphs=0,b=0,f=0,rc=0;vector<int> exdeg(N+1),rodeg(N+1),exconn(N+1);vector<unsigned> goodconn;
  for(unsigned m=0;m<(1u<<N);m++)if(!seen[m]){
   graphs++;bool eg=true,rg=true;
   for(auto&pm:pp){unsigned im=0;vector<int> images;
    for(int e=0;e<N;e++)if(m>>e&1){im|=1u<<pm.ep[e];images.push_back(pm.ep[e]);}
    seen[im]=true;
    if(im==m){int r=__builtin_popcount(m&pm.rev)%2,ep=0;for(int i=0;i<(int)images.size();i++)for(int j=i+1;j<(int)images.size();j++)ep^=(images[i]>images[j]);if(r)rg=false;if(r^ep)eg=false;}
   }
   int k=__builtin_popcount(m);if(eg){ext++;exdeg[k]++;}if(rg){roy++;rodeg[k]++;}
   unsigned reach=n?1:0,prev;
   do{prev=reach;for(int e=0;e<N;e++)if(m>>e&1){auto [i,j]=ed[e];if(reach>>i&1)reach|=1u<<j;if(reach>>j&1)reach|=1u<<i;}}while(prev!=reach);
   bool conn=n&&reach==(1u<<n)-1;
   if(conn&&eg){exconn[k]++;if(k%2)f++;else b++;goodconn.push_back(m);}if(conn&&rg)rc++;
  }
  cout<<" {\"n\":"<<n<<",\"graphs\":"<<graphs<<",\"exterior_types\":"<<ext<<",\"royle_types\":"<<roy<<",\"bosonic_connected\":"<<b<<",\"fermionic_connected\":"<<f<<",\"royle_connected\":"<<rc<<",\"exterior_degrees\":[";
  for(int k=0;k<=N;k++){if(k)cout<<',';cout<<exdeg[k];}cout<<"],\"royle_degrees\":[";for(int k=0;k<=N;k++){if(k)cout<<',';cout<<rodeg[k];}cout<<"],\"connected_exterior_degrees\":[";for(int k=0;k<=N;k++){if(k)cout<<',';cout<<exconn[k];}cout<<"]}"<<(n<7?",":"")<<'\n';
 }
 cout<<"]\n";
}
