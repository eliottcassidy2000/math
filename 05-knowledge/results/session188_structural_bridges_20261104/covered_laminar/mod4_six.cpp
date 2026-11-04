#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
int main(){
 const int n=6,N=1<<(n*(n-1)/2),P=1<<(2*n);int edge[n][n],ec=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edge[i][j]=edge[j][i]=ec++;
 vector<int> basis,idx(P,-1);for(int m=0;m<P;m++)if(__builtin_popcount(m&63)<=2&&__builtin_popcount(m>>6)<=2){idx[m]=basis.size();basis.push_back(m);}int K=basis.size();
 vector<unsigned char> table(K*N,0);array<int,8>w; iota(w.begin(),w.end(),0);long words=0,terms=0;
 do{bool bad=false; vector<pair<int,int>> polys={{0,1}};int pos=0,neg=0;
  for(int k=0;k<7;k++){int u=w[k],v=w[k+1];if(u<n&&v<n){if(u<v)pos|=1<<edge[u][v];else neg|=1<<edge[u][v];}
   else if(u>=n&&v>=n){if(u!=n){bad=true;break;}}
   else{int bit=1<<((u<n)?(v-n)*n+u:(u-n)*n+v);vector<pair<int,int>> pp;
    for(auto [m,z]:polys){if(u<n)pp.emplace_back(m|bit,z);else{pp.emplace_back(m,z);pp.emplace_back(m|bit,-z);}}polys.swap(pp);}}
  if(bad){continue;} words++;
  for(auto [m,z]:polys){int b=idx[m];if(b<0){cerr<<"degree error";return 2;}for(int s=neg;;s=(s-1)&neg){int c=z*((__builtin_popcount(s)&1)?-1:1);auto &v=table[b*N+(pos|s)];v=(v+c+4)&3;terms++;if(s==0)break;}}
 }while(next_permutation(w.begin(),w.end()));
 for(int b=0;b<K;b++)for(int bit=1;bit<N;bit<<=1)for(int m=0;m<N;m++)if(m&bit)table[b*N+m]=(table[b*N+m]+table[b*N+(m^bit)])&3;
 map<string,vector<int>>classes;for(int m=0;m<N;m++){string s;for(int b=0;b<K;b++)s+=(char)table[b*N+m];classes[s].push_back(m);}
 int collisions=0,pairs=0;ofstream out("mod4_six_collisions.json");out<<"[";bool first=true;for(auto &kv:classes)if(kv.second.size()>1){collisions++;pairs+=kv.second.size()*(kv.second.size()-1)/2;if(!first)out<<",";first=false;out<<"[";for(int j=0;j<(int)kv.second.size();j++){if(j)out<<",";out<<kv.second[j];}out<<"]";}out<<"]\n";
 ofstream binary("mod4_six_vectors.bin",ios::binary);binary.write((char*)table.data(),table.size());
 cout<<"{\"n\":6,\"tournaments\":"<<N<<",\"basis_coordinates\":"<<K<<",\"surviving_hamilton_words\":"<<words<<",\"signed_monomial_contributions\":"<<terms<<",\"response_mod4_classes\":"<<classes.size()<<",\"nontrivial_classes\":"<<collisions<<",\"collision_pairs\":"<<pairs<<"}\n";
}
