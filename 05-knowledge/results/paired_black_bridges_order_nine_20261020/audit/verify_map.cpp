#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
constexpr int n=9,m=28;
struct Rec {int u;uint32_t t;uint64_t mass,blue,H,aut;int sc,strong,c3;uint64_t degree;};
array<int,n> arc,path; vector<pair<int,int>> chords; vector<uint32_t> fibre; uint64_t hp;
uint32_t encode(){uint32_t t=0;for(int k=0;k<m;k++){auto [i,j]=chords[k];if(!(arc[path[i]]&(1<<path[j]))) t|=1u<<k;}return t;}
void visit(int used,int depth){if(depth==n){hp++;uint32_t t=encode();fibre.push_back(t);return;}for(int v=0;v<n;v++)if(!(used>>v&1)&&(depth==0||(arc[path[depth-1]]>>v&1))){path[depth]=v;visit(used|(1<<v),depth+1);}}
int main(){for(int j=n-1;j>=2;j--)for(int i=0;i<j-1;i++)chords.push_back({i,j});assert(chords.size()==m);vector<Rec> rows;ifstream records("audit/expected_nodes.tsv");Rec r;while(records>>r.u>>r.t>>r.mass>>r.blue>>r.H>>r.aut>>r.sc>>r.strong>>r.c3>>r.degree)rows.push_back(r);vector<int32_t> mp(1u<<m);ifstream data("census_n9.map",ios::binary);data.read((char*)mp.data(),mp.size()*4);assert(data.gcount()==streamsize(mp.size()*4));vector<uint64_t> mass(rows.size());vector<uint32_t> minimum(rows.size(),UINT32_MAX);for(uint32_t t=0;t<mp.size();t++){int u=mp[t];assert(u>=0&&u<(int)rows.size());mass[u]++;minimum[u]=min(minimum[u],t);}uint64_t totalH=0;int sc=0,strong=0;array<int,m> sig;for(int k=0;k<m;k++){auto[i,j]=chords[k];sig[k]=find(chords.begin(),chords.end(),pair<int,int>{n-1-j,n-1-i})-chords.begin();}auto refl=[&](uint32_t t){uint32_t rt=0;while(t){int k=__builtin_ctz(t);rt|=1u<<sig[k];t&=t-1;}return rt;};
 for(auto&r:rows){assert(mass[r.u]==r.mass&&minimum[r.u]==r.t);arc.fill(0);for(int i=0;i<n-1;i++)arc[i]|=1<<(i+1);for(int k=0;k<m;k++){auto[i,j]=chords[k];if(r.t>>k&1)arc[j]|=1<<i;else arc[i]|=1<<j;}int cyc=84;for(int a:arc){int d=__builtin_popcount((unsigned)a);cyc-=d*(d-1)/2;}assert(cyc==r.c3);array<int,n> reach=arc;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];bool isstrong=true;for(int i=0;i<n;i++)if((reach[i]|(1<<i))!=511)isstrong=false;assert(isstrong==bool(r.strong));strong+=isstrong;
 fibre.clear();hp=0;visit(0,0);assert(hp==r.H);totalH+=hp;sort(fibre.begin(),fibre.end());fibre.erase(unique(fibre.begin(),fibre.end()),fibre.end());uint64_t ordinaryMass=fibre.size();bool self=binary_search(fibre.begin(),fibre.end(),refl(r.t));assert(self==bool(r.sc));sc+=self;assert(hp%ordinaryMass==0&&hp/ordinaryMass==r.aut);for(size_t i=0;i<ordinaryMass;i++)fibre.push_back(refl(fibre[i]));sort(fibre.begin(),fibre.end());fibre.erase(unique(fibre.begin(),fibre.end()),fibre.end());assert(fibre.size()==r.mass);uint64_t blue=0;for(uint32_t t:fibre){assert(mp[t]==r.u);blue+=refl(t)==t;}assert(blue==r.blue);assert(r.degree==(r.mass-r.blue)/2);if(r.u%10000==0)cerr<<"Verified representative "<<r.u<<"\n";
 }cout<<"{\"tilings\":"<<mp.size()<<",\"merged_classes\":"<<rows.size()<<",\"self_converse\":"<<sc<<",\"ordinary_classes\":"<<2*rows.size()-sc<<",\"hamiltonian_paths_reenumerated\":"<<totalH<<",\"strong_merged_classes\":"<<strong<<",\"all_checks\":\"passed\"}\n";
}
