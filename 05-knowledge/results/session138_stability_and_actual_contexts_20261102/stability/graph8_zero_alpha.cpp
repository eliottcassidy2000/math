#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdint>
using namespace std;
vector<pair<int,int>> es;int ei[8][8];vector<uint32_t> M3,M4;
void matching(int vertices,int k,uint32_t mask,vector<uint32_t>&out){if(!k){out.push_back(mask);return;}if(__builtin_popcount((unsigned)vertices)<2*k)return;int i=__builtin_ctz((unsigned)vertices);int tail=vertices^(1<<i);matching(tail,k,mask,out);for(int j=i+1;j<8;j++)if(tail>>j&1)matching(tail^(1<<j),k-1,mask|(1u<<ei[i][j]),out);}
int eval(uint32_t g){int d[8]={},m=__builtin_popcount(g);for(int e=0;e<28;e++)if(g>>e&1){d[es[e].first]++;d[es[e].second]++;}int w=0;for(int x:d)w+=x*(x-1)/2;int m2=m*(m-1)/2-w,m3=0,m4=0,p21=0;for(auto a:M3)m3+=(g&a)==a;for(auto a:M4)m4+=(g&a)==a;
for(int e=0;e<28;e++)if(g>>e&1){auto [i,j]=es[e];for(int v=0;v<8;v++)if(v!=i&&v!=j){int q=d[v]-((g>>ei[v][i])&1)-((g>>ei[v][j])&1);p21+=q*(q-1)/2;}}
return 41*m-32*m2+8*m3-4*m4+8*p21;}
int main(){for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){ei[i][j]=ei[j][i]=es.size();es.push_back({i,j});}matching(255,3,0,M3);matching(255,4,0,M4);
vector<pair<int,int>> e7;for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)e7.push_back({i,j});int idx7[7][7];for(int e=0;e<21;e++){auto[i,j]=e7[e];idx7[i][j]=idx7[j][i]=e;}
vector<array<uint32_t,21>> perms;array<int,7>p={0,1,2,3,4,5,6};do{array<uint32_t,21>a;for(int e=0;e<21;e++){auto[i,j]=e7[e];a[e]=1u<<idx7[p[i]][p[j]];}perms.push_back(a);}while(next_permutation(p.begin(),p.end()));
vector<uint8_t> seen(1u<<21);int reps=0,covered=0,minimum=100000,zeros=0;uint64_t cases=0;uint32_t witness=0;
for(uint32_t g=0;g<(1u<<21);g++)if(!seen[g]){reps++;for(auto &a:perms){uint32_t image=0;for(uint32_t b=g;b;b&=b-1)image|=a[__builtin_ctz(b)];if(!seen[image]){seen[image]=1;covered++;}}
uint32_t base=0;for(int e=0;e<21;e++)if(g>>e&1)base|=1u<<ei[e7[e].first][e7[e].second];for(int attach=0;attach<128;attach++){uint32_t h=base;for(int i=0;i<7;i++)if(attach>>i&1)h|=1u<<ei[i][7];int z=eval(h);if(z<minimum){minimum=z;witness=h;}zeros+=z==0;cases++;}}
cout<<"{\"n\":8,\"seven_representatives\":"<<reps<<",\"seven_covered\":"<<covered<<",\"augmentations\":"<<cases<<",\"matching3_patterns\":"<<M3.size()<<",\"matching4_patterns\":"<<M4.size()<<",\"minimum\":"<<minimum<<",\"zero_augmentations\":"<<zeros<<",\"witness\":"<<witness<<"}\n";
}
