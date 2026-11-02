#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_set>
#include <cstdint>
using namespace std;
int e8[8][8],e9[9][9];vector<pair<int,int>> es8,es9;vector<uint64_t> pat3,pat4;
void match(int vs,int k,uint64_t g,vector<uint64_t>&out){if(!k){out.push_back(g);return;}if(__builtin_popcount((unsigned)vs)<2*k)return;int i=__builtin_ctz((unsigned)vs),tail=vs^(1<<i);match(tail,k,g,out);for(int j=i+1;j<9;j++)if(tail>>j&1)match(tail^(1<<j),k-1,g|(1ull<<e9[i][j]),out);}
uint32_t canon8(uint32_t g){int d[8]={};for(int i=0;i<28;i++)if(g>>i&1){d[es8[i].first]++;d[es8[i].second]++;}array<int,8>p;for(int i=0;i<8;i++)p[i]=i;sort(p.begin(),p.end(),[&](int a,int b){return d[a]!=d[b]?d[a]<d[b]:a<b;});vector<pair<int,int>>blocks;for(int i=0;i<8;){int j=i+1;while(j<8&&d[p[j]]==d[p[i]])j++;blocks.push_back({i,j});i=j;}uint32_t best=~0u;
 auto rec=[&](auto&&self,int k)->void{if(k==(int)blocks.size()){uint32_t h=0;for(int i=0;i<28;i++){auto[a,b]=es8[i];h|=((g>>e8[p[a]][p[b]])&1)<<i;}best=min(best,h);return;}auto[a,b]=blocks[k];do{self(self,k+1);}while(next_permutation(p.begin()+a,p.begin()+b));};rec(rec,0);return best;}
array<int,3> eval9(uint64_t g){int d[9]={},m=__builtin_popcountll(g);for(int e=0;e<36;e++)if(g>>e&1){d[es9[e].first]++;d[es9[e].second]++;}int W=0;for(int z:d)W+=z*(z-1)/2;int m2=m*(m-1)/2-W,m3=0,m4=0,p21=0;for(uint64_t p:pat3)m3+=(g&p)==p;for(uint64_t p:pat4)m4+=(g&p)==p;
 for(int e=0;e<36;e++)if(g>>e&1){auto[i,j]=es9[e];for(int v=0;v<9;v++)if(v!=i&&v!=j){int z=d[v]-((g>>e9[v][i])&1)-((g>>e9[v][j])&1);p21+=z*(z-1)/2;}}
 int base=81*m-64*m2+16*m3-4*m4+16*p21;return{base,W,m};}
int main(){for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){e8[i][j]=e8[j][i]=es8.size();es8.push_back({i,j});}for(int i=0;i<9;i++)for(int j=i+1;j<9;j++){e9[i][j]=e9[j][i]=es9.size();es9.push_back({i,j});}match(511,3,0,pat3);match(511,4,0,pat4);
 vector<pair<int,int>> es7;int e7[7][7];for(int i=0;i<7;i++)for(int j=i+1;j<7;j++){e7[i][j]=e7[j][i]=es7.size();es7.push_back({i,j});}vector<array<uint32_t,21>>perms;array<int,7>p={0,1,2,3,4,5,6};do{array<uint32_t,21>a;for(int e=0;e<21;e++){auto[i,j]=es7[e];a[e]=1u<<e7[p[i]][p[j]];}perms.push_back(a);}while(next_permutation(p.begin(),p.end()));
 vector<uint8_t>seen(1u<<21);unordered_set<uint32_t> reps8;int reps7=0,covered=0;uint64_t aug8=0;
 for(uint32_t g=0;g<(1u<<21);g++)if(!seen[g]){reps7++;for(auto&a:perms){uint32_t h=0;for(uint32_t b=g;b;b&=b-1)h|=a[__builtin_ctz(b)];if(!seen[h]){seen[h]=1;covered++;}}uint32_t base=0;for(int e=0;e<21;e++)if(g>>e&1)base|=1u<<e8[es7[e].first][es7[e].second];for(int att=0;att<128;att++){uint32_t h=base;for(int i=0;i<7;i++)if(att>>i&1)h|=1u<<e8[i][7];reps8.insert(canon8(h));aug8++;}}
 vector<uint32_t> sorted8(reps8.begin(),reps8.end());sort(sorted8.begin(),sorted8.end());cerr<<"reps8 "<<sorted8.size()<<"\n";
 int min19=100000,min32=100000,z19=0,z32=0,an=0,ad=1;uint64_t cases=0,witness=0;
 for(uint32_t g:sorted8){uint64_t base=0;for(int e=0;e<28;e++)if(g>>e&1)base|=1ull<<e9[es8[e].first][es8[e].second];for(int att=0;att<256;att++){uint64_t h=base;for(int i=0;i<8;i++)if(att>>i&1)h|=1ull<<e9[i][8];auto a=eval9(h);int D19=a[0]+19*a[1],D32=a[0]+32*a[1];if(D19<min19){min19=D19;witness=h;}min32=min(min32,D32);z19+=D19==0;z32+=D32==0;if(a[1]&&(int64_t)(-a[0])*ad>(int64_t)an*64*a[1]){an=-a[0];ad=64*a[1];}cases++;}}
 cout<<"{\"n\":9,\"reps7\":"<<reps7<<",\"covered7\":"<<covered<<",\"augmentations8\":"<<aug8<<",\"reps8\":"<<sorted8.size()<<",\"augmentations9\":"<<cases<<",\"m3_patterns\":"<<pat3.size()<<",\"m4_patterns\":"<<pat4.size()<<",\"minimum_alpha19over64\":"<<min19<<",\"zeros_alpha19over64\":"<<z19<<",\"minimum_alpha1over2\":"<<min32<<",\"zeros_alpha1over2\":"<<z32<<",\"max_required_alpha_numerator\":"<<an<<",\"max_required_alpha_denominator\":"<<ad<<",\"witness\":"<<witness<<"}\n";
}
