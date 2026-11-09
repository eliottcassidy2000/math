#include <cstdio>
#include <cstdint>
#include <map>
#include <array>
#include <algorithm>
#include <vector>
#include <cassert>
struct Rec {int rep=0,h=0,reciprocal=0;bool strong=false;long long labelled=0;};
int main(){std::map<uint64_t,Rec> classes; std::vector<std::array<int,6>> perms; std::array<int,6> p={0,1,2,3,4,5};do {perms.push_back(p);}while(std::next_permutation(p.begin(),p.end()));constexpr int n=6,N=64; long long hist[721]={}; long long strong7=0,strong21=0; int first7=-1,first21=-1;for(int code=0;code<14348907;code++) {int z=code,in[n]={},out[n]={};for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){int s=z%3;z/=3;if(s!=1){out[i]|=1<<j;in[j]|=1<<i;}if(s!=0){out[j]|=1<<i;in[i]|=1<<j;}}uint16_t dp[N][n]={};for(int j=0;j<n;j++)dp[1<<j][j]=1;for(int s=1;s<N;s++)for(int js=s;js;js&=js-1){int j=__builtin_ctz((unsigned)js),p=s^(1<<j);if(!p)continue;int v=0;for(int is=p&in[j];is;is&=is-1)v+=dp[p][__builtin_ctz((unsigned)is)];dp[s][j]=v;}int h=0;for(int j=0;j<n;j++)h+=dp[N-1][j];hist[h]++;if(h==7||h==21){int reach=1;for(int k=0;k<n;k++){int r=reach;for(int s=reach;s;s&=s-1)r|=out[__builtin_ctz((unsigned)s)];reach=r;}bool strong=reach==63;reach=1;for(int k=0;k<n;k++){int r=reach;for(int s=reach;s;s&=s-1)r|=in[__builtin_ctz((unsigned)s)];reach=r;}strong&=reach==63;int hd=0;for(auto p:perms){bool ok=true;for(int i=0;i<5;i++)if(!(out[p[i]]>>p[i+1]&1)){ok=false;break;}hd+=ok;}assert(h==hd);
uint64_t canon=~uint64_t(0);for(auto p:perms){uint64_t v=0;for(int i=0;i<6;i++)for(int j=0;j<6;j++)v=(v<<1)|(out[p[i]]>>p[j]&1);if(v<canon)canon=v;}
auto &rec=classes[canon];if(!rec.labelled){rec.rep=code;rec.h=h;rec.strong=strong;for(int i=0;i<6;i++)for(int j=i+1;j<6;j++)rec.reciprocal+=(out[i]>>j&1)&&(out[j]>>i&1);}rec.labelled++;
if(strong){if(h==7){strong7++;if(first7<0)first7=code;}else{strong21++;if(first21<0)first21=code;}}}}
FILE*g=fopen("classes6.csv","w");fprintf(g,"canonical_binary_integer,representative_code,hamiltonian_paths,strong,reciprocal_pairs,labelled_count\n");for(auto [c,r]:classes)fprintf(g,"%llu,%d,%d,%d,%d,%lld\n",(unsigned long long)c,r.rep,r.h,r.strong,r.reciprocal,r.labelled);fclose(g);FILE*f=fopen("results6.csv","w");fprintf(f,"hamiltonian_paths,labelled_count\n");for(int h=0;h<=720;h++)if(hist[h])fprintf(f,"%d,%lld\n",h,hist[h]);fclose(f);printf("h7=%lld strong7=%lld first_strong7=%d h21=%lld strong21=%lld first_strong21=%d\n",hist[7],strong7,first7,hist[21],strong21,first21);}
