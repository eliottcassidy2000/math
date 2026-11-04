#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;
struct Term{int require, forbid,pos,neg;};
bool pairmodule(int n,int t){
 int bit[8][8]={}; int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)bit[i][j]=bit[j][i]=k++;
 auto edge=[&](int i,int j){return i<j? bool(t>>bit[i][j]&1):!bool(t>>bit[i][j]&1);};
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){bool yes=true;for(int w=0;w<n;w++)if(w!=i&&w!=j&&edge(i,w)!=edge(j,w))yes=false;if(yes)return true;}
 return false;
}
int main(){
 for(int n=2;n<=6;n++){
 int arc[8][8]={},k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)arc[i][j]=arc[j][i]=k++;
 vector<int> idx(1<<(2*n),-1),mons;for(int m=0;m<(1<<(2*n));m++)if(__builtin_popcount(unsigned(m&((1<<n)-1)))<=2&&__builtin_popcount(unsigned(m>>n))<=2){idx[m]=mons.size();mons.push_back(m);}
 vector<Term> terms;vector<int> p(n+2);iota(p.begin(),p.end(),0);
 do {Term t{0,0,0,0};bool legal=true;for(int j=1;j<n+2;j++){
 int u=p[j-1],v=p[j];
 if(u<n&&v<n){if(u<v)t.require|=1<<arc[u][v];else t.forbid|=1<<arc[u][v];}
 else if(u>=n&&v>=n){if(u>v){legal=false;break;}}
 else {int port=u>=n?u:v, old=u<n?u:v;int b=1<<(old+(port-n)*n);if(u<n)t.pos|=b;else t.neg|=b;}
 }if(legal)terms.push_back(t);
 }while(next_permutation(p.begin(),p.end()));
 map<vector<int>,array<int,3>> seen;long long equal=0,opposite=0;int paircount=0;
 for(int t=0;t<(1<<k);t++){
 vector<int> c(mons.size());for(auto &z:terms)if((t&z.require)==z.require&&!(t&z.forbid)){
 int sub=z.neg;while(true){c[idx[z.pos|sub]]+=(__builtin_popcount(unsigned(sub))%2)?-1:1;if(!sub)break;sub=(sub-1)&z.neg;}}
 bool mod=pairmodule(n,t);paircount+=mod;
 auto it=seen.find(c);if(it!=seen.end()){
 equal++;if(bool(it->second[1])!=mod){opposite++;cout<<"COLLISION n="<<n<<" labels="<<it->second[0]<<","<<t<<endl;}
 it->second[2]++;
 }else seen.emplace(move(c),array<int,3>{t,int(mod),1});
 }
 cout<<"n="<<n<<" tournaments="<<(1<<k)<<" coefficients="<<mons.size()<<" distinct_full_responses="<<seen.size()<<" repeat_responses="<<equal<<" pair_module_tournaments="<<paircount<<" cross_structure_collisions="<<opposite<<endl;
 }
}
