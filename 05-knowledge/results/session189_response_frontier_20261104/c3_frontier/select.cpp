#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;using I=long long;
I E[512][9],S[512][9];
void endpoints(const vector<int>&r,I d[512][9]){int n=r.size(),N=1<<n;fill(&d[0][0],&d[0][0]+512*9,0);for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(m>>v&1){int q=m^(1<<v);I z=q?0:1;for(int u=0;u<n;u++)if(q>>u&1&&r[u]>>v&1)z+=d[q][u];d[m][v]=z;}}
I skew(int i,int j,int F){I z=0;int base=F^(1<<i)^(1<<j);for(int sub=base;;sub=(sub-1)&base){int a=sub|(1<<i),b=F^a;z+=E[a][i]*S[b][j]-S[a][i]*E[b][j];if(!sub)break;}return z;}
bool module(const vector<int>&r,int m){for(int v=0;v<(int)r.size();v++)if(!(m>>v&1))if((r[v]&m)!=0&&(r[v]&m)!=m)return false;return true;}
bool middle(const vector<int>&r,int m,int v){int k=__builtin_popcount((unsigned)m);if(k<3||k%2==0||!(m>>v&1)||__builtin_popcount((unsigned)(r[v]&m))!=k/2)return false;int scores=0;for(int x=0;x<(int)r.size();x++)if(m>>x&1)scores|=1<<__builtin_popcount((unsigned)(r[x]&m));return scores==(1<<k)-1;}
int main(int argc,char**argv){int n=9,F=511;ifstream in(argv[1]);string line;getline(in,line);cout<<"{\"candidates\":[";int total=0,primes=0,anchors=0,out=0;while(getline(in,line)){istringstream ss(line);int id;ss>>id;uint64_t code=stoull(line.substr(line.rfind('|')+1));vector<int>r(n),op(n);int bit=35;for(int q=1;q<n;q++)for(int p=0;p<q;p++){if(code>>bit&1)r[p]|=1<<q;else r[q]|=1<<p;bit--;}total++;bool prime=true;for(int m=1;m<F;m++)if(__builtin_popcount((unsigned)m)>1&&module(r,m)){prime=false;break;}if(!prime)continue;primes++;endpoints(r,E);for(int v=0;v<n;v++)op[v]=F^r[v]^(1<<v);endpoints(op,S);bool all=true;for(int v=0;v<n;v++)if((E[F][v]+S[F][v])%2)all=false;if(!all){anchors++;continue;}if(out++)cout<<",";cout<<"{\"node\":"<<id<<",\"code\":"<<code<<",\"rows\":[";for(int v=0;v<n;v++){if(v)cout<<",";cout<<r[v];}cout<<"],\"starts\":[";for(int v=0;v<n;v++){if(v)cout<<",";cout<<S[F][v];}cout<<"],\"ends\":[";for(int v=0;v<n;v++){if(v)cout<<",";cout<<E[F][v];}cout<<"]}";}cout<<"],\"total_merged_nodes\":"<<total<<",\"prime_nodes\":"<<primes<<",\"p_nonzero\":"<<anchors<<",\"p_zero\":"<<out<<"}\n";}