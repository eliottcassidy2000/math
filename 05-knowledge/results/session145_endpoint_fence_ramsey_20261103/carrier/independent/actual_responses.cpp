#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;
using I=long long;
I checks=0;
void need(bool b,const char*s){++checks;if(!b)throw runtime_error(s);}
I count_h(const vector<int>&a){int n=a.size(),lim=1<<n;I dp[512][9]{};for(int v=0;v<n;v++)dp[1<<v][v]=1;for(int m=1;m<lim;m++)for(int v=0;v<n;v++)if(dp[m][v]){int todo=a[v]&~m;while(todo){int w=__builtin_ctz((unsigned)todo);todo&=todo-1;dp[m|(1<<w)][w]+=dp[m][v];}}I s=0;for(int v=0;v<n;v++)s+=dp[lim-1][v];return s;}
vector<int> source(int n,int m){vector<int>a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(m>>k&1)a[i]|=1<<j;else a[j]|=1<<i;return a;}
void birth(vector<int>&a,int cut){int n=a.size();for(int v=0;v<n;v++)if(!(cut>>v&1))a[v]|=1<<n;a.push_back(cut);}
void walsh(vector<I>&a){for(int k=1;k<(int)a.size();k*=2)for(int b=0;b<(int)a.size();b+=2*k)for(int j=b;j<b+k;j++){I p=a[j],q=a[j+k];a[j]=p+q;a[j+k]=p-q;}}
int main(int argc,char**argv){if(argc!=2)throw runtime_error("need coefficient map");ifstream f(argv[1]);int m,key,sgn;vector<int> keys(32768,-1),signs(32768);while(f>>m>>key>>sgn){keys[m]=key;signs[m]=sgn;}need(count_if(keys.begin(),keys.end(),[](int k){return k>=0;})==1486,"forest census");vector<vector<int>>contexts={{},{0},{1},{2},{4},{8},{16},{32},{63},{21},{42},{0,1},{1,2},{21,85},{42,42},{63,127},{17,69},{54,37}};cout<<"{\"n\":6,\"sources\":32768,\"queries\":[";int qi=0;I all_graphs=0;for(auto cuts:contexts){vector<I>v(32768);I d1=0,d2=0,mean_sum=0,second_sum=0;for(int a=0;a<32768;a++){auto adj=source(6,a);for(int c:cuts)birth(adj,c);I h=count_h(adj);v[a]=h;int delta=((__builtin_popcount((unsigned)(a&10376))%2)?-1:1)-((__builtin_popcount((unsigned)(a&5392))%2)?-1:1);d1+=delta*h;d2+=delta*h*h;mean_sum+=h;second_sum+=h*h;all_graphs++;}need(d1==0,"new positive laws have same mean");walsh(v);vector<I>coeff(571);vector<bool>seen(571);for(int a=0;a<32768;a++){if(keys[a]<0)need(v[a]==0,"forbidden support");else {int k=keys[a];I q=signs[a]*v[a];if(seen[k])need(coeff[k]==q,"endpoint parity coefficient identity");else {seen[k]=true;coeff[k]=q;}}}if(qi++)cout<<",";cout<<"{\"cuts\":[";for(int j=0;j<(int)cuts.size();j++){if(j)cout<<",";cout<<cuts[j];}cout<<"],\"H_sum\":"<<mean_sum<<",\"H2_sum\":"<<second_sum<<",\"density_delta_H\":"<<d1<<",\"density_delta_H2\":"<<d2<<"}";}cout<<"],\"actual_graphs_counted\":"<<all_graphs<<",\"checks\":"<<checks<<"}\n";}
