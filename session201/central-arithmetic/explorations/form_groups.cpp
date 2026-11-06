#include <bits/stdc++.h>
#include "forms.h"
using namespace std;
int main(int argc,char**argv){int N=argc>1?atoi(argv[1]):1000001;int cap=(N-1)/6;vector<int>g(cap+1);g[1]=1;for(int u=2;u<=cap;u++){if(u%2==0)g[u]=g[u/2];else{int v=(u+1)/2;g[u]=g[v]+v-v%2;}}
 vector<uint64_t>hit(N+1,0);for(int L=2;3LL*L*L+1<=N;L++){long long left[]={g[L],L,L%2,1};for(int f=0;f<40;f++){long long c[4]={};for(int i=0;i<4;i++)for(int j=0;j<4;j++)c[j]+=left[i]*forms[f][i][j];for(int K=L;3LL*L*K+1<=N;K++){long long n=c[0]*g[K]+c[1]*K+c[2]*(K%2)+c[3];if(n<=N)hit[n]|=1ULL<<f;}}}
 vector<pair<string,uint64_t>> masks={{"N00",1ULL},{"Nall",15ULL},{"NM",255ULL},{"C3",((1ULL<<40)-1)^255ULL},{"all",(1ULL<<40)-1}};
 for(int f=0;f<40;f++)masks.push_back({"f"+to_string(f),1ULL<<f});
 for(auto [name,mask]:masks){int cnt=0,last=0;vector<int>first;for(int n=3;n<=N;n+=2)if(!(hit[n]&mask)){cnt++;last=n;if(first.size()<35)first.push_back(n);}cout<<name<<" gaps "<<cnt<<" last "<<last<<" first";for(int n:first)cout<<" "<<n;cout<<"\n";}
}
