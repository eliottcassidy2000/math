#include <bits/stdc++.h>
#include "forms.h"
using namespace std;
int main(int argc,char**argv){int N=argc>1?atoi(argv[1]):1000001,B=argc>2?atoi(argv[2]):32;int cap=(N-1)/6;
 vector<int> g(cap+1);g[1]=1;for(int u=2;u<=cap;u++){if(u%2==0)g[u]=g[u/2];else {int v=(u+1)/2;g[u]=g[v]+v-v%2;}}
 vector<uint16_t> minimum(N+1,0);
 for(int L=2;L<=B;L++){long long left[]={g[L],L,L%2,1};
  for(int f=0;f<40;f++){long long c[4]={};for(int i=0;i<4;i++)for(int j=0;j<4;j++)c[j]+=left[i]*forms[f][i][j];
   for(int K=L;3LL*L*K+1<=N;K++){long long n=c[0]*g[K]+c[1]*K+c[2]*(K%2)+c[3];if(n<=N&&minimum[n]==0)minimum[n]=L;}}
  if((L&(L-1))==0||L==B){int count=0,last=0;vector<int> first;for(int n=3;n<=N;n+=2)if(!minimum[n]){count++;last=n;if(first.size()<35)first.push_back(n);}cout<<"B "<<L<<" gaps "<<count<<" last "<<last<<" first";for(int n:first)cout<<" "<<n;cout<<"\n";}
 }
 map<int,int>frequency;for(int n=3;n<=N;n+=2)frequency[minimum[n]]++;
 cout<<"min_L_hist";for(auto [l,c]:frequency)cout<<" "<<l<<":"<<c;cout<<"\n";
}
