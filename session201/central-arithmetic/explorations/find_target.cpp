#include<bits/stdc++.h>
#include "forms.h"
using namespace std;
int main(int argc,char**argv){int n=atoi(argv[1]);int cap=(n-1)/6;vector<int>g(cap+1);g[1]=1;for(int u=2;u<=cap;u++){if(u%2==0)g[u]=g[u/2];else{int x=(u+1)/2;g[u]=g[x]+x-x%2;}}
 for(int L=2;3LL*L*L+1<=n;L++){long long a[]={g[L],L,L%2,1};for(int f=0;f<40;f++){long long c[4]={};for(int i=0;i<4;i++)for(int j=0;j<4;j++)c[j]+=a[i]*forms[f][i][j];for(int K=L;3LL*L*K+1<=n;K++)if(c[0]*g[K]+c[1]*K+c[2]*(K%2)+c[3]==n){cout<<n<<" "<<L<<" "<<K<<" "<<f<<" gL "<<g[L]<<" gK "<<g[K]<<" row";for(auto v:c)cout<<" "<<v;cout<<"\n";return 0;}}}
 cout<<"ABSENT\n";
}
