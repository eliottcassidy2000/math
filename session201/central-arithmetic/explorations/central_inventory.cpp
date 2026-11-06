#include<bits/stdc++.h>
#include "central_forms.h"
using namespace std;
int main(int argc,char**argv){int N=argc>1?atoi(argv[1]):1000001, maxlen=argc>2?atoi(argv[2]):5;
 vector<int>g(N+1);g[1]=1;for(int u=2;u<=N;u++){if(u%2==0)g[u]=g[u/2];else{int t=(u+1)/2;g[u]=g[t]+t-t%2;}}
 vector<uint8_t>minimum(N+1,255);
 for(int f=0;f<FORMCOUNT;f++){if(length[f]>maxlen)continue;
  for(int u=2;1LL*u*u*lower4[f]<=4LL*N;u++){long long v[]={g[u],u%2,u,1},row[4]={};for(int i=0;i<4;i++)for(int j=0;j<4;j++)row[j]+=v[i]*forms[f][i][j];
   for(int k=u;1LL*u*k*lower4[f]<=4LL*N;k++){long long n=row[0]*g[k]+row[1]*(k%2)+row[2]*k+row[3];if(n<=N)minimum[n]=min<int>(minimum[n],length[f]);}}
 }
 for(int s=0;s<=maxlen;s++){int count=0,last=0;vector<int>first;for(int n=3;n<=N;n+=2)if(minimum[n]>s){count++;last=n;if(first.size()<50)first.push_back(n);}cout<<"len "<<s<<" gaps "<<count<<" last "<<last<<" first";for(int n:first)cout<<" "<<n;cout<<"\n";}
}
