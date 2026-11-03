#include <bits/stdc++.h>
using namespace std;
long count(int m,unsigned x,unsigned y,bool xy){long d[3][3][2]={};d[0][0][0]=1;for(int i=0;i<m;i++){long e[3][3][2]={};int a=x>>i&1,b=y>>i&1;for(int u=0;u<3;u++)for(int v=0;v<3;v++)for(int w=0;w<2;w++){auto t=d[u][v][w];if(!t)continue;e[a?1:2][v][w]+=t; if(v || a)e[u][b?1:2][w]+=t; if(w||b)e[u][v][1]+=t;}memcpy(d,e,sizeof d);}long s=0;for(int u: {0,2})for(int v:{0,2})for(int w=0;w<2;w++)if(v||xy)s+=d[u][v][w];return s;}
int main(){mt19937_64 rng(293429);for(long k=0;k<4000000;k++){int m=12+k%8;unsigned mask=(1u<<m)-1,x=rng()&mask,y=rng()&mask; bool xy=rng()&1;long h=count(m,x,y,xy)+count(m,y,x,!xy);if(h==80269){cout<<"found m="<<m<<" x="<<x<<" y="<<y<<" xy="<<xy<<" k="<<k<<endl;return 0;}if(k%1000000==0)cerr<<k<<endl;}cout<<"none\n";}
