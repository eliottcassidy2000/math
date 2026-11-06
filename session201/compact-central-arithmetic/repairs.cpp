#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include "central_forms.h"
int main(){const int N=2134117;std::vector<int> g(N+1);g[1]=1;for(int v=2;v<=N;v++){if(v%2==0)g[v]=g[v/2];else{int h=(v+1)/2;g[v]=g[h]+h-h%2;}}
std::vector<uint8_t> best(N+1,255);std::vector<int> index(N+1),left(N+1),right(N+1);
for(int f=0;f<FORMCOUNT;f++)for(int u=2;1LL*u*u*lower4[f]<=4LL*N;u++){long long s[]={g[u],u%2,u,1},r[4]={};for(int i=0;i<4;i++)for(int j=0;j<4;j++)r[j]+=s[i]*forms[f][i][j];for(int v=u;1LL*u*v*lower4[f]<=4LL*N;v++){long long n=r[0]*g[v]+r[1]*(v%2)+r[2]*v+r[3];if(n<=N&&lengths[f]<best[n]){best[n]=lengths[f];index[n]=f;left[n]=u;right[n]=v;}}}
std::cout<<"{\"cap\":"<<N<<",\"gaps5\":[";bool first=true;for(int n=3;n<=N;n+=2)if(best[n]>5){if(!first)std::cout<<",";first=false;std::cout<<n;}std::cout<<"],\"repairs\":[";first=true;for(int n=3;n<=N;n+=2)if(best[n]>5&&best[n]<255){if(!first)std::cout<<",";first=false;std::cout<<"["<<n<<","<<index[n]<<","<<left[n]<<","<<right[n]<<","<<int(best[n])<<"]";}std::cout<<"]}\n";}
