#include <algorithm>
#include <cassert>
#include <cstdio>
using namespace std;
int g[7][7];
int h(int const* vs,int n){int d[128][7]={};for(int v=0;v<n;v++)d[1<<v][v]=1;for(int m=1;m<1<<n;m++)for(int v=0;v<n;v++)if(d[m][v])for(int u=0;u<n;u++)if(!(m>>u&1) && g[vs[v]][vs[u]])d[m|1<<u][u]+=d[m][v];int a=0;for(int v=0;v<n;v++)a+=d[(1<<n)-1][v];return a;}
bool strong(const int* vs,int n){for(int i=0;i<n;i++){int s=1<<i;for(int l=1;l<n;l++)for(int j=0;j<n;j++)if(s>>j&1)for(int k=0;k<n;k++)if(g[vs[j]][vs[k]])s|=1<<k;if(s!=(1<<n)-1)return false;}return true;}
void edge(int i,int j,int x){g[i][j]=x;g[j][i]=1-x;}
int main(){edge(0,1,1);edge(1,2,1);edge(2,0,1);int A[]={0,1,2,3,4},B[]={0,1,2,5,6},S[]={0,1,2,3,4,5,6};int lo=1000,atlo=0,total=0;
for(int m=0;m<1<<18;m++){int k=0;for(int j=3;j<7;j++)for(int i=0;i<j;i++)edge(i,j,(m>>k++)&1);
if(!strong(A,5)||!strong(B,5))continue;
int e=h(S,7)-h(A,5)-h(B,5)+3;++total;if(e<lo){lo=e;atlo=0;}if(e==lo)++atlo;
}
printf("%d strong-core unions; minimum mixed weight %d (%d attaining)\n",total,lo,atlo);
assert(total==153664 && lo==10 && atlo==72);
}
