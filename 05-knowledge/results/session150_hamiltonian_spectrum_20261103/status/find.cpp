#include <iostream>
#include <random>
using namespace std;
int main(){mt19937_64 r(17); bool f[3]={}; int goals[3]={49,63,343}; for(int i=0;i<2000000 && !(f[0]&&f[1]&&f[2]);i++){int n=8; unsigned long long code=r()&((1ULL<<28)-1); int out[8]={},k=0; for(int a=0;a<n;a++)for(int b=a+1;b<n;b++,k++) if((code>>k)&1)out[a]|=1<<b;else out[b]|=1<<a; int dp[256][8]={};for(int a=0;a<n;a++)dp[1<<a][a]=1;for(int s=1;s<256;s++)for(int v=0;v<n;v++)if(s>>v&1){int t=out[v]&~s;while(t){int w=__builtin_ctz(t);t&=t-1;dp[s|1<<w][w]+=dp[s][v];}} int h=0;for(int a=0;a<n;a++)h+=dp[255][a];bool strong=true; for(int a=0;a<n;a++){int seen=1<<a, prev=0;while(seen!=prev){prev=seen;for(int b=0;b<n;b++)if(seen>>b&1)seen|=out[b];}if(seen!=255)strong=false;} for(int j=0;j<3;j++)if(strong&&!f[j]&&h==goals[j]){f[j]=1;cout<<h<<" "<<code<<"\n";}}
}
