#include <vector>
#include <iostream>
#include <cstdint>
using namespace std;
using U=uint64_t;
vector<int> birth(vector<int> a,int cut){int n=a.size();for(int i=0;i<n;i++)if(cut>>i&1)a[i]|=1<<n;a.push_back(((1<<n)-1)^cut);return a;}
pair<U,U> fields(const vector<int>&a,int w){int n=a.size(),N=1<<n;vector<U> dp(N*n),h(N);h[0]=1;for(int s=1;s<N;s++){for(int v=0;v<n;v++)if(s>>v&1){int r=s^(1<<v);U val=r?0:1;for(int u=0;u<n;u++)if((r>>u&1)&&(a[u]>>v&1))val+=dp[r*n+u];dp[s*n+v]=val;h[s]+=val;}}U f=0;for(int s=1;s<N;s++)if(s>>w&1)f+=dp[s*n+w]*h[(N-1)^s];return {h[N-1],f};}
bool strong(const vector<int>&a){for(int v=0;v<(int)a.size();v++){int seen=1<<v,old=0;while(old!=seen){old=seen;for(int u=0;u<(int)a.size();u++)if(seen>>u&1)seen|=a[u];}if(seen!=(1<<(int)a.size())-1)return false;}return true;}
int main(){vector<vector<int>> sources={{56,33,51,54,2,16},{56,1,19,22,2,30}};for(auto a:sources){auto[h,f]=fields(a,0);int n=a.size();int hits[2]={};cout<<"SOURCE "<<h<<" "<<f<<"\n";for(int x=0;x<(1<<n);x++){auto b=birth(a,x);for(int y=0;y<(1<<(n+1));y++){auto c=birth(b,y);if(!strong(c))continue;auto[H,F]=fields(c,0);for(int k=0;k<2;k++)if(H==3*h+2+2*k && F==3*f){hits[k]++;cout<<"HIT "<<2+2*k<<" "<<x<<" "<<y<<" ROWS";for(auto z:c)cout<<" "<<z;cout<<"\n";}}}cout<<"COUNTS "<<hits[0]<<" "<<hits[1]<<"\n";}}
