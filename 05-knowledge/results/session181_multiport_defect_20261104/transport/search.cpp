#include <bits/stdc++.h>
using namespace std;
using ll=long long;
pair<ll,ll> response(vector<int> r){int n=r.size(),N=1<<n;vector<array<ll,8>> dp(N);vector<ll> h(N);h[0]=1;for(int s=1;s<N;s++){for(int v=0;v<n;v++)if(s>>v&1){int t=s^(1<<v);if(!t)dp[s][v]=1;else for(int u=0;u<n;u++)if(r[u]>>v&1)dp[s][v]+=dp[t][u];h[s]+=dp[s][v];}}ll f=0;for(int s=1;s<N;s++)f+=dp[s][0]*h[(N-1)^s];return {h.back(),f};}
vector<int> birth(vector<int>r,int mask){int n=r.size();for(int i=0;i<n;i++)if(mask>>i&1)r[i]|=1<<n;r.push_back(((1<<n)-1)^mask);return r;}
int main(){vector<vector<int>> seeds={{56,33,51,54,2,16},{56,1,19,22,2,30}};for(auto r:seeds){auto [h,f]=response(r);int hits=0;cout<<"SOURCE "<<h<<" "<<f<<"\n";for(int a=0;a<64;a++)for(int b=0;b<128;b++){auto s=birth(birth(r,a),b);auto [H,F]=response(s);if(F==3*f&&(H==3*h+2||H==3*h+4)){cout<<"HIT "<<a<<" "<<b<<" "<<H<<" "<<F<<" ROWS";for(int x:s)cout<<" "<<x;cout<<"\n";hits++;}}cout<<"TESTED 8192 HITS "<<hits<<"\n";}}
