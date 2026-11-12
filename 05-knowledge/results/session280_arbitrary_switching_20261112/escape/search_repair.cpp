#include <bits/stdc++.h>
using namespace std; using ll=long long;
int n,N; vector<int> a; vector<array<ll,12>> dp;vector<ll> hs;
array<ll,3> profile(){
 fill(dp.begin(),dp.end(),array<ll,12>{});fill(hs.begin(),hs.end(),0);hs[0]=1;
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int s=1;s<N;s++){for(int v=0;v<n;v++){auto x=dp[s][v];hs[s]+=x;if(!x)continue;int t=a[v]&~s;while(t){int j=__builtin_ctz(t);t&=t-1;dp[s|1<<j][j]+=x;}}}
 ll h=hs[N-2],H=hs[N-1],V=0;for(int s=1;s<N;s+=2)V+=dp[s][0]*hs[(N-1)^s];return {H,h,V};}
int main(int argc,char**argv){if(argc!=3){cerr<<"usage: search_repair t number_of_flips\n";return 2;}int t=atoi(argv[1]),k=atoi(argv[2]);if(t<1||t>5||k<0||k>5){cerr<<"require1<=t<=5,0<=flips<=5\n";return 2;}n=t+7;N=1<<n;dp.resize(N);hs.resize(N);a.assign(n,0);vector<int>T={52,13,24,33,42,6};int m=t+1;
 a[0]|=1<<1;for(int i=1;i<=m;i++){if(i>1)a[i]|=1;for(int j=i+1;j<=m;j++)a[i]|=1<<j;}
 auto id=[&](int i){return i?m+i:0;};for(int i=0;i<6;i++)for(int j=0;j<6;j++)if(T[i]>>j&1)a[id(i)]|=1<<id(j);
 for(int i=1;i<=m;i++)for(int j=m+1;j<n;j++)a[i]|=1<<j;
 vector<pair<int,int>> es;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)es.push_back({i,j});vector<int> chosen;ll wanted=45+57*(1<<t),vtarget=(1<<(t+1))*57-15;long tested=0,match=0;
 function<void(int,int)>dfs=[&](int from,int left){if(!left){tested++;auto p=profile();if(p[0]==wanted && p[1]==15 && p[2]==vtarget){cout<<"MATCH ";for(int q:chosen)cout<<es[q].first<<","<<es[q].second<<" ";cout<<" rows ";for(int x:a)cout<<x<<",";cout<<endl;match++;}return;}for(int q=from;q<=int(es.size())-left;q++){auto[i,j]=es[q];a[i]^=1<<j;a[j]^=1<<i;chosen.push_back(q);dfs(q+1,left-1);chosen.pop_back();a[i]^=1<<j;a[j]^=1<<i;}};
 dfs(0,k);cerr<<"tested "<<tested<<" matches "<<match<<endl;
}
