#include <bits/stdc++.h>
using namespace std;
int aa[10],memo[1024][10], n;
int dfs(int mask,int last){if(mask==((1<<n)-1))return 1;int &v=memo[mask][last];if(v!=-1)return v;v=0;for(int j=0;j<n;j++)if(!(mask>>j&1)&&(aa[last]>>j&1))v+=dfs(mask|1<<j,j);return v;}
int countpaths(vector<int> a){n=a.size();copy(a.begin(),a.end(),aa);memset(memo,-1,sizeof memo);int s=0;for(int i=0;i<n;i++)s+=dfs(1<<i,i);return s;}
vector<int> ear(vector<int> a,int m){int n=a.size();a.push_back(m);for(int i=0;i<n;i++)if(!(m>>i&1))a[i]|=1<<n;return a;}
int main(){vector<int>a(5);int k=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++,k++)if(759>>k&1)a[i]|=1<<j;else a[j]|=1<<i;a=ear(a,2);for(auto [pm,target]: vector<pair<int,int>>{{50,283},{24,287},{24,289}}){auto b=ear(a,pm);int pairs=0,hits=0;for(int u=0;u<128;u++){auto c=ear(b,u);if(countpaths(c)>target)continue;for(int v=0;v<128;v++){auto d=ear(c,v);if(countpaths(d)>target)continue;pairs++;for(int w=0;w<512;w++)if(countpaths(ear(d,w))==target)hits++;}}cout<<"mask="<<pm<<" target="<<target<<" pairs="<<pairs<<" direct_masks="<<pairs*512<<" hits="<<hits<<endl;}}
