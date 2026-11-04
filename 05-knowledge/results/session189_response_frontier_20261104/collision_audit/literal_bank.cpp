#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
using namespace std;
int main(){int cases;cin>>cases;while(cases--){int id,n;cin>>id>>n;vector<int>rows(n);for(int &x:rows)cin>>x;int x=n,y=n+1;vector<int>mask{0};for(int s=1;s<(1<<n);s++)if(__builtin_popcount((unsigned)s)==1)mask.push_back(s);for(int s=1;s<(1<<n);s++)if(__builtin_popcount((unsigned)s)==2)mask.push_back(s);int L=mask.size();vector<int>ix(1<<n,-1);for(int j=0;j<L;j++)ix[mask[j]]=j;vector<long long>co(L*L);vector<int>p(n+2);iota(p.begin(),p.end(),0);
do{bool ok=true;int am=0,cm=0;vector<pair<int,int>>neg;for(int j=0;j<n+1;j++){int u=p[j],v=p[j+1];if(u<n&&v<n){if(!(rows[u]>>v&1)){ok=false;break;}}else if(u>=n&&v>=n){if(u!=x||v!=y){ok=false;break;}}else if(u<n){if(v==x)am|=1<<u;else cm|=1<<u;}else neg.push_back({u,v});}if(!ok)continue;for(int s=0;s<(1<<neg.size());s++){int a=am,c=cm;for(int j=0;j<(int)neg.size();j++)if(s>>j&1){if(neg[j].first==x)a|=1<<neg[j].second;else c|=1<<neg[j].second;}if(ix[a]<0||ix[c]<0){cerr<<"bidegree violation";return 2;}co[ix[a]*L+ix[c]]+=(__builtin_popcount((unsigned)s)&1)?-1:1;}}while(next_permutation(p.begin(),p.end()));
cout<<id<<" "<<L;for(auto z:co)cout<<" "<<z;cout<<"\n";cerr<<"literal case "<<id<<" complete\n";}}
