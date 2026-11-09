#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;
int main(int argc,char**argv){int n=argc>1?atoi(argv[1]):6,m=n*(n-1)/2; vector<int>pw(m+1,1); for(int i=0;i<m;i++)pw[i+1]=pw[i]*3; int N=pw[m]; int ix[6][6]; vector<pair<int,int>>pairs;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ix[i][j]=ix[j][i]=pairs.size();pairs.push_back({i,j});}
vector<int>p(n);iota(p.begin(),p.end(),0);vector<vector<int>>perms;vector<int>h(N,0);do{perms.push_back(p);int c=N-1;for(int k=1;k<n;k++){int a=p[k-1],b=p[k];c-=(2-(a<b?0:1))*pw[ix[a][b]];}h[c]++;}while(next_permutation(p.begin(),p.end()));
// Upward zeta: a path graph contributes to every oriented supergraph.
for(int k=0;k<m;k++)for(int b=0;b<N;b+=3*pw[k])for(int a=0;a<pw[k];a++){int v=h[b+a+2*pw[k]];h[b+a]+=v;h[b+a+pw[k]]+=v;}
vector<int>s=h; // Multiaffinity: f(both)=f(forward)+f(backward)-f(empty).
for(int k=0;k<m;k++)for(int b=0;b<N;b+=3*pw[k])for(int a=0;a<pw[k];a++)s[b+a+2*pw[k]]=s[b+a]+s[b+a+pw[k]]-s[b+a+2*pw[k]];
map<pair<int,int>,long long>joint,by_missing;map<int,long long>hist;vector<char>seen(N,0);ofstream reps("classes"+to_string(n)+".csv");reps<<"code,H,missing_pairs,orbit_size\n";map<int,int>cls;int maxdiff=-1,maxcode=0;long long eq=0;
for(int c=0;c<N;c++){assert(s[c]>=h[c]);joint[{h[c],s[c]}]++;hist[h[c]]++;int t=c,d[15],missing=0;for(int k=0;k<m;k++){d[k]=t%3;t/=3;missing+=d[k]==2;}by_missing[{h[c],missing}]++;if(h[c]==s[c])eq++;if(s[c]-h[c]>maxdiff){maxdiff=s[c]-h[c];maxcode=c;}
if((h[c]==7||h[c]==21)&&!seen[c]){int orbit=0;for(auto &q:perms){int cc=0;for(int k=0;k<m;k++){auto [i,j]=pairs[k];int a=q[i],b=q[j],v=d[k];if(v!=2&&a>b)v=1-v;cc+=v*pw[ix[a][b]];}if(!seen[cc]){seen[cc]=1;orbit++;}assert(h[cc]==h[c]);}cls[h[c]]++;reps<<c<<","<<h[c]<<","<<missing<<","<<orbit<<"\n";}}
ofstream out("hist"+to_string(n)+".csv");out<<"H,count\n";for(auto [a,b]:hist)out<<a<<","<<b<<"\n";ofstream j("joint"+to_string(n)+".csv");j<<"H_oriented,H_semicomplete,count\n";for(auto [a,b]:joint)j<<a.first<<","<<a.second<<","<<b<<"\n";ofstream z("missing"+to_string(n)+".csv");z<<"H,missing_pairs,count\n";for(auto [a,b]:by_missing)z<<a.first<<","<<a.second<<","<<b<<"\n";
cout<<"n="<<n<<" N="<<N<<" H7="<<hist[7]<<" H21="<<hist[21]<<" classes7="<<cls[7]<<" classes21="<<cls[21]<<" equal="<<eq<<" max_difference="<<maxdiff<<" max_code="<<maxcode<<"\n";
// Binary arrays used only by independent controls, reproducible and not versioned.
ofstream bin("counts"+to_string(n)+".bin",ios::binary);bin.write((char*)h.data(),N*sizeof(int));bin.write((char*)s.data(),N*sizeof(int));
}
