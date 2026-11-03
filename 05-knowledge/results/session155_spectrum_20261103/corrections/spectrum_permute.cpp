#include <algorithm>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
int main(){for(int n=1;n<=6;++n){int m=n*(n-1)/2;set<int> counts; for(unsigned g=0;g<(1u<<m);++g){bool a[6][6]{};int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit){if(g>>bit&1) a[j][i]=true; else a[i][j]=true;} vector<int> p(n); iota(p.begin(),p.end(),0);int H=0;do{bool ok=true;for(int k=1;k<n;++k)if(!a[p[k-1]][p[k]]){ok=false;break;}H+=ok;}while(next_permutation(p.begin(),p.end())); counts.insert(H);}cout<<n<<':';for(int h:counts)cout<<' '<<h;cout<<'\n';}}
