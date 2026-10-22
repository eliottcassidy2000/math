#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
int main(){int n;string bits;while(cin>>n>>bits){vector<int>s(n),p(n);for(int&i:s)cin>>i;int a[10][10]={};int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){a[i][j]=bits[k++]=='1'?1:-1;a[j][i]=-a[i][j];}iota(p.begin(),p.end(),0);long long aut=0,anti=0,c=0,H=0,blue=0,f=0;do{bool au=true,an=true,path=true;for(int i=0;i<n;i++)for(int j=i+1;j<n&&(au||an);j++){au&=a[p[i]][p[j]]==a[i][j];an&=a[p[i]][p[j]]==-a[i][j];}aut+=au;anti+=an;if(au){bool comm=true;for(int i=0;i<n;i++)comm&=p[s[i]]==s[p[i]];c+=comm;}for(int i=1;i<n;i++)path&=a[p[i-1]][p[i]]==1;if(path){H++;bool sym=true,fixed=true;for(int i=0;i<n;i++){fixed&=s[p[i]]==p[n-1-i];for(int j=i+1;j<n;j++)sym&=a[p[n-1-i]][p[n-1-j]]==-a[p[i]][p[j]];}blue+=sym;f+=fixed;}}while(next_permutation(p.begin(),p.end()));cout<<aut<<' '<<anti<<' '<<c<<' '<<H<<' '<<blue<<' '<<f<<'\n';}}
