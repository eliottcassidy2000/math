#include <bits/stdc++.h>
using namespace std;
int H(vector<int> a) {int n=a.size();vector<array<int,16>> d(1<<n);for(int i=0;i<n;i++)d[1<<i][i]=1;for(int m=1;m<(1<<n);m++)for(int i=0;i<n;i++)if(int x=d[m][i])for(int s=a[i]&~m;s;s&=s-1){int b=s&-s;d[m|b][__builtin_ctz(b)]+=x;}return accumulate(d.back().begin(),d.back().end(),0);}
int main(){vector<vector<int>> ss={{94,88,194,20,4,95,24,123}};int h=63;for(int depth=0;depth<5;depth++){vector<vector<int>>out;h*=3;for(auto a:ss){int n=a.size();for(int m=1;m<(1<<n)-1;m++){auto b=a;b.push_back(((1<<n)-1)^m);for(int i=0;i<n;i++)if(m>>i&1)b[i]|=1<<n;if(H(b)==h){out.push_back(b);cout<<h<<" "<<m<<" :";for(auto x:b)cout<<x<<" ";cout<<"\n";}}}cout<<"hit "<<out.size()<<endl;if(out.empty())break;ss=out;}}
