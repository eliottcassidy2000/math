#include <bits/stdc++.h>
using namespace std; using I=unsigned long long;
int main(int argc,char**argv){ifstream f(argv[1]);string line;map<string,I> target;int edges=0;
while(getline(f,line)){istringstream in(line);string w,v,kind,map;I h,k;int a,b;if(!(in>>w>>h>>v>>k>>kind>>map>>a>>b))throw runtime_error("parse");for(auto [s,c]:vector<pair<string,I>>{{w,h},{v,k}}){if(target.count(s)&&target[s]!=c)throw runtime_error("inconsistent");target[s]=c;}edges++;}
for(auto [w,h]:target){int n=w.size()+2;if(n>12)throw runtime_error("bound");vector<int>out(n);for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){bool backward=b>=2&&a<2&&((w[b-2]-'0')&(1<<a));out[backward?b:a]|=1<<(backward?a:b);}
int N=1<<n;vector<I>d(N*n);for(int mask=1;mask<N;mask++)for(int a=0;a<n;a++)if(mask>>a&1){int rest=mask^(1<<a);if(!rest)d[mask*n+a]=1;else for(int q=rest&out[a];q;q&=q-1)d[mask*n+a]+=d[rest*n+__builtin_ctz(q)];}I got=0;for(int a=0;a<n;a++)got+=d[(N-1)*n+a];if(got!=h)throw runtime_error(w);}
cout<<"PASS independent reverse-path recurrence: "<<target.size()<<" graphs, "<<edges<<" edges\n";}
