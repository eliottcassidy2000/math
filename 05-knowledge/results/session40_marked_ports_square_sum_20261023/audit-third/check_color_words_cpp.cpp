#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <functional>
using namespace std;
bool traceable(int n){return n==1||n==15||n==16||n==17||n==23||n>=25;}
int color(int a,int b){int n=max(a,b);if(traceable(n))return 1;return (int(sqrt(2*n-1))-int(sqrt(n))==2)?2:3;}
int main(){for(int n=33;n<=36;n++){
 vector<vector<int>>adj(n+1);for(int a=1;a<=n;a++)for(int b=1;b<=n;b++)if(a!=b){int q=sqrt(a+b);if(q*q==a+b)adj[a].push_back(b);}
 vector<int>p(n);long long total=0,nodes=0;map<string,int>words;set<string>hist;vector<int>first,second;map<string,vector<int>>one;
 function<void(int,unsigned long long)> dfs=[&](int k,unsigned long long seen){nodes++;if(k==n){if(p[0]>p[n-1])return;total++;string key;key.push_back(p[0]);key.push_back(p[n-1]);int count[4]={};for(int i=0;i<n-1;i++){int c=color(p[i],p[i+1]);key.push_back(c);count[c]++;}words[key]++;if(words[key]==1)one[key]=p;else if(first.empty()){first=one[key];second=p;}string h;h.push_back(p[0]);h.push_back(p[n-1]);for(int c=1;c<=3;c++)h.push_back(count[c]);hist.insert(h);return;}
 for(int v:adj[p[k-1]])if(!(seen>>(v-1)&1)){p[k]=v;dfs(k+1,seen|(1ULL<<(v-1)));}};
 for(int s=1;s<=n;s++){p[0]=s;dfs(1,1ULL<<(s-1));}
 int maximum=0;for(auto z:words)maximum=max(maximum,z.second);
 cout<<n<<" "<<total<<" "<<words.size()<<" "<<hist.size()<<" "<<maximum<<" "<<nodes<<"\n";
 if(!first.empty()){for(int x:first)cout<<x<<",";cout<<"\n";for(int x:second)cout<<x<<",";cout<<"\n";}
 }}
