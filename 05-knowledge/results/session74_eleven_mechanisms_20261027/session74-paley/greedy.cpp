#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;const int n=11,all=(1<<n)-1;
array<int,n> base,cur,path;int dp[1<<n],lastv[1<<n],depthWanted;long checked,eligible; vector<int> foundorder;
int root(int x,int*p){while(p[x]!=x)x=p[x];return x;}
int L(array<int,n>t){int a[n];for(int v=0;v<n;v++)a[v]=__builtin_popcount(t[v]);sort(a,a+n);int s=0;for(int v=0;v<n;v++)s+=abs(a[v]-v);return s/2;}
bool test(){checked++;for(int v=0;v<n;v++)cur[v]=all^(1<<v)^base[v];for(int i=0;i<n-1;i++){int a=path[i],b=path[i+1];cur[a]|=1<<b;cur[b]&=~(1<<a);}if(L(cur)!=15-depthWanted)return false;eligible++;dp[0]=0;for(int s=1;s<=all;s++){dp[s]=100;int todo=s;while(todo){int v=__builtin_ctz(todo);todo&=todo-1;int r=s^(1<<v);int val=dp[r]+__builtin_popcount(cur[v]&r);if(val<dp[s]){dp[s]=val;lastv[s]=v;}}}if(dp[all]>16-depthWanted)return false;vector<int>order;int s=all;while(s){int v=lastv[s];order.push_back(v);s^=1<<v;}reverse(order.begin(),order.end());int p[n];iota(p,p+n,0);int c=n,m=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(cur[order[j]]>>order[i]&1){int x=root(order[i],p),y=root(order[j],p);if(x!=y){p[x]=y;c--;}m++;}if(m-n+c!=6-depthWanted)return false;foundorder=order;cout<<"depth "<<depthWanted<<" checked "<<checked<<" eligible "<<eligible<<" L "<<L(cur)<<" f "<<dp[all]<<" c "<<c<<" beta "<<m-n+c<<"\npath";for(int v:path)cout<<" "<<v;cout<<"\norder";for(int v:order)cout<<" "<<v;cout<<"\nrows";for(int v:cur)cout<<" "<<v;cout<<"\n"<<flush;return true;}
bool rec(int d,int used){if(d==n)return test();int poss=d?base[path[d-1]]&~used:all;while(poss){int v=__builtin_ctz(poss);poss&=poss-1;path[d]=v;if(rec(d+1,used|(1<<v)))return true;}return false;}
int main(){int sq=0;for(int x=1;x<n;x++)sq|=1<<(x*x%n);for(int x=0;x<n;x++)for(int y=0;y<n;y++)if(sq>>((y-x+n)%n)&1)base[x]|=1<<y;for(depthWanted=1;depthWanted<=6;depthWanted++){checked=eligible=0;if(!rec(0,0)){cout<<"NO_GREEDY_EXTENSION depth "<<depthWanted<<" checked "<<checked<<" eligible "<<eligible<<"\n";break;}base=cur;}}
