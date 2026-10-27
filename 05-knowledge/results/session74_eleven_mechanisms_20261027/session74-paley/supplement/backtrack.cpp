#include <algorithm>
#include <array>
#include <iostream>
#include <set>
#include <vector>
using namespace std;const int n=11,all=(1<<n)-1;using T=array<int,n>;long tests=0,nodes=0;const long CAP=5000000;bool aborted=false;vector<T> winpaths,winrows;set<T> seen;
int L(const T&t){int a[n];for(int v=0;v<n;v++)a[v]=__builtin_popcount(t[v]);sort(a,a+n);int s=0;for(int v=0;v<n;v++)s+=abs(a[v]-v);return s/2;}
T tau(const T&t,const T&p){T u;for(int v=0;v<n;v++)u[v]=all^(1<<v)^t[v];for(int i=0;i<n-1;i++){int a=p[i],b=p[i+1];u[a]|=1<<b;u[b]&=~(1<<a);}return u;}
bool legal(const T&t,const T&p){for(int i=0;i<n-1;i++)if(!(t[p[i]]>>p[i+1]&1))return false;return true;}
bool solve(const T&t,int want);
bool paths(const T&t,T&p,int d,int used,int want){if(aborted)return false;if(d==n){if(++tests>CAP){aborted=true;return false;}T u=tau(t,p);if(L(u)!=want-1)return false;if(want>1)for(int x:u)if(!x||__builtin_popcount(x)==10)return false;if(solve(u,want-1)){winpaths.push_back(p);winrows.push_back(u);return true;}return false;}int poss=d?t[p[d-1]]&~used:all;while(poss){int v=__builtin_ctz(poss);poss&=poss-1;p[d]=v;if(paths(t,p,d+1,used|(1<<v),want))return true;}return false;}
bool solve(const T&t,int want){if(want==0)return true;if(!seen.insert(t).second)return false;nodes++;if(want==1){T p;for(int i=0;i<n;i++)p[i]=i;sort(p.begin(),p.end(),[&](int a,int b){return pair<int,int>{__builtin_popcount(t[a]),a}<pair<int,int>{__builtin_popcount(t[b]),b};});for(int a=0;a<2;a++){for(int b=0;b<2;b++){if(legal(t,p)){T u=tau(t,p);if(L(u)==0){winpaths.push_back(p);winrows.push_back(u);return true;}}swap(p[n-2],p[n-1]);}swap(p[0],p[1]);}return false;}T p;return paths(t,p,0,0,want);}
int main(){T t={16,2037,513,1015,932,709,21,1093,741,193,893};int startL=L(t);bool ok=solve(t,startL);cout<<"found "<<ok<<" start_L "<<startL<<" tests "<<tests<<" nodes "<<nodes<<" cap "<<CAP<<" aborted "<<aborted<<"\n";reverse(winpaths.begin(),winpaths.end());reverse(winrows.begin(),winrows.end());for(int i=0;i<(int)winpaths.size();i++){cout<<"target_L "<<L(winrows[i])<<"\npath";for(int v:winpaths[i])cout<<" "<<v;cout<<"\nrows";for(int v:winrows[i])cout<<" "<<v;cout<<"\n";}}
