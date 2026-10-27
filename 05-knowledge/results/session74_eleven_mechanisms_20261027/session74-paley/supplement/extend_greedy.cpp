#include <algorithm>
#include <array>
#include <iostream>
using namespace std;const int n=11,all=(1<<n)-1;array<int,n> base={16,1533,1801,1617,964,541,421,45,681,195,1009},cur,path;int want;long checked,eligible;
int L(array<int,n>t){int a[n];for(int v=0;v<n;v++)a[v]=__builtin_popcount(t[v]);sort(a,a+n);int s=0;for(int v=0;v<n;v++)s+=abs(a[v]-v);return s/2;}
bool test(){checked++;for(int v=0;v<n;v++)cur[v]=all^(1<<v)^base[v];for(int i=0;i<n-1;i++){int a=path[i],b=path[i+1];cur[a]|=1<<b;cur[b]&=~(1<<a);}if(L(cur)!=want)return false;eligible++;if(want)for(int v:cur)if(!v||__builtin_popcount(v)==n-1)return false;cout<<"target_L "<<want<<" checked "<<checked<<" eligible "<<eligible<<"\npath";for(int v:path)cout<<" "<<v;cout<<"\nrows";for(int v:cur)cout<<" "<<v;cout<<"\n"<<flush;return true;}
bool rec(int d,int used){if(d==n)return test();int poss=d?base[path[d-1]]&~used:all;while(poss){int v=__builtin_ctz(poss);poss&=poss-1;path[d]=v;if(rec(d+1,used|(1<<v)))return true;}return false;}
int main(){for(want=8;want>=0;want--){checked=eligible=0;if(!rec(0,0)){cout<<"NO_GREEDY_EXTENSION target_L "<<want<<" checked "<<checked<<" eligible "<<eligible<<"\n";return 0;}base=cur;}}
