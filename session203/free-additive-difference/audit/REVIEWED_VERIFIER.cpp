#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=long long;using G=vector<unsigned>;
void need(bool ok,const char*s){if(!ok)throw runtime_error(s);}
struct Data{vector<I>h;vector<vector<I>>start,end;};
Data profile(const G&a){int n=a.size(),N=1<<n;Data z{vector<I>(N),vector<vector<I>>(N,vector<I>(n)),vector<vector<I>>(N,vector<I>(n))};z.h[0]=1;for(int i=0;i<n;i++)z.start[1<<i][i]=z.end[1<<i][i]=1;for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(S>>i&1){int T=S^(1<<i);if(T)for(int j=0;j<n;j++)if(T>>j&1){if(a[i]>>j&1)z.start[S][i]+=z.start[T][j];else z.end[S][i]+=z.end[T][j];}z.h[S]+=z.start[S][i];}return z;}
G glue(const G&core,const G&u){int c=core.size(),n=u.size();G a=core;a.resize(c+n-1);for(int j=1;j<n;j++){int x=c+j-1;for(int i=0;i<c;i++)if(i==0?(u[0]>>j&1):true)a[i]|=1u<<x;else a[x]|=1u<<i;for(int k=1;k<j;k++){int y=c+k-1;if(u[j]>>k&1)a[x]|=1u<<y;else a[y]|=1u<<x;}}return a;}
G endpoint_source(int s){int m=0;while((1<<m)<=s)m++;G a(m+3);int x=m+1,y=m+2;for(int i=1;i<=m;i++){for(int j=i+1;j<=m;j++)a[i]|=1u<<j;a[i]|=1u<<x;a[y]|=1u<<i;if(s>>(m-i)&1)a[0]|=1u<<i;else a[i]|=1u;}a[x]|=(1u<<y)|1u;a[y]|=1u;return a;}
int main(){int all=0,free=0;I maxh=0;for(int n=1;n<=4;n++)for(int code=0;code<(1<<(n*(n-1)/2));code++){G u(n);int b=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){if(code>>b++&1)u[i]|=1u<<j;else u[j]|=1u<<i;}auto z=profile(u);int full=(1<<n)-1;I h=z.h[full],s=z.start[full][0],L=0;for(int S=1;S<=full;S++)if(S&1)L+=z.end[S][0]*z.h[full^S];I base=profile(glue({2,0},u)).h.back(),parent=profile(glue({2,0,11,19,7},u)).h.back(),child=profile(glue({2,8,3,5},u)).h.back();need(base==h+L-s,"edge profile");need(parent==3*base,"triangle parent");need(child==parent+2*s,"additive identity");all++;}
 for(int s=1;s<=127;s++){G u=endpoint_source(s);auto z=profile(u);need(z.start.back()[0]==s,"free endpoint");I base=profile(glue({2,0},u)).h.back(),parent=profile(glue({2,0,11,19,7},u)).h.back(),child=profile(glue({2,8,3,5},u)).h.back();need(parent==3*base&&child-parent==2*s,"free additive");need(child>parent&&parent%2==1&&child%2==1,"odd positive outputs");maxh=max(maxh,child);free++;}
 cout<<"{\"status\":\"PASS\",\"all_labelled_root_zero_inputs_orders_1_through_4\":"<<all<<",\"freely_prescribed_corrections_2s_for_s_1_through_127\":"<<free<<",\"largest_child_count\":"<<maxh<<",\"coverage_claim\":false}\n";
}
