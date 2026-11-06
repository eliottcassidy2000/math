#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=long long;using G=vector<unsigned>;
void need(bool ok,const char*s){if(!ok)throw runtime_error(s);}
struct Data{vector<I>h;vector<vector<I>>start,end;};
Data profile(const G&a){int n=a.size(),N=1<<n;Data z{vector<I>(N),vector<vector<I>>(N,vector<I>(n)),vector<vector<I>>(N,vector<I>(n))};z.h[0]=1;for(int i=0;i<n;i++)z.start[1<<i][i]=z.end[1<<i][i]=1;for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(S>>i&1){int T=S^(1<<i);if(T)for(int j=0;j<n;j++)if(T>>j&1){if(a[i]>>j&1)z.start[S][i]+=z.start[T][j];else z.end[S][i]+=z.end[T][j];}z.h[S]+=z.start[S][i];}return z;}
array<I,8> data(const G&a){auto z=profile(a);int N=1<<a.size(),F=N-1,T=F^1;I L=0,R=0,beta=0;for(int S=1;S<N;S++){L+=z.end[S][0]*z.h[F^S];R+=z.start[S][0]*z.h[F^S];}for(int S=T;;S=(S-1)&T){beta+=z.h[S]*z.h[T^S];if(!S)break;}I h=z.h[F],s=z.start[F][0],e=z.end[F][0],d=z.h[T];return {h+L-s,e,R+beta,d,s,h,L,R};}
G inflate(G a,int v,G x){int n=a.size(),m=x.size();vector<vector<int>>slot(n);int z=0;for(int i=0;i<n;i++)for(int j=0;j<(i==v?m:1);j++)slot[i].push_back(z++);G b(z);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i]>>j&1)for(int q:slot[i])for(int r:slot[j])b[q]|=1u<<r;for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(x[i]>>j&1)b[slot[v][i]]|=1u<<slot[v][j];return b;}
array<I,5> covers(G a){auto z=profile(a);int N=1<<a.size();vector<I>d(N);d[0]=1;array<I,5>p={I(a.empty()),0,0,0,0};for(int k=1;k<=4;k++){vector<I>e(N);for(int S=1;S<N;S++)for(int T=S;T;T=(T-1)&S)e[S]+=z.h[T]*d[S^T];d=e;p[k]=d.back();}return p;}
int main(){G a={102,0,10,83,135,30,182,47},b={142,0,34,22,71,155,47,94};int cases=0;for(int n=0;n<=4;n++)for(unsigned code=0;code<(1u<<(n*(n-1)/2));code++){G x(n);int bit=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(code>>bit++&1)x[i]|=1u<<j;else x[j]|=1u<<i;auto p=covers(x);auto aa=data(inflate(a,5,x)),bb=data(inflate(b,3,x));I B=93*p[0]+253*p[1]+204*p[2]+48*p[3],V=300*p[0]+1022*p[1]+1178*p[2]+520*p[3]+72*p[4],d=11*p[0]+29*p[1]+21*p[2]+4*p[3],s=5*p[0]+16*p[1]+15*p[2]+4*p[3],h=37*p[0]+111*p[1]+95*p[2]+24*p[3],L=61*p[0]+158*p[1]+124*p[2]+28*p[3],R=100*p[0]+382*p[1]+486*p[2]+236*p[3]+36*p[4];need(aa==array<I,8>{B,0,V,d,s,h,L,R},"A profile");need(bb==array<I,8>{B,0,V,d,s+2*(p[1]+p[2]),h+2*(p[1]+p[2]),L,R},"B profile");cases++;}
cout<<"{\"status\":\"PASS\",\"all_labeled_module_inputs_orders_0_through_4\":"<<cases<<",\"literal_context_graphs\":"<<2*cases<<",\"coverage_claim\":false}\n";}
