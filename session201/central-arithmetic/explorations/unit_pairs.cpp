#include<bits/stdc++.h>
using namespace std;
struct S{long long D,E,U,d;};S J(S v){return{v.D-v.E+v.U,v.D-v.E,v.D+v.U,v.D};}
long long K(S v,S w){return v.D*w.d+v.d*w.D+v.E*w.E+v.U*w.U-v.E*w.d-v.d*w.E-v.U*w.d-v.d*w.U+v.d*w.d;}
int main(int argc,char**argv){int N=argc>1?atoi(argv[1]):1000001;vector<int>g(N+1);g[1]=1;for(int u=2;u<=N;u++){if(u%2==0)g[u]=g[u/2];else{int t=(u+1)/2;g[u]=g[t]+t-t%2;}}
 vector<unsigned char>hit(N+1,0);
 for(int u=2;1LL*u*u<=2LL*(N-1);u++){S v={g[u],u%2,u,1};for(int k=u;1LL*u*k<=2LL*(N-1);k++){S w={g[k],k%2,k,1};for(int a=0;a<2;a++)for(int b=0;b<2;b++){long long n=K(a?J(v):v,b?J(w):w);if(n<=N)hit[n]|=1<<(a*2+b);}}}
 for(int mask:{1,3,7,15}){int count=0,last=0;vector<int>f;for(int n=3;n<=N;n+=2)if(!(hit[n]&mask)){count++;last=n;if(f.size()<40)f.push_back(n);}cout<<mask<<" gaps "<<count<<" last "<<last<<" first";for(int n:f)cout<<" "<<n;cout<<"\n";}
}
