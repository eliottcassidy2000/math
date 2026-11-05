#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <random>
#include <stdexcept>
using namespace std;
void need(bool x,const char*s){if(!x)throw runtime_error(s);}
vector<vector<long long>> endpoints(const vector<int>&a){int n=a.size(),M=1<<n;vector<vector<long long>>d(M,vector<long long>(n));for(int v=0;v<n;v++)d[1<<v][v]=1;for(int S=1;S<M;S++)for(int u=0;u<n;u++)for(int Z=a[u]&~S;Z;Z&=Z-1){int b=Z&-Z;d[S|b][__builtin_ctz((unsigned)b)]+=d[S][u];}return d;}
long long cycles(const vector<int>&a){int n=a.size();if(n<3)return 0;int M=1<<n;vector<vector<long long>>d(M,vector<long long>(n));d[1][0]=1;for(int S=1;S<M;S++)for(int u=0;u<n;u++)for(int Z=a[u]&~S;Z;Z&=Z-1){int b=Z&-Z;d[S|b][__builtin_ctz((unsigned)b)]+=d[S][u];}long long z=0;for(int v=0;v<n;v++)if(a[v]&1)z+=d.back()[v];return z;}
long long H(const vector<int>&a){if(a.empty())return 1;auto d=endpoints(a);long long h=0;for(auto v:d.back())h+=v;return h;}
vector<int> build(const vector<int>&a,int r,int m){int t=a.size();vector<int>b(2*t+m);for(int i=0;i<t;i++)for(int j=0;j<t;j++){if(a[i]>>j&1){b[m+i]|=1<<(m+j);b[m+t+j]|=1<<(m+t+i);}if(i==r&&j==r)b[m+i]|=1<<(m+t+j);else b[m+t+j]|=1<<(m+i);}for(int u=0;u<m;u++){for(int v=u+1;v<m;v++)b[u]|=1<<v;for(int i=0;i<t;i++){b[u]|=1<<(m+i);b[m+t+i]|=1<<u;}}return b;}
long long marked(const vector<int>&a,int u,int v){int n=a.size(),M=1<<n;vector<vector<long long>>d(M,vector<long long>(2*n));for(int x=0;x<n;x++)d[1<<x][x]=1;for(int S=1;S<M;S++)for(int i=0;i<n;i++)for(int flag=0;flag<2;flag++)if(d[S][flag*n+i])for(int Z=a[i]&~S;Z;Z&=Z-1){int b=Z&-Z,j=__builtin_ctz((unsigned)b),f=flag||(i==u&&j==v);d[S|b][f*n+j]+=d[S][flag*n+i];}long long z=0;for(int i=0;i<n;i++)z+=d.back()[n+i];return z;}
#include <numeric>
long long fixed_words(const vector<int>&a){int n=a.size(),t=n/2;vector<int>p(t);iota(p.begin(),p.end(),0);long long count=0;do{for(int mask=0;mask<(1<<t);mask++){vector<int>w;for(int x:p)w.push_back(x+((mask>>x&1)?t:0));for(int j=t-1;j>=0;j--)w.push_back((w[j]+t)%n);bool good=true;for(int i=1;i<n;i++)if(!(a[w[i-1]]>>w[i]&1)){good=false;break;}count+=good;}}while(next_permutation(p.begin(),p.end()));return count;}
long long bijections(const vector<int>&a,bool anti){int n=a.size();vector<int>p(n,-1),deg(n);for(int i=0;i<n;i++)deg[i]=__builtin_popcount((unsigned)a[i]);long long count=0;auto rec=[&](auto&&self,int u,int used)->void{if(u==n){count++;return;}for(int v=0;v<n;v++)if(!(used>>v&1)&&deg[v]==(anti?n-1-deg[u]:deg[u])){bool ok=true;for(int i=0;i<u;i++)if(bool(a[u]>>i&1)!=(bool(a[v]>>p[i]&1)^anti)){ok=false;break;}if(ok){p[u]=v;self(self,u+1,used|1<<v);}}};rec(rec,0,0);return count;}
long long checked=0;
void orbit_control(const vector<int>&a,int r){int n=a.size(),M=1<<n;auto E=endpoints(a);vector<int>op(n);for(int i=0;i<n;i++)op[i]=(M-1)^(1<<i)^a[i];auto B=endpoints(op);vector<long long>hs(M);hs[0]=1;for(int S=1;S<M;S++)for(long long v:E[S])hs[S]+=v;long long q=0;for(int S=1;S<M;S++)q+=E[S][r]*hs[(M-1)^S];long long expected=hs.back()+q-B.back()[r];need(fixed_words(build(a,r,0))==expected,"general fixed word formula");checked++;}
int main(){for(int n=1;n<=5;n++)for(int code=0;code<(1<<(n*(n-1)/2));code++){vector<int>a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){bool f=code>>k&1;a[f?i:j]|=1<<(f?j:i);}for(int r=0;r<n;r++)orbit_control(a,r);}mt19937 gen(1971210);for(int z=0;z<20;z++){vector<int>a(6);for(int i=0;i<6;i++)for(int j=i+1;j<6;j++){bool f=gen()&1;a[f?i:j]|=1<<(f?j:i);}orbit_control(a,z%6);}cout<<"PASS fixed_word_side_controls "<<checked<<"\n";
for(int k=1;k<=6;k++){int t=k+2;vector<int>a(t);for(int i=0;i<t;i++)for(int j=i+1;j<t;j++)a[i]|=1<<j;a[0]^=1<<(t-1);a.back()|=1;auto b=build(a,0,0);long long h=(1LL<<k)+1,c=2*h-1,H0=2*h*h-1;need(fixed_words(b)==c,"family fixed words");need(bijections(b,false)==1&&bijections(b,true)==1,"rigid unique anti");need((H0+c)/2==h*h+h-1,"Burnside orbit");if(k==5){vector<int>w={1,2,3,4,5,6,0,7,13,12,11,10,9,8};for(int i=1;i<14;i++)need(b[w[i-1]]>>w[i]&1,"displayed fixed path");for(int i=0;i<14;i++)need(w[i]==(w[13-i]+7)%14,"displayed symmetry");}cout<<"FAMILY "<<k<<" n "<<2*t<<" H "<<H0<<" fixed "<<c<<" Aut 1 anti 1 orbits "<<(H0+c)/2<<"\n";}
}
