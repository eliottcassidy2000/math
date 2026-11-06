#include <bits/stdc++.h>
using namespace std;
using I=long long; using G=vector<int>;
void need(bool q,const string&s){if(!q)throw runtime_error(s);}
struct Counts {
 G a; int n,N,F; vector<I> en,st,H;
 Counts(G g):a(g),n(g.size()),N(1<<n),F(N-1),en(N*n),st(N*n),H(N){
  for(int i=0;i<n;i++){need(a[i]>=0&&a[i]<N&&!(a[i]>>i&1),"row");for(int j=0;j<i;j++)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"orientation");}
  H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;
  for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}
 }
 I h()const{return H[F];}
 I C(int S)const{I z=0;int A=S;while(true){z+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}return z;}
 I L(int v,int S)const{I z=0;for(int A=S;;A=(A-1)&S){z+=en[A*n+v]*H[S^A];if(!A)break;}return z;}
 I R(int v,int S)const{I z=0;for(int A=S;;A=(A-1)&S){z+=st[A*n+v]*H[S^A];if(!A)break;}return z;}
 array<I,9> profile(int v)const{int A=F^(1<<v);I h=H[F],s=st[F*n+v],e=en[F*n+v],l=L(v,F),r=R(v,F),d=H[A],b=C(A),B=h+l-s,D=r+b-e-B;return {h,s,e,l,r,d,b,B,D};}
};
G graph(int n,int mask){G a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){if(mask>>k&1)a[i]|=1<<j;else a[j]|=1<<i;}return a;}
G cyclic(const G&t){int m=t.size();G a(m+2);a[0]=2;for(int i=0;i<m;i++){a[1]|=1<<(i+2);a[i+2]=(t[i]<<2)|1;}return a;}
G zmove(G a,int v){int n=a.size();for(int i=0;i<n;i++)if(i!=v)a[i]|=1<<n;a.push_back(1<<v);return a;}
bool strong(G a){int n=a.size(),F=(1<<n)-1;for(int i=0;i<n;i++)a[i]|=1<<i;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(a[i]>>k&1)a[i]|=a[k];return all_of(a.begin(),a.end(),[&](int x){return x==F;});}
