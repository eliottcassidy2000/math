// Independent backward recurrences; input n,k,rows,ports. Checked uint64.
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace std; using I=uint64_t;
void add(I& a,I b){if(b>numeric_limits<I>::max()-a)throw overflow_error("sum");a+=b;}
I mul(I a,I b){if(a && b>numeric_limits<I>::max()/a)throw overflow_error("product");return a*b;}
int main(){int n,k;while(cin>>n>>k){if(n<1||n>20||k<1||k>n)throw invalid_argument("n/ports");
 vector<unsigned>T(n),in(n);for(auto &a:T)cin>>a;vector<int>p(k);for(auto&a:p){cin>>a;if(a<0||a>=n)throw invalid_argument("port");}
 int f=(1<<n)-1;for(int i=0;i<n;i++){if(T[i]&~unsigned(f) || T[i]>>i&1)throw invalid_argument("row");for(int j=0;j<n;j++)if(i!=j){if(((T[i]>>j)&1)==((T[j]>>i)&1))throw invalid_argument("tournament");if(T[j]>>i&1)in[i]|=1<<j;}}
 vector<I>e(size_t(f+1)*n),s(size_t(f+1)*n),h(f+1);h[0]=1;
 for(int m=1;m<=f;m++)for(unsigned q=m;q;q&=q-1){int v=countr_zero(q),z=m^(1<<v);I &ev=e[size_t(m)*n+v],&sv=s[size_t(m)*n+v];
  if(!z)ev=sv=1;else{for(unsigned b=z&in[v];b;b&=b-1)add(ev,e[size_t(z)*n+countr_zero(b)]);for(unsigned b=z&T[v];b;b&=b-1)add(sv,s[size_t(z)*n+countr_zero(b)]);}add(h[m],ev);}
 I P=0;vector<I>E(k),S(k),W(k*k);
 for(int m=0;m<=f;m++){int c=f^m;add(P,mul(h[m],h[c]));for(int i=0;i<k;i++){add(E[i],mul(e[size_t(m)*n+p[i]],h[c]));add(S[i],mul(s[size_t(m)*n+p[i]],h[c]));for(int j=0;j<k;j++)add(W[i*k+j],mul(e[size_t(m)*n+p[i]],s[size_t(c)*n+p[j]]));}}
 cout<<h[f]<<' '<<P;for(int i=0;i<k;i++)cout<<' '<<e[size_t(f)*n+p[i]]<<' '<<s[size_t(f)*n+p[i]]<<' '<<h[f^(1<<p[i])]<<' '<<E[i]<<' '<<S[i];for(I a:W)cout<<' '<<a;cout<<'\n';
}}
