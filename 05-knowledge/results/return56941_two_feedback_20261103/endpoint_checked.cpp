// Previous packet's independent forward DP, with explicit integer-overflow guards.
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <limits>
#include <stdexcept>
using namespace std;using I=int64_t;
void add(I&x,I y){I z;if(__builtin_add_overflow(x,y,&z)||z<0)throw overflow_error("endpoint sum");x=z;}
I product(I x,I y){__int128 z=(__int128)x*y;if(z<0||z>numeric_limits<I>::max())throw overflow_error("endpoint product");return I(z);}
int main(){int n,r,b;while(cin>>n>>r>>b){if(n<1||n>20||r<0||b<0||r>=n||b>=n)throw invalid_argument("declared endpoint universe");vector<unsigned>T(n);for(auto&x:T)cin>>x;int full=(1<<n)-1;
 vector<I>e((full+1)*n),s((full+1)*n),h(full+1);h[0]=1;
 for(int i=0;i<n;i++)e[(1<<i)*n+i]=s[(1<<i)*n+i]=1;
 for(int m=1;m<=full;m++)for(unsigned z=m;z;z&=z-1){int v=countr_zero(z);I a=e[m*n+v],c=s[m*n+v];
  add(h[m],a);for(unsigned q=T[v]&~m;q;q&=q-1)add(e[(m|(q&-q))*n+countr_zero(q)],a);
  for(unsigned q=(full^T[v]^(1<<v))&~m;q;q&=q-1)add(s[(m|(q&-q))*n+countr_zero(q)],c);}
 vector<int>p={r,b};I P=0;vector<I>E(2),S(2),W(4);
 for(int m=0;m<=full;m++){add(P,product(h[m],h[full^m]));for(int i=0;i<2;i++){add(E[i],product(e[m*n+p[i]],h[full^m]));add(S[i],product(s[m*n+p[i]],h[full^m]));for(int j=0;j<2;j++)add(W[2*i+j],product(e[m*n+p[i]],s[(full^m)*n+p[j]]));}}
 cout<<h[full]<<' '<<P;for(int i=0;i<2;i++)cout<<' '<<e[full*n+p[i]]<<' '<<s[full*n+p[i]]<<' '<<h[full^(1<<p[i])]<<' '<<E[i]<<' '<<S[i];for(I x:W)cout<<' '<<x;cout<<'\n';
}}
