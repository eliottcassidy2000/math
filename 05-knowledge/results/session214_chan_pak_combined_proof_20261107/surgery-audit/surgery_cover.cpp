#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using I=int64_t;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
struct Polynomial {I constant;vector<I>linear;vector<vector<I>>quad;};
Polynomial cover_poly(const vector<int>&a){int n=a.size(),N=1<<n,F=N-1;vector<vector<I>> E(N,vector<I>(n)),S=E;for(int v=0;v<n;v++)E[1<<v][v]=S[1<<v][v]=1;
 for(int mask=1;mask<N;mask++)for(int v=0;v<n;v++)if((mask>>v)&1){int rest=mask^(1<<v);if(!rest)continue;for(int u=0;u<n;u++)if((rest>>u)&1){if((a[u]>>v)&1)E[mask][v]+=E[rest][u];if((a[v]>>u)&1)S[mask][v]+=S[rest][u];}}
 vector<vector<I>> W(n,vector<I>(n));for(int mask=1;mask<F;mask++){int rest=F^mask;for(int i=0;i<n;i++)if(E[mask][i])for(int j=0;j<n;j++)if(S[rest][j])W[i][j]+=E[mask][i]*S[rest][j];}
 Polynomial p{0,vector<I>(n),vector<vector<I>>(n,vector<I>(n))};for(int i=0;i<n;i++){p.constant+=E[F][i];p.linear[i]=S[F][i]-E[F][i];for(int j=0;j<n;j++){p.linear[i]+=W[j][i];p.quad[i][j]=-W[i][j]-W[j][i];}}
 return p;
}
I direct(const vector<int>&a){int n=a.size(),N=1<<n;vector<vector<I>>d(N,vector<I>(n));for(int v=0;v<n;v++)d[1<<v][v]=1;for(int mask=1;mask<N;mask++)for(int v=0;v<n;v++)if(d[mask][v]){int bits=a[v]&~mask;while(bits){int b=bits&-bits;bits-=b;d[mask|b][__builtin_ctz(b)]+=d[mask][v];}}I z=0;for(I x:d[N-1])z+=x;return z;}
vector<int>parent(int L,int r){vector<int>a(L+1);for(int i=0;i<L;i++){for(int j=i+1;j<L;j++)a[i]|=1<<j;if(i!=r)a[i]|=1<<L;}if(r>=0)a[L]|=1<<r;return a;}
I eval(const Polynomial&p,int mask){I z=p.constant;for(int i=0;i<(int)p.linear.size();i++)if((mask>>i)&1){z+=p.linear[i];for(int j=0;j<i;j++)if((mask>>j)&1)z+=p.quad[i][j];}return z;}
int main(){const int target=7477;set<int>attained;I cases=0,controls=0,apex_cases=0;vector<pair<int,I>>by_order;
 for(int L=1;L<=15;L++){I start=cases;for(int r=0;r<L;r++){auto a=parent(L,r);auto p=cover_poly(a);for(int ba=0;ba<2;ba++){I h=p.constant+ba*p.linear[L];vector<I>d(L);for(int i=0;i<L;i++)d[i]=p.linear[i]+ba*p.quad[i][L];int old=0;for(int step=0;step<(1<<L);step++){int mask=step^(step>>1);if(step){int i=__builtin_ctz(mask^old),sign=((mask>>i)&1)?1:-1;h+=sign*d[i];for(int j=0;j<L;j++)if(j!=i)d[j]+=sign*p.quad[i][j];}need(h>0&&h%2==1,"positive odd count");need(h!=target,"found forbidden value");if(h<target)attained.insert(h);if(L<=6){vector<int>b=a;for(int i=0;i<=L;i++)if(!((mask|(ba<<L))>>i&1))b[i]|=1<<(L+1);b.push_back(mask|(ba<<L));need(h==direct(b),"direct DP control");controls++;}cases++;old=mask;}}}by_order.push_back({L,cases-start});}
 // With only one apex, enumerate all binary cuts on a chain of length <=16.
 for(int L=0;L<=16;L++){for(int b=0;b<(1<<L);b++){I N=0;for(int i=0;i<L;i++)if((b>>i)&1)for(int j=i+1;j<L;j++)if(!((b>>j)&1))N+=(j-i==1)?1:(1LL<<(j-i-2));I h=1+2*N;need(h!=target,"one-apex forbidden value");if(h<target)attained.insert(h);apex_cases++;if(L<=8){vector<int>a(L+1);for(int i=0;i<L;i++){for(int j=i+1;j<L;j++)a[i]|=1<<j;if(!((b>>i)&1))a[i]|=1<<L;}a[L]=b;need(h==direct(a),"one-apex DP control");controls++;}}}
 vector<int>positive={510,1020,1016,1520,2016,1984,896,768,1536,1033,199};need(direct(positive)==7477,"two-exit positive control");
 cout<<"{\"status\":\"PASS\",\"target\":7477,\"two_apex_case_count\":"<<cases<<",\"one_apex_case_count\":"<<apex_cases<<",\"direct_DP_controls\":"<<controls<<",\"two_exit_positive_control_H\":7477,\"per_chain_length\":[";bool first=true;for(auto[L,c]:by_order){if(!first)cout<<",";first=false;cout<<"["<<L<<","<<c<<"]";}cout<<"],\"missing_odds_below7477\":[";first=true;for(int h=1;h<target;h+=2)if(!attained.count(h)){if(!first)cout<<",";first=false;cout<<h;}cout<<"]}\n";
}
