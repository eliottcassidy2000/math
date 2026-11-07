#include <bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,9>;
P direct(const vector<int>&a){int n=a.size(),N=1<<n,F=N-1;vector<I>en(N*n),st(N*n),H(N),C(N);H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;
for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}
for(int S=0;S<N;S++){int A=S;while(true){C[S]+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}}
I L=0,R=0,T=0;for(int S=1;S<N;S++){L+=en[S*n]*H[F^S];R+=st[S*n]*H[F^S];T+=st[S*n]*C[F^S];}return{H[F],st[F*n],en[F*n],L,R,H[F^1],C[F],C[F^1],T};}
int main(){ifstream in("BANK9.txt");int total;in>>total;vector<vector<int>>bank;for(int i=0;i<total;i++){int n;in>>n;vector<int>a(n);for(auto&v:a)in>>v;bank.push_back(a);}ifstream cand("CANDIDATES10.txt");int ix,mask,checks=0;P expect;while(cand>>ix>>mask){for(auto&x:expect)cand>>x;auto&a=bank[ix];vector<int>g{mask<<1};for(int i=0;i<(int)a.size();i++)g.push_back((a[i]<<1)|(1-((mask>>i)&1)));if(direct(g)!=expect)throw runtime_error("literal profile mismatch");checks++;}cout<<"{\"status\":\"PASS\",\"literal_root10_profiles\":"<<checks<<"}\n";}
