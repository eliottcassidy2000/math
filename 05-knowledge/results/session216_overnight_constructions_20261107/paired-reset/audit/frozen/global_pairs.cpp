#include<bits/stdc++.h>
using namespace std;using I=long long;struct Rec{int n,ix;vector<int>a;};
int main(){map<pair<I,I>,array<Rec,2>>seen;long long graphs=0,hits=0;ofstream all("GLOBAL_COUNTS.tsv");
for(int n=1;n<=9;n++){ifstream in("../session216-carriers-order9/BANK"+to_string(n)+".txt");int total;in>>total;for(int ix=0;ix<total;ix++){int nn;in>>nn;vector<int>a(n);for(auto&x:a)in>>x;int N=1<<n,F=N-1;vector<I>dp(N*n),H(N);H[0]=1;for(int v=0;v<n;v++)dp[(1<<v)*n+v]=1;for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if((A>>w&1)&&(a[w]>>v&1))dp[S*n+v]+=dp[A*n+w];H[S]+=dp[S*n+v];}I C=0;for(int A=0;A<N;A++)C+=H[A]*H[F^A];all<<n<<' '<<ix<<' '<<H[F]<<' '<<C<<'\n';graphs++;auto k=make_pair(H[F],C);auto&z=seen[k];z[n%2]={n,ix,a};if(z[0].n&&z[1].n){hits++;if(hits<=20){cout<<"MATCH "<<H[F]<<" "<<C<<"\n";for(auto&r:z){cout<<r.n<<" "<<r.ix;for(auto v:r.a)cout<<" "<<v;cout<<"\n";}}}}cerr<<"order "<<n<<" graphs "<<graphs<<" hits "<<hits<<"\n";}
}
