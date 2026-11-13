#include <bits/stdc++.h>
using namespace std;
using ull=unsigned long long;
vector<int> R,P;
int n,N,M; vector<ull> dp; ull target;
ull H(const vector<int>&A){
 int n=A.size(),F=(1<<n)-1;
 vector<int> in(n);
 for(int j=0;j<n;j++) for(int i=0;i<n;i++) if(A[i]>>j&1) in[j]|=1<<i;
 vector<ull> sum(n+1);
 for(int S=1;S<=F;S++){
   ull z=0;
   for(int q=S;q;q&=q-1){
    int j=__builtin_ctz(q),U=S^(1<<j);
    ull h=0;
    if(!U)h=1;
    else for(int v=in[j]&U;v;v&=v-1) h+=dp[U*n+__builtin_ctz(v)];
    dp[S*n+j]=h;z+=h;
   }
   if(z>target) return z;
   if(S==F)return z;
 }
 return 0;
}
int main(){
 cin>>n>>target;N=n+2;M=(1<<N)-1;R.resize(n);P.resize(n);
 for(int&i:R)cin>>i;for(int&i:P)cin>>i;
 dp.resize(N*(1<<N));
 for(int w=0;w<(1<<(n+1));w++){
  vector<int>A(R);A.resize(N);
  for(int i=0;i<n;i++){
   if(w>>i&1){A[n]|=1<<i;A[P[i]]|=1<<(n+1);}
   else{A[i]|=1<<n;A[n+1]|=1<<P[i];}
  }
  if(w>>n&1)A[n]|=1<<(n+1);else A[n+1]|=1<<n;
  if(H(A)!=target)continue;
  auto good=[&](const vector<int>& A){
   int S=1;
   for(int k=0;k<N;k++)for(int i=0;i<N;i++)if(S>>i&1)S|=A[i];
   return S==M;
  };
  vector<int> D(N);
  for(int i=0;i<N;i++)for(int j=0;j<N;j++)if(A[i]>>j&1)D[j]|=1<<i;
  if(!good(A)||!good(D))continue;
  cout<<w<<" : ";for(int x:A)cout<<x<<" ";cout<<endl;
 }
}
