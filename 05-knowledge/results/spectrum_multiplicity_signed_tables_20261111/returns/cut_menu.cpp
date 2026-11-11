// Exact one-/two-birth cut menu from endpoint and ordered-cover data.
#include <iostream>
#include <vector>
using namespace std;using I=long long;
vector<I> ends(const vector<int>&a){
 int n=a.size(),M=(1<<n)-1;vector<I>d((size_t)(M+1)*n);
 for(int i=0;i<n;++i)d[(size_t)(1<<i)*n+i]=1;
 for(int m=1;m<=M;++m)for(int i=0;i<n;++i)if(I x=d[(size_t)m*n+i]){
  int todo=a[i]&~m&M;
  while(todo){int b=todo&-todo;todo-=b;d[(size_t)(m|b)*n+__builtin_ctz((unsigned)b)]+=x;}
 }
 return d;
}
int main(){
 int mode,n,root;I low,high;cin>>mode>>n>>root>>low>>high;
 vector<int>a(n),rev(n);for(auto&x:a)cin>>x;
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(a[j]>>i&1)rev[i]|=1<<j;
 auto e=ends(a);int M=(1<<n)-1;vector<I>hs(M+1);hs[0]=1;
 for(int m=1;m<=M;++m)for(int i=0;i<n;++i)hs[m]+=e[(size_t)m*n+i];
 I H=hs[M];cout<<"{\"source_H\":"<<H<<",\"hits\":[";
 if(H>high){cout<<"],\"monotonicity_pruned\":true}\n";return 0;}
 auto s=ends(rev);vector<I>F(n),K(n*n);
 for(int m=1;m<=M;++m)for(int i=0;i<n;++i)if(I x=e[(size_t)m*n+i]){
  F[i]+=x*hs[M^m];
  if(m<M)for(int j=0;j<n;++j)K[i*n+j]+=x*s[(size_t)(M^m)*n+j];
 }
 int hits=0;
 for(int cut=0;cut<=M;++cut){
  I h=mode?H:0;
  for(int i=0;i<n;++i){
   bool in=cut>>i&1;
   if(mode){if(in){h+=e[(size_t)M*n+i]+F[i];for(int j=0;j<n;++j)if(!(cut>>j&1))h+=K[i*n+j];}}
   else if(in)h+=s[(size_t)M*n+i];
   else{h+=e[(size_t)M*n+i];for(int j=0;j<n;++j)if(cut>>j&1)h+=K[i*n+j];}
  }
  if(h==low||h==high){if(hits++)cout<<',';cout<<"{\"cut\":"<<cut<<",\"H\":"<<h<<",\"retains_indegree_one\":"<<((mode||!(cut>>root&1))?"true":"false")<<'}';}
 }
 cout<<"],\"monotonicity_pruned\":false}\n";
}
