// Exact quadratic one-newborn-cut search on a named source.
// Input: n root target, followed by n outgoing masks.
#include <iostream>
#include <vector>
using namespace std;using I=long long;
vector<I> ends(const vector<int>& rows){
 int n=rows.size(),M=(1<<n)-1;vector<I>d((M+1)*n);
 for(int i=0;i<n;++i)d[(1<<i)*n+i]=1;
 for(int m=1;m<=M;++m)for(int i=0;i<n;++i)if(I x=d[m*n+i]){
  int a=rows[i]&~m&M;
  while(a){int b=a&-a;a-=b;int j=__builtin_ctz((unsigned)b);d[(m|b)*n+j]+=x;}
 }
 return d;
}
int main(){
 int n,r;I target;cin>>n>>r>>target;vector<int>a(n),rev(n);
 for(auto&x:a)cin>>x;
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(a[j]>>i&1)rev[i]|=1<<j;
 auto e=ends(a),s=ends(rev);int M=(1<<n)-1;vector<I>W(n*n);I H=0;for(int i=0;i<n;++i)H+=e[M*n+i];
 for(int mask=1;mask<M;++mask){int other=M^mask;
  for(int i=0;i<n;++i)if(I x=e[mask*n+i])for(int j=0;j<n;++j)if(I y=s[other*n+j])W[i*n+j]+=x*y;
 }
 cout<<"{\"source_H\":"<<H<<",\"target\":"<<target<<",\"cuts\":"<<(M+1)<<",\"hits\":[";int hits=0;
 I closest=1LL<<60,nearH=0;int nearCut=0;
 for(int cut=0;cut<=M;++cut){I h=0;
  for(int i=0;i<n;++i){if(cut>>i&1)h+=s[M*n+i];else{h+=e[M*n+i];for(int j=0;j<n;++j)if(cut>>j&1)h+=W[i*n+j];}}
  if(llabs(h-target)<closest){closest=llabs(h-target);nearH=h;nearCut=cut;}
  if(h==target){if(hits++)cout<<',';cout<<"{\"cut\":"<<cut<<",\"preserves_root_indegree_one\":"<<((cut>>r&1)?"false":"true")<<'}';}
 }
 cout<<"],\"nearest_H\":"<<nearH<<",\"nearest_cut\":"<<nearCut<<",\"nearest_distance\":"<<closest<<"}\n";
}
