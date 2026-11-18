#include <vector>
#include <iostream>
#include <array>
using namespace std;
using z=unsigned long long;
array<z,4> reads(const vector<unsigned>& a,int r){
 int n=a.size();unsigned f=(1u<<n)-1;
 vector<z>D((f+1)*n),H(f+1),S(f+1);vector<unsigned>rev(n);
 for(int j=0;j<n;j++)rev[j]=f^(1u<<j)^a[j],D[(1u<<j)*n+j]=1;
 H[0]=1;
 auto extend=[&](const vector<unsigned>&A, bool h){
  for(unsigned s=1;s<=f;s++)for(unsigned u=s;u;u&=u-1){
    int j=__builtin_ctz(u);auto c=D[s*n+j];if(!c)continue;
    if(h)H[s]+=c;
    for(unsigned w=A[j]&~s;w;w&=w-1)D[(s|(w&-w))*n+__builtin_ctz(w)]+=c;
  }
 };
 extend(a,true); auto e=D[f*n+r];
 fill(D.begin(),D.end(),0);
 for(int j=0;j<n;j++)D[(1u<<j)*n+j]=1;
 extend(rev,false);z w=0;
 for(unsigned s=1;s<=f;s++)if(s>>r&1)w+=D[s*n+r]*H[f^s];
 return {H[f],e,w,H[f^(1<<r)]};
}
int main(){int n,r;while(cin>>n>>r){vector<unsigned>a(n);for(auto&j:a)cin>>j;
for(auto j:reads(a,r))cout<<j<<" ";cout<<endl;}}
