#include<bits/stdc++.h>
using namespace std;
using V=vector<int>;using I=long long;
pair<I,I> prof(const V&a){
 int n=a.size(),M=1<<n;
 vector<array<I,20>>d(M);vector<I> h(M);h[0]=1;
 for(int j=0;j<n;j++)d[1<<j][j]=1;
 for(int s=1;s<M;s++)for(int b=s;b;b&=b-1){
  int i=__builtin_ctz(b);I v=d[s][i];if(!v)continue;h[s]+=v;
  for(int p=a[i]&~s;p;p&=p-1)d[s|(p&-p)][__builtin_ctz(p)]+=v;
 }
 I P=0;for(int s=0;s<M;s++)P+=h[s]*h[M-1-s];return {h.back(),P};
}
V app(V a,int x){int n=a.size();for(int i=0;i<n;i++)if(!(x>>i&1))a[i]|=1<<n;a.push_back(x);return a;}
int main(){
 V a={174,124,120,240,65,80,1,118},b={174,124,120,240,65,16,33,118};
 int count=0;
 for(int q=0;q<256;q++){
  auto x=app(a,q),y=app(b,q);
  for(int s=0;s<256;s++){
   auto xx=app(x,s),yy=app(y,s);auto [h,p]=prof(xx);auto [g,Q]=prof(yy);
   if(g-h==10&&Q-p==-8){count++;cout<<q<<" "<<s<<" "<<h<<" "<<p<<endl;}
  }
 }
 cerr<<count<<endl;
}
