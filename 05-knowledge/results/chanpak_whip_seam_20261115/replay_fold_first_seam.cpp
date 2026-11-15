#include <bits/stdc++.h>
using namespace std;using L=long long;
struct Data{int n,F;vector<L> h,s;Data(const vector<int>&R):n(R.size()),F((1<<n)-1),h(F+1),s(n){
 vector<L> p((F+1)*n,0);h[0]=1;
 for(int m=1;m<=F;m++)for(int x=m;x;x&=x-1){
  int i=__builtin_ctz(unsigned(x)),t=m^(1<<i);
  L z=0;
  if(!t)z=1;
  else for(int y=t&R[i];y;y&=y-1)z+=p[t*n+__builtin_ctz(unsigned(y))];
  p[m*n+i]=z;h[m]+=z;
 }
 for(int m=1;m<=F;m++)for(int x=m;x;x&=x-1){int i=__builtin_ctz(unsigned(x));s[i]+=p[m*n+i]*h[F^m];}
}};
vector<int> add(vector<int> r,int u){
 int n=r.size();for(int i=0;i<n;i++)if(!(u>>i&1))r[i]|=1<<n;r.push_back(u);return r;
}
int main(){
 vector<int> R={48,29,25,1,8,30},T={18,20,17,7,8};
 vector<vector<int>> bases={R,R,T,R};int mm[]={25,16,1,1}, dd[]={1,2,2,2};
 using S=tuple<int,int,int,int>;
 for(int j=3;j<=5;j++){
  map<pair<L,L>,map<L,vector<S>>> bank;
  int cnt=0;
  for(int i=0;i<4;i++){
   auto r=bases[i];int m=mm[i];
   for(int u=0;u<dd[i]+j-3;u++){m|=1<<r.size();r=add(r,0);}
   r=add(r,0);int n=r.size(); // old guard mark 0
   for(int e=0;e<(1<<n);e++){
    auto a=add(r,e);Data P(a);
    for(int b=0;b<=1;b++){
     int m2=m|(b<<n);L S0=0;
     for(int k=0;k<=n;k++)if(m2>>k&1)S0+=P.s[k];
     Data Q(add(a,m2));++cnt;
     bank[{P.h.back(),S0}][Q.h.back()].emplace_back(i,e,b,n);
    }
   }
  }
  int z=0,r=0;
  for(auto&[I,A]:bank)for(auto&[x,L]:A)if(!A.count(x-2)){
   int w=0,B=0;
   while(A.count(x+2*w)){for(auto[i,e,b,n]:A[x+2*w])B|=1<<b;w++;}
   if(w>z|| (B==3&&w>r)){
    cout<<"j "<<j<<" run "<<w<<" q "<<B<<" profile "<<I.first<<" "<<I.second<<" ";
    for(int k=0;k<w;k++){cout<<x+2*k<<":";for(auto[i,e,b,n]:A[x+2*k])cout<<i<<","<<e<<","<<b<<"|";cout<<";";}
    cout<<endl;
   }
   z=max(z,w);if(B==3)r=max(r,w);
  }
  cout<<"complete j="<<j<<" candidates="<<cnt<<" states="<<bank.size()<<" max="<<z<<" mixed="<<r<<endl;
 }
}
