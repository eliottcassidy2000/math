#include<bits/stdc++.h>
using namespace std;
using V=vector<int>;using I=long long;
constexpr int N=14,M=1<<N;
V app(V a,int q){int n=a.size();for(int i=0;i<n;i++)if(!(q>>i&1))a[i]|=1<<n;a.push_back(q);return a;}
struct Box {
 int n,U;
 vector<I> h,p;
 I D[2][N][N]={};
 Box(V a):n(a.size()),U(1<<n),h(U),p(U){
  static I d[M][N];
  auto grow=[&](){for(int s=1;s<U;s++)for(int b=s;b;b&=b-1){
   int i=__builtin_ctz(b);I v=d[s][i];if(!v)continue;
   for(int t=a[i]&~s;t;t&=t-1)d[s|(t&-t)][__builtin_ctz(t)]+=v;
  }};
  memset(d,0,U*sizeof(d[0]));for(int i=0;i<n;i++)d[1<<i][i]=1;grow();h[0]=1;
  for(int s=1;s<U;s++)h[s]=accumulate(d[s],d[s]+n,0ll);
  for(int s=0;s<U;s++)for(int t=s;;t=(t-1)&s){p[s]+=h[t]*h[s^t];if(!t)break;}
  for(int i=0;i<n;i++){
   memset(d,0,U*sizeof(d[0]));d[1<<i][i]=1;grow();
   for(int s=0;s<U;s++)if((s>>i&1)&&!(__builtin_popcount(s)&1))for(int j=0;j<n;j++)if(d[s][j]){
    D[0][i][j]+=2*d[s][j]*h[U-1-s];
    D[1][i][j]+=4*d[s][j]*p[U-1-s];
   }
  }
 }
};
int main(){
 V a0={174,124,120,240,65,80,1,118};V b0=a0;b0[5]^=64;b0[6]^=32;
 for(auto [ca,cb]:vector<pair<int,int>>{{43,75}}){
  auto a=app(app(app(app(a0,ca),cb),612),392),b=app(app(app(app(b0,ca),cb),612),392);
  for(int q=0;q<4096;q++){
   Box A(app(a,q)),B(app(b,q));
   int n=A.n,C=1<<(n-1);
   static I p[2][M], E[2][N][N];memset(p,0,sizeof(p));memset(E,0,sizeof(E));
   for(int c=0;c<2;c++){
    I sums[N]={};
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
     int d=B.D[c][i][j]-A.D[c][i][j];sums[i]+=d;
     E[c][i][j]+=d;E[c][j][i]+=d;
    }
    for(int s=1;s<C;s++){
     int u=s&-s,i=__builtin_ctz(u),t=s^u;I g=p[c][t]+sums[i];
     for(int k=t;k;k&=k-1)g-=E[c][i][__builtin_ctz(k)];
     p[c][s]=g;
    }
   }
   for(int s=0;s<C;s++)if(B.h.back()-A.h.back()+p[0][s]==10 && 2*(B.p.back()-A.p.back())+p[1][s]==-8){
    I h=A.h.back(),P=2*A.p.back();
    for(int i=0;i<n;i++)if(s>>i&1)for(int j=0;j<n;j++)if(!(s>>j&1)){h+=A.D[0][i][j];P+=A.D[1][i][j];}
    cout<<ca<<" "<<cb<<" "<<q<<" "<<s<<" "<<h<<" "<<P<<endl;
   }
  }
 }
}
