#include<bits/stdc++.h>
using namespace std;
using V=vector<int>;
struct Data{
 int n,F;V P,h; Data(const V&a):n(a.size()),F((1<<n)-1),P((F+1)*n),h(F+1){
 h[0]=1;
 for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){
   int i=__builtin_ctz(b),R=S^(1<<i),v=!R;
   for(int x=R&a[i];x;x&=x-1)v+=P[R*n+__builtin_ctz(x)];
   P[S*n+i]=v;h[S]+=v;
 }}
};
V join(V a,const V&b){
 int n=a.size(),added=((1<<b.size())-1)<<n;
 for(auto&x:a)x|=added;
 for(int x:b)a.push_back(x<<n);
 return a;
}
int odd(int a){return a>>__builtin_ctz(a);}
int sums(const Data &d){int p=0;for(int i=0;i<=d.F;i++)p+=d.h[i]*d.h[d.F^i];return p;}
map<int,int> minimum,matches;
void examine(V a,Data &d){

 int n=d.n,F=d.F, h=d.h[F];V X(F+1),E(F+1),w(n),e(n),ev(n*n);X[0]=h;
 map<pair<int,int>,map<int,int>> curves;
 for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){int i=__builtin_ctz(b);int t=d.P[S*n+i]*d.h[F^S];w[i]+=t;if(__builtin_popcount(S)&1)e[i]+=t;}
 V p((F+1)*n);
 for(int st=0;st<n;st++){
  fill(p.begin(),p.end(),0);p[(1<<st)*n+st]=1;
  for(int S=1<<st;S<=F;S++)if(S>>st&1)for(int b=S;b;b&=b-1){
   int k=__builtin_ctz(b),v=p[S*n+k];if(!v)continue;
   if((__builtin_popcount(S)&1)==0)ev[st*n+k]+=v*d.h[F^S];
   for(int u=(F^S)&a[k];u;u&=u-1)p[(S|(u&-u))*n+__builtin_ctz(u)]+=v;
  }
 }
 V cut(F+1),tot(n);
 for(int i=0;i<n;i++)for(int k=0;k<n;k++)tot[i]+=ev[i*n+k];
 curves[{h,0}][h]=0;
 for(int q=1;q<=F;q++){
   int i=__builtin_ctz(q),s=q&(q-1);X[q]=X[s]+w[i];E[q]=E[s]+e[i];cut[q]=cut[s]+tot[i];
   for(int b=s;b;b&=b-1){int k=__builtin_ctz(b);cut[q]-=ev[i*n+k]+ev[k*n+i];}
   int H=h+2*cut[q];
   curves[{X[q],E[q]}].try_emplace(H,q);
 }
 for(auto& [p,m]:curves)for(auto& [H,q]:m)if(m.count(H+2)&&m.count(H+4)){
   int &best=minimum[n];
   if(!best||H<best){
    best=H;
    cout<<"window "<<n<<" "<<H<<" E "<<p.second<<" O "<<p.first-p.second<<" g "<<h<<" Q "<<q<<" "<<m[H+2]<<" "<<m[H+4]<<" rows";
    for(int x:a)cout<<" "<<x;cout<<"\n";
   }
   ++matches[n];
 }
}
int main(int argc,char **argv){
 int cap=argc>1?atoi(argv[1]):10;
 for(int n=6;n<=cap;n++){
  int k=n<9?2:3,t=n-3*k,bits=3*k*(k-1)/2+k*t+t*(t-1)/2;
  for(int mask=0;mask<(1<<bits);mask++){
   V a(n);int at=0;
   auto val=[&](){return(mask>>at++)&1;};
   auto put=[&](int i,int j,bool b){if(b)a[i]|=1<<j;else a[j]|=1<<i;};
   for(int i=0;i<k;i++)for(int z=0;z<3;z++)put(3*i+z,3*i+(z+1)%3,1);
   for(int i=0;i<k;i++)for(int j=i+1;j<k;j++)for(int d=0;d<3;d++){
     bool b=val();for(int z=0;z<3;z++)put(3*i+z,3*j+(z+d)%3,b);
   }
   for(int i=0;i<k;i++)for(int j=3*k;j<n;j++){
    bool b=val();for(int z=0;z<3;z++)put(3*i+z,j,b);
   }
   for(int i=3*k;i<n;i++)for(int j=i+1;j<n;j++)put(i,j,val());
   Data d(a);examine(a,d);
  }
  cout<<"done "<<n<<" orientations "<<(1<<bits)<<" matches "<<matches[n]<<" minimum "<<minimum[n]<<endl;
 }
}
