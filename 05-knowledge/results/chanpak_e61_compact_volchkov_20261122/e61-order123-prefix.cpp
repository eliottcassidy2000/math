#include<bits/stdc++.h>
using namespace std;
using I=__int128_t;
long long P[65536],B[65536];
I p(I z){
 I V=0;
 while(z>=65536){
  int r=z&65535;z>>=16;
  V+=(65535-B[r])*z+P[r]+(r<32768&&(z&1));
 }
 return V+P[int(z)];
}
I mod(I n,I q){I r=n%q;return r<0?r+q:r;}
I floorq(I a,I b){return (a-mod(a,b))/b;}
I ceilq(I a,I b){return -floorq(-a,b);}
string out(I n){if(n<0)return "-"+out(-n);string r;do{r+='0'+n%10;n/=10;}while(n);reverse(r.begin(),r.end());return r;}
I rev(I x,int k){
 I out=0;int j=0;
 for(;j<k;j+=16){out=(out<<16)|B[(int)(x&65535)];x>>=16;}
 return out>>(j-k);
}
int main(int argc,char**argv){
 if(argc!=2 || atoi(argv[1])<1 || atoi(argv[1])>110) return 2;
 for(int i=0;i<65536;i++){
  P[i]=P[i>>1]+(i&1?0:(i>>1)+((i>>1)&1));
  for(int j=0;j<16;j++)B[i]|=((i>>j)&1)<<(15-j);
 }
 I N=(I(1)<<60)*((I(1)<<61)-1)-1,X=I(1)<<atoi(argv[1]);
 const array<I,5> charts[]={{81,25,43,68,115},{1021,125,915,940,112}};
 for(auto [c,h,k0,k1,k]:charts){
  int w=k-1-atoi(argv[1]);I A=I(1)<<w;
  I J=c*A*X+h*(A-1)*X+h*p(X);
  I Lo=max<I>(0,ceilq(J+k0-N,h*X)),Hi=min(A-1,floorq(J+k1+(c+h)*(A-1)+h-N,h*X));
  I eligible=0,hits=0;
  for(I q=Lo;q<=Hi;q++){
   I r=rev(q,w),y=A*X+r,K=(r&1)?k1:k0;
   if(mod(K+c*y-N,2*h))continue;
   ++eligible;
   I f=p(X)+(A-1-q)*X+p(r)+(r<A/2&&(X&1));
   if(K+c*y+h*f==N){++hits;cout<<"HIT "<<out(y)<<endl;}
  }
  cout<<"chart "<<out(c)<<" prefix bits "<<atoi(argv[1])+1<<" reverse window "<<out(max<I>(0,Hi-Lo+1))<<" congruence "<<out(eligible)<<" hits "<<out(hits)<<endl;
 }
}
