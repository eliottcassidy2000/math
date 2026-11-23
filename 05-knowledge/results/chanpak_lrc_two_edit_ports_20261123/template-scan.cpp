#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
using V=vector<int>;using Trip=array<int,3>;
int height(const V& t){
 const int m=t.size(),n=m+2,N=1<<n,a=m,b=m+1;
 uint32_t r[16]={};
 for(int i=0;i<m;i++){
  r[i]=((1<<m)-1)^((1<<(i+1))-1);
  for(int ap=0;ap<2;ap++)
   if(t[i]>>ap&1)r[m+ap]|=1<<i;else r[i]|=1<<(m+ap);
 }
 r[b]|=1<<a;
 static int dp[(1<<12)*12]; fill(dp,dp+N*n,0);
 for(int i=0;i<n;i++)dp[(1<<i)*n+i]=1;
 for(int s=1;s<N-1;s++)for(int z=s;z;z&=z-1){
  int i=__builtin_ctz(z),v=dp[s*n+i];if(v)
   for(int w=r[i]&~s;w;w&=w-1){int j=__builtin_ctz(w);dp[(s|1<<j)*n+j]+=v;}
 }
 int h=0; for(int j=0;j<n;j++)h+=dp[(N-1)*n+j];return h;
}
Trip coeff(int x,int y,int z){
 assert((z-3*y+2*x)%6==0);int C=(z-3*y+2*x)/6;
 assert((y-x-6*C)%2==0);int B=(y-x-6*C)/2;
 return {x-2*B-3*C,B,C};
}
int main(int argc,char**argv){
 const int cap=argc>1?stoi(argv[1]):6;
 if(cap<0||cap>7)return 2;
 for(int k=0;k<=cap;k++){
  int m=k+1,L=1<<(2*m);
  vector<int> h(L);V t(m);
  for(int w=0;w<L;w++){
   for(int j=0;j<m;j++)t[j]=w>>(2*j)&3;
   h[w]=height(t);
  }
  using Record=array<Trip,4>;
  map<Trip,int> lower,certain;
  map<Trip,array<int,4>> rep;
  map<Trip,Record> data;
  int flagged=0;uint64_t universe=0;
  for(int w=0;w<L;w++)for(int gap=0;gap<m;gap++){
   for(int i=0;i<2*m;i++)if(i/2!=gap)
    for(int j=i+1;j<2*m;j++)if(j/2!=gap){
     universe++;
     const array<int,4> ww={w,w^(1<<i),w^(1<<j),w^(1<<i)^(1<<j)};
     const int top=7*h[w];
     if(h[ww[3]]!=top||h[ww[1]]<=h[w]||h[ww[2]]<=h[w]||h[ww[1]]>=top||h[ww[2]]>=top)continue;
     auto H=[&](int word,int width){
      V v;for(int q=0;q<m;q++)v.insert(v.end(),q==gap?width:1,word>>(2*q)&3);
      return height(v);
     };
     const int p2=H(w,2);if(H(ww[3],2)!=7*p2)continue;
     const int p3=H(w,3);if(H(ww[3],3)!=7*p3)continue;
     auto C=coeff(h[w],p2,p3);
     auto B=coeff(h[ww[1]],H(ww[1],2),H(ww[1],3));
     auto D=coeff(h[ww[2]],H(ww[2],2),H(ww[2],3));
     Trip upper={7*C[0],7*C[1],7*C[2]};
     Record U={C,B,D,upper};
     auto positive=[](Trip low,Trip high){
      Trip d;for(int z=0;z<3;z++)d[z]=high[z]-low[z];
      // Express in nonempty-slot readings D, Q-2D, R-3Q+3D.
      return d[0]+2*d[1]+3*d[2]>0 && d[1]+3*d[2]>=0 && d[2]>=0;
     };
     lower[C]++;
     if(positive(C,B)&&positive(C,D)&&positive(B,upper)&&positive(D,upper)){
      certain[C]++;
      if(!rep.count(C)){rep[C]={w,gap,i,j};data[C]=U;}
     } else {
      flagged++;
      if(!rep.count(C)){rep[C]={w,gap,i,j};data[C]=U;}
     }
    }
  }
  cout<<"k="<<k<<" universe="<<universe<<" undecided="<<flagged<<endl;
  for(auto [abc,count]:lower){
   auto [w,g,i,j]=rep.at(abc);
   cout<<"  "<<abc[0]<<","<<abc[1]<<","<<abc[2]<<" raw="<<count<<" proven="<<certain[abc]<<" rep ";
   for(int z=0;z<m;z++)cout<<(w>>(2*z)&3);
   cout<<" gap,contact-bits="<<g<<","<<i<<","<<j<<" tables";
   for(auto v:data.at(abc))cout<<" ("<<v[0]<<","<<v[1]<<","<<v[2]<<")";
   cout<<endl;
  }
 }
}
