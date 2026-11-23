#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
using V=vector<int>; using Trip=array<int,3>; using State=array<int,12>;
int idx(int f,int m,int a){return (f*3+m)*2+a;}
State advance(const State& old,int tag,int from){
 State s{};
 bool first=tag>>from&1,second=tag>>(1-from)&1;
 for(int f=0;f<2;f++)for(int m=0;m<3;m++)for(int a=0;a<2;a++){
  int z=old[idx(f,m,a)];if(!z)continue;
  s[idx(!first,m,a)]+=z;
  if(m||first)s[idx(f,second?1:2,a)]+=z;
  if(a||second)s[idx(f,m,1)]+=z;
 }
 return s;
}
int read(const State& s,int from){
 return s[idx(1,2,0)]+s[idx(1,2,1)]+(from?s[idx(1,0,0)]+s[idx(1,0,1)]:0);
}
int height(const V& tags){
 State p{},q{}; p[idx(1,0,0)]=q[idx(1,0,0)]=1;
 for(int t:tags){p=advance(p,t,0);q=advance(q,t,1);}
 return read(p,0)+read(q,1);
}
void allwords(int m, int step,int word,const State& p,const State& q,V& h){
 if(step==m){h[word]=read(p,0)+read(q,1);return;}
 for(int t=0;t<4;t++)allwords(m,step+1,word|(t<<(2*step)),advance(p,t,0),advance(q,t,1),h);
}
Trip coeff(int x,int y,int z){
 assert((z-3*y+2*x)%6==0);int C=(z-3*y+2*x)/6;
 assert((y-x-6*C)%2==0);int B=(y-x-6*C)/2;
 return {x-2*B-3*C,B,C};
}
int main(int argc,char**argv){
 const int lo=argc>1?stoi(argv[1]):7, hi=argc>2?stoi(argv[2]):9;
 if(lo<2||hi<lo||hi>11)return 2;
 for(int k=lo;k<=hi;k++){
  const int m=k+1,L=1<<(2*m); V h(L);
  State p{};p[idx(1,0,0)]=1;allwords(m,0,0,p,p,h);
  using Record=array<Trip,4>;
  map<Trip,int> lower;
  map<Trip,array<int,4>> rep;
  map<Trip,Record> data;
  int passing=0,exact=0,r_curved=0,r_based=0,r_rank3=0;
  const uint64_t universe=uint64_t(L)*m*(2*k)*(2*k-1)/2;
  for(int w=0;w<L;w++)for(int i=0;i<2*m;i++)for(int j=i+1;j<2*m;j++){
   const array<int,4> ww={w,w^(1<<i),w^(1<<j),w^(1<<i)^(1<<j)};
   const int top=7*h[w];
   if(h[ww[3]]!=top||h[ww[1]]<=h[w]||h[ww[2]]<=h[w]||h[ww[1]]>=top||h[ww[2]]>=top)continue;
   for(int gap=0;gap<m;gap++)if(i/2!=gap&&j/2!=gap){
    passing++;
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
     return d[0]+2*d[1]+3*d[2]>0 && d[1]+3*d[2]>=0 && d[2]>=0;
    };
    assert(positive(C,B)&&positive(C,D)&&positive(B,upper)&&positive(D,upper));
    lower[C]++;exact++;
    if(C[2]){r_based++;if((8*C[0]-B[0]-D[0])%6 || (8*C[1]-B[1]-D[1])%6 || (8*C[2]-B[2]-D[2])%6)r_rank3++;}
    if(8*C[2]-B[2]-D[2]){
      if(!r_curved) {
        cout<<"first R curvature "<<k<<": ";
        for(int t=0;t<m;t++)cout<<(w>>(2*t)&3);
        cout<<" gap bits "<<gap<<" "<<i<<" "<<j;
        for(auto v:U)cout<<" ("<<v[0]<<","<<v[1]<<","<<v[2]<<")";
        cout<<endl;
      }
      r_curved++;
    }
    if(!rep.count(C)){rep[C]={w,gap,i,j};data[C]=U;}
   }
  }
  cout<<"k="<<k<<" universe="<<universe<<" width1="<<passing<<" exact="<<exact<<" R_lower="<<r_based<<" R_curvature="<<r_curved<<" R_lower_d_not_coeffwise_3="<<r_rank3<<endl;
  for(auto [abc,count]:lower){
   auto [w,g,i,j]=rep.at(abc);
   cout<<"  "<<abc[0]<<","<<abc[1]<<","<<abc[2]<<" count="<<count<<" rep ";
   for(int z=0;z<m;z++)cout<<(w>>(2*z)&3);
   cout<<" gap,contact-bits="<<g<<","<<i<<","<<j<<" tables";
   for(auto v:data.at(abc))cout<<" ("<<v[0]<<","<<v[1]<<","<<v[2]<<")";
   cout<<endl;
  }
 }
}
