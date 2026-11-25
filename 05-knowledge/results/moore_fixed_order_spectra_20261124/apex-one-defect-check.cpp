// Exact scan: add the same single core reversal at both ends of two apex seams.
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
using Z=unsigned __int128;
using Mask=unsigned long long;
struct Host {
 int a,b,n,w;
 Mask u,v;
 Host(Mask u0,Mask v0,int w0):w(w0),u(u0),v(v0){
  a=u?64-__builtin_clzll(u):0;b=v?64-__builtin_clzll(v):0;n=a+b+1;
 }
 bool bit(int i)const{
   if(i<a)return u>>(a-1-i)&1;
   if(i>a)return !(v>>(i-a-1)&1);
   return w;
 }
 bool arc(int i,int j)const{
   if(i==n)return i!=j && bit(j);
   if(j==n)return !bit(i);
   return i<j;
 }
 // Count all paths using the designated directed block i->j, not rejecting
 // it when it is the arc we plan to add.
 Z block(int i,int j)const{
  Z tot=0;
  const int place[3]={1,3,9};
  for(int rev=0;rev<2;rev++){
   std::array<int,2> in={!rev?i:n,!rev?n:i};
   std::array<int,2> out={!rev?j:n,!rev?n:j};
   std::array<Z,27> z{};z[0]=1;
   for(int x=0;x<n;x++)if(x!=i&&x!=j){
    std::array<Z,27> next{};
    for(int s=0;s<27;s++)if(z[s]!=0)for(int t=0;t<3;t++){
      int d=s/place[t]%3;
      if(d==0&&t>0&&!arc(out[t-1],x))continue;
      const int k=t==2||arc(x,in[t])?1:2;
      next[s+(k-d)*place[t]]+=z[s];
    }
    z.swap(next);
   }
   for(int s=0;s<27;s++)
    if(s%3!=2&&s/3%3!=2&&s/9%3!=2&&(s/3%3||arc(out[0],in[1])))tot+=z[s];
  }
  return tot;
 }
 Z jump(int i,int j)const{return block(j,i)-block(i,j);}
};
Z g(Mask u){
 Z h=0;for(int i=1;i<63;i++)if((u>>i)&1)for(int j=0;j<i;j++)if(!(u>>j&1))
   h+=Z(1)<<(i-j==1?1:i-j-1);
 return h;
}
Z height(const Host& t){return 1+g(t.u)+g(t.v)+2*Z(t.u)*t.v+
 (t.w? t.v+t.v%2:t.u+t.u%2);}
std::string show(Z x){
 std::string s;
 do{ s+='0'+x%10; x/=10; }while(x);
 std::reverse(s.begin(),s.end());return s;
}
int main(){
 // Independent direct permutation control on small cases.
 for(Mask u=0;u<8;u++)for(Mask v=0;v<8;v++)for(int w=0;w<2;w++){
  Host T(u,v,w);
  for(int i=0;i<T.n;i++)for(int j=i+1;j<T.n;j++){
   int h=0,h1=0;std::vector<int> p;
   for(int i0=0;i0<=T.n;i0++)p.push_back(i0);
   do{
     bool ok=true,ok1=true;
     for(unsigned k=1;k<p.size();k++){
      int a=p[k-1],b=p[k];bool x=T.arc(a,b);
      ok&=x;ok1&=x^((a==i&&b==j)||(a==j&&b==i));
     }
     h+=ok;h1+=ok1;
   }while(std::next_permutation(p.begin(),p.end()));
   assert(height(T)==h && T.jump(i,j)==Z(h1-h));
  }
 }
 const Mask x=1ULL<<59;
 const Z E=(Z(1)<<60)*((Z(1)<<61)-1);
 for(Mask q:{x-1,x})for(int offset:{-1,0})for(int after:{1,2}){
  std::pair<Mask,Mask> uv{2*q+offset,2*q+after};
  Host L(uv.first,uv.second,0),U(uv.first,uv.second,1);
  Z h=height(L),h1=height(U); assert(h1==h+2);
  int hits=0,lower=0,upper=0,same=0;
  for(int i=0;i<L.n;i++)for(int j=i+1;j<L.n;j++){
    Z d=L.jump(i,j),d1=U.jump(i,j);
    if(h+d==E-1)lower++;
    if(h1+d1==E+1)upper++;
    if(d==d1)same++;
    if(h+d==E-1&&h1+d1==E+1)hits++;
  }
  assert(hits==0 && lower==0 && upper==0);
  std::cout<<"flanks "<<uv.first<<","<<uv.second<<" order "<<L.n+1<<": midpoint error "<<(h+1<E ? "-"+show(E-h-1):show(h+1-E))
           <<", lower mod 2^128 hits "<<lower<<", upper mod 2^128 hits "<<upper
           <<", both "<<hits<<", unit-jump preserving mod 2^128 "<<same<<std::endl;
 }
}
