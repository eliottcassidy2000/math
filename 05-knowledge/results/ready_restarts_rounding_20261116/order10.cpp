#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using A=std::array<unsigned,10>;
using S=std::array<int,1024>;
S subsets(const A& a,int n,int root=-1) {
 int d[1024][10]={}; S sum{};
 for(int i=0;i<n;++i) if(root<0 || root==i)d[1<<i][i]=1;
 for(int s=1;s<(1<<n);++s)for(int i=0;i<n;++i)if(int v=d[s][i]){
  sum[s]+=v;
  for(unsigned E=a[i]&~s;E;E&=E-1){unsigned b=E&-E;d[s|b][__builtin_ctz(b)]+=v;}
 }
 sum[0]=root<0; return sum;
}
A append(A a,int n,unsigned v){
 a[n]=v;
 for(int i=0;i<n;++i) if(!(v>>i&1))a[i]|=1u<<n;
 return a;
}
int main(int argc,char**argv) {
 assert(argc==2); std::ifstream f(argv[1]); std::string s; std::getline(f,s);
 unsigned nb=0, bases=0, scalar=0, good=0, lower=0;
 while(std::getline(f,s)){
  std::istringstream ss(s); unsigned long long code; unsigned h; char c;
  if(!(ss>>code>>c>>h) || c!=',')return 2;
  if(h>27) continue;
  ++nb; A A8{};
  for(int j=0;j<8;++j)A8[j]=(code>>(8*j))&255;
   const auto prefix=subsets(A8,8);
  if (prefix[255]!=h)return 3;
  for(unsigned v=0;v<256;++v){
   A A9=append(A8,8,v);
   if (subsets(A9,9)[511]!=27) continue;
   ++bases;
   for(unsigned w=0;w<512;++w){
    A A10=append(A9,9,w);
    auto H=subsets(A10,10);
    if(H[1023]!=127 && H[1023]!=125) continue;
    ++scalar;
    auto R=subsets(A10,10,9); int r=0,t=0;
    for(unsigned a=512;a<=1023;++a){
     int v=R[a]*H[1023^a];
     if(__builtin_popcount(a)%2) t+=v; else r+=v;
    }
    if(r==62&&t==89){
     if(H[1023]==127) ++good;
     else ++lower;
    }
   }
  }
 }
 if(nb!=103 || bases!=272 || scalar!=470 || lower!=11 || good!=0)return 4;
 std::cout<<"n8 classes H<=27: "<<nb<<"; n9 H27 base presentations: "<<bases
          <<"; total scalar 125 or127: "<<scalar<<"; exact profile H125: "<<lower<<"; exact profile H127: "<<good<<"\n";
}
