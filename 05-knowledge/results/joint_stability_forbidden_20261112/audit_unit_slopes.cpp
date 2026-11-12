#include <array>
#include <fstream>
#include <iostream>
using I=long long;
int main(int argc,char**argv){
 if(argc!=2)return 2;std::ifstream in(argv[1]);if(!in)return 3;
 constexpr int n=10,N=1024,full=N-1;
 std::array<unsigned,n> row{},pred{};std::array<std::array<I,n>,N> st{},en{};
 std::array<I,N> et{};std::array<I,n> U{};I expected,sources=0,marks=0,totalU=0;
 while(in>>expected){
  for(auto &r:row)if(!(in>>r))return 4;int count;if(!(in>>count))return 5;
  pred.fill(0);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(row[j]>>i&1)pred[i]|=1u<<j;
  st={};en={};et.fill(0);U.fill(0);
  for(int i=0;i<n;i++)st[1<<i][i]=en[1<<i][i]=1;
  for(unsigned m=1;m<N;m++){
   unsigned todo=m;
   while(todo){unsigned b=todo&-todo;todo-=b;int v=__builtin_ctz(b);unsigned rest=m^b;
    if(rest){unsigned a=rest&pred[v];while(a){unsigned bb=a&-a;a-=bb;en[m][v]+=en[rest][__builtin_ctz(bb)];}
     a=rest&row[v];while(a){unsigned bb=a&-a;a-=bb;st[m][v]+=st[rest][__builtin_ctz(bb)];}}
    et[m]+=en[m][v];
   }
  }
  I H=et[full];if(H!=expected||H<=1){std::cerr<<"H mismatch "<<sources<<'\n';return 6;}
  for(int j=0;j<n;j++)U[j]=st[full][j];
  for(unsigned m=1;m<full;m++)if(et[m]){unsigned a=full^m;while(a){unsigned bb=a&-a;a-=bb;int j=__builtin_ctz(bb);U[j]+=et[m]*st[full^m][j];}}
  for(int z=0;z<count;z++){
   unsigned O;I value;if(!(in>>O>>value))return 7;I actual=0;
   for(int i=0;i<n;i++)if(O>>i&1)actual+=U[i];
   I tax=0;
   for(unsigned mask=1;mask<full;mask++){
    I a=0,b=0;unsigned t=mask&O;while(t){unsigned bb=t&-t;t-=bb;a+=en[mask][__builtin_ctz(bb)];}
    if(!a)continue;t=(full^mask)&O;while(t){unsigned bb=t&-t;t-=bb;b+=st[full^mask][__builtin_ctz(bb)];}tax+=a*b;
   }
   if(actual!=value||H-tax!=1){std::cerr<<"interface mismatch "<<sources<<' '<<O<<'\n';return 8;}
   marks++;totalU+=actual;
  }
  sources++;
 }
 std::cout<<"{\"sources\":"<<sources<<",\"marked_interfaces\":"<<marks<<",\"sum_U\":"<<totalU<<",\"all_H_delta_U_equal\":true}\n";
 return 0;
}
