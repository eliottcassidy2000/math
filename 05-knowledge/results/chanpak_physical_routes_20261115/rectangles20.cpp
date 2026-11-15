// Complete one-apex additive-mask inverse at fixed width 20.
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>
using namespace std;
int main(){
 const int width=20, size=1<<width;
 vector<int> h(size,1);
 vector<vector<int>> pos(size+20);
 for(int x=0;x<size;x++){
  if(x)h[x]=h[x/2]+((x%2)?0:x/2+(x/2)%2);
  pos.at(h[x]).push_back(x);
 }
 int count=0,cube=0;set<tuple<int,int,int,int,int>> reduced;
 for(int w=0;w<size;w++){
  int base=h[w];
  for(int x:pos.at(base+2))for(int y:pos.at(base+4)){
   int A=x^w,B=y^w;
   if((A&B)||h[w^A^B]!=base+6)continue;
   ++count;
   for(int z:pos.at(base+8)){
    int C=z^w;
    if((C&(A|B))==0 && h[w^A^C]==base+10 && h[w^B^C]==base+12 && h[w^A^B^C]==base+14)
     ++cube;
   }
   int t=w,a=A,b=B,n=width;
   while(((t|a|b)>>(n-1)&1)==0)--n;
   while((t&1)&&!(a&1)&&!(b&1)){
    t>>=1; a>>=1;b>>=1;--n;
   }
   reduced.emplace(n,t,a,b,base);
  }
 }
 cout<<count<<" width-20 rectangles, "<<reduced.size()<<" reduced, "<<cube<<" +2/+4/+8 additive cubes\n";
 for(auto [n,w,A,B,H]:reduced)cout<<n<<" "<<w<<" "<<A<<" "<<B<<" "<<H<<"\n";
 if(count!=56 || reduced.size()!=6 || cube)return 1;
}
