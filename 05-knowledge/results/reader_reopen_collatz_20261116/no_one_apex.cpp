#include <cstdint>
#include <iostream>
int main(){
 for(std::int64_t H:{7747LL,16974403LL}){
  int cnt=0;
  for(int n=2; (1LL<<(n-2))+1<=H; n++){
   if((1LL<<n)<H)continue;
   int hits=0;
   auto go=[&](auto&&self,int i,std::int64_t A,int p,std::int64_t h)->void{
      if(i==n-1){hits+=h+A+p==H;return;}
      self(self,i+1,2*A+1,1,h);
      self(self,i+1,2*A,0,h+A+p);
   };
   go(go,1,1,1,1); cnt+=hits;
   std::cout<<H<<" length "<<n<<" hits "<<hits<<"\n";
  }
  if(cnt)return 1;
 }
}
