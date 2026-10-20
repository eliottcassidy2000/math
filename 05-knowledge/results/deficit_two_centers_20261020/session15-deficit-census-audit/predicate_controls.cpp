#define main census_main
#include "independent_scan.cpp"
#undef main
bool cuts_strong(const std::vector<U>&a,U full){
 for(U s=(full-1)&full;s;s=(s-1)&full){U other=full^s;bool edge=false;
  for(int i=0;i<(int)a.size();i++)if((s>>i)&1)edge|=(a[i]&other)!=0;
  if(!edge)return false;
 }return full!=0;
}
int main(){try{
 for(int n=1;n<=6;n++){
  uint64_t total=uint64_t(1)<<(n*(n-1)/2),strong=0;
  for(uint64_t mask=0;mask<total;mask++){
   std::vector<U>a(n);int bit=0;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,bit++)if((mask>>bit)&1)a[i]|=U(1)<<j;else a[j]|=U(1)<<i;
   U full=(U(1)<<n)-1;bool x=closure_strong(a,full),y=cuts_strong(a,full);
   need(x==y,"closure/cut disagreement");strong+=x;
   if(n==4){std::array<int,4>s={0,1,2,3};need(fourcycle(a,s)==y,"4-cycle/strongness disagreement");
    bool noextreme=true;for(U row:a)noextreme&=(std::popcount(row)==1||std::popcount(row)==2);
    need(noextreme==y,"degree proxy fails on four vertices");}
  }
  std::cout<<"{\"order\":"<<n<<",\"labelled\":"<<total<<",\"strong\":"<<strong<<",\"closure_matches_all_cuts\":true}"<<std::endl;
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
