#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
using U=std::uint32_t;
U step(U x,int p,int sig,bool w){
 const int d=(p-1)/2;
 int prev=0,here=x%p; x/=p;
 U out=0,weight=1;
 for(int s=1;s<=d;++s){
  int next=s==d?sig*here:int(x%p);
  out+=weight*((s*s*here+(w?s*(next-prev):0)+p*p)%p);
  prev=here;here=next;x/=p;weight*=p;
 }
 return out;
}
int det(int p,int sig){
 int a=0,b=1;
 for(int i=1;i<=(p-1)/2;++i){
  int diag=i*i+(i==(p-1)/2?sig*i:0);
  int c=(diag*b+i*(i-1)*a)%p;
  a=b;b=c;
 }
 return b;
}
bool prime(int p){
 for(int d=2;d*d<=p;++d)if(p%d==0)return false;
 return p>1;
}
int main(){
 const std::vector<int> primes={3,5,7,11,13};
 const std::vector<std::string> shortwords={"","WWWWPW","WPPWWPP","WPPPWPPPPWPPPWWWPWWP","WPWPWPWWWPPPWWPWWPW"};
 for(unsigned i=0;i<primes.size();++i){
  int p=primes[i];U size=1,start=0;
  for(int s=1;s<=(p-1)/2;++s){start+=size*(s*s%p);size*=p;}
  std::vector<int> dist(size,-1);std::vector<U> queue{start};queue.reserve(size);dist[start]=0;
  int longest=0;
  for(unsigned j=0;j<queue.size();++j){
   U x=queue[j];
   for(bool w:{false,true}){
    U y=step(x,p,1,w);
    if(dist[y]!=-1)continue;
    dist[y]=dist[x]+1;longest=std::max(longest,dist[y]);queue.push_back(y);
   }
  }
  assert(queue.size()==size-1&&dist[0]==-1&&det(p,1)!=0);
  U x=1;
  for(char ch:shortwords[i])x=step(x,p,1,ch=='W');
  assert(x==size/p);
  std::cout<<p<<" / "<<queue.size()<<" / longest "<<longest
           <<" / maximum order "<<2+2*longest<<" / determinant "<<det(p,1)<<'\n';
 }
 for(U x=0;x<343;++x){
  int z=(6*(x%7)+5*((x/7)%7)+2*(x/49))%7;
  assert(step(step(x,7,-1,true),7,-1,true)==U(z*50+7*((3*z)%7)));
  assert(step(z*50+7*((3*z)%7),7,-1,true)==U(((4*z)%7)*50+7*((12*z)%7)));
 }
 for(int sig:{1,-1}){
  std::cout<<"singular up to 999, parity "<<sig<<":";
  for(int p=3;p<1000;p+=2)if(prime(p)&&!det(p,sig))std::cout<<' '<<p;
  std::cout<<'\n';
 }
}
