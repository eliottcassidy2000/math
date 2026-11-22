#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using I=__int128;using U=unsigned __int128;
using namespace std;
I P(int k){return I(1)<<k;}
int n,Width,Run;long long c,h,k0,k1;string mask;
uint64_t pairs=0; const I N=P(121)-P(60)-1; int hits=0;
int minZero=1, minTerminal=1;
void go(int b,int remaining,I y,I f){
 if (b==Run){
   if(remaining==0){
    ++pairs;
    if((y&1?k1:k0)+c*y+h*f==N){
     I g=0;for(I t=y;t;t>>=1)if(!(t&1))g+=(t>>1)+((t>>1)&1);
     assert(f==g);
     cout<<"exact hit "<<mask<<"\n";++hits;
    }
   }
   return;
 }
 int low=b==1?minZero:(b==Run-1&&b%2==0?minTerminal:1);
 int lo=b==Run-1?max(remaining,low):low;
 int hi=remaining-(Run-b-1);
 for(int g=lo;g<=hi;++g){
  if(b&1)go(b+1,remaining-g,y<<g,f+(P(g)-1)*y+1);
  else go(b+1,remaining-g,(y+1)*P(g)-1,f);
 }
}
int main(int argc,char**argv){
 const int order=argc>1?stoi(argv[1]):123;
 const int cap=argc>2?stoi(argv[2]):6;
 const int floor=argc>3?stoi(argv[3]):6;
 int examined=0;
 while(cin>>n>>c>>h>>k0>>k1>>mask){
   Width=order-7-n;
   if(Width<=0||c*P(Width-1)+min(k0,k1)>N)continue;
   if((c+h)*(P(Width)-1)+max(k0,k1)<N)continue;
   ++examined;
   I maxY=(N-min(k0,k1))/c, Y=P(Width-1),room=N-min(k0,k1)-c*Y;
   minZero=1; minTerminal=1;
   while(minZero<Width&&maxY<Y+(Y>>(minZero+1)))++minZero;
   int t=0;
   while(t+1<Width&&(Y>>(t+1))>room/h)++t;
   minTerminal=max(1,t);
   for(Run=floor;Run<=cap;++Run){
      if(t&&Run%2==0)continue;
      const int maxa=maxY<Y+(Y>>1)?1:Width-Run+1;
      for(int a=1;a<=maxa&&a<=Width-Run+1;++a)go(1,Width-a,P(a)-1,0);
   }
 }
 cerr<<"order="<<order<<" cap="<<cap<<" floor="<<floor<<" charts="<<examined<<" tested="<<pairs<<" hits="<<hits<<"\n";
}
