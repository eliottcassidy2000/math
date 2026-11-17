#include <array>
#include <iostream>
#include <cstdint>
using namespace std;
using ll=long long;
struct S {ll d,xp,xm,yp,ym;};
uint64_t hits[17]={}; bool once[17]={};
S add(S a, int q, int x, int y,ll z){
 return {q*a.d+(x-y)*z+((1-y)*(a.xp-a.xm)-(1-x)*(a.yp-a.ym))/2,
 (q+1)*a.xp+x*z,(q-1)*a.xm+x*z,(q+1)*a.yp+y*z,(q-1)*a.ym+y*z};
}
void run(int n,int X,int Y,S a,S b){
 if(a.d==5 && b.d==-2){
   ++hits[n];if(!once[n]){once[n]=true;cout<<"first "<<n<<" "<<X<<" "<<Y<<endl;}
 }
 if(n==12)return;
 for(int x=0;x<2;x++)for(int y=0;y<2;y++)run(n+1,X|(x<<n),Y|(y<<n),add(a,1,x,y,1),add(b,2,x,y,1LL<<n));
}
int main(){run(0,0,0,{},{});
for(int n=0;n<=12;n++)cout<<n<<" "<<hits[n]<<"\n";}
