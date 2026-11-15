#include <bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
  const int X=59999;
  int w[X+1]={};
  for(int m=1;m<=X;m++){int k=m/2;w[m]=w[k]+((m&1)?k+(k&1):0);}
  array<int,X+1> masks={},smallmasks={},alt={},plus={};
  for(int a=1;a<argc;a++){
    ifstream s(argv[a]);int n,h,p;string code,sc;
    while(s>>n>>h>>p>>code>>sc){
     int q=p+2*h;
     if(q>=X) continue;
     for(int m=1,cap=(X-h)/(q-h);m<=cap;++m){
       int y=h+m*q-w[m]*h;
       if (y<=X){
         masks[y]|=1<<n;
         if(m<64)smallmasks[y]|=1<<n;
         if(m==1)plus[y]|=1<<n;
       }
     }
     for(int m:{3,6,1024}){
      int y=m*q-h;
      if(y<=X)alt[y]|=1<<n;
     }
    }
  }
  for(int y=69;y<=X;y+=2)cout<<y<<" "<<masks[y]<<" "<<smallmasks[y]<<" "<<alt[y]<<" "<<plus[y]<<"\n";
}
