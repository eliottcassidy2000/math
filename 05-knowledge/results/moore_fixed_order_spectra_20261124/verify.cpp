#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using namespace std;
int main(int argc,char**argv) {
 int n=stoi(argv[1]),full=(1<<n)-1,F[2048][11],rows=0;
 bool fill=argc>2;string s;if(!fill)getline(cin,s);
 for(int i=0;i<n;i++)F[0][i]=0;
 while(getline(cin,s)){
  auto i=s.find(','),j=s.find(',',i+1);
  int H=stoi(s.substr(0,i));
  uint64_t mask=stoull(fill?s.substr(i+1,j-i-1):s.substr(j+1));
  int out[11]={};int k=0;
  for(int a=0;a<n;a++)for(int b=a+1;b<n;b++,k++)
   if(mask>>k&1)out[b]|=1<<a;else out[a]|=1<<b;
  for(int S=1;S<=full;S++)for(int B=S;B;B&=B-1){
    int v=__builtin_ctz(B),L=S^(1<<v),z=L==0;
    for(int a=L&~out[v];a;a&=a-1)z+=F[L][__builtin_ctz(a)];
    F[S][v]=z;
  }
  int x=0;for(int v=0;v<n;v++)x+=F[full][v];
  assert(x==H);rows++;
 }
 cout<<"order "<<n<<" checked "<<rows<<" witnesses by independent endpoint DP\n";
}
