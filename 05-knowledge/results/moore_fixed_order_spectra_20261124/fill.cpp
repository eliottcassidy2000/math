#include <array>
#include <cassert>
#include <cstdio>
#include <random>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;
int main(int argc,char**argv){
 int n=11,full=(1<<n)-1,remain=0;
 bool want[15746]={},hit[15746]={};
 ifstream data(argv[1]);string s;getline(data,s);
 while(getline(data,s))hit[stoi(s)]=true;
 for(int h=1;h<15746;h+=2)if(!hit[h]&&h!=7&&h!=21)want[h]=true,remain++;
 mt19937_64 rng(20261124); int F[2048][11];for(int i=0;i<n;i++)F[0][i]=0;
 for(int t=0;t<1000000&&remain;t++){
  int in[11]={}; int p=(t%3==0)?64:(t%3==1?70:76);
  uint64_t key=0;int k=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){
    if(rng()%100<p)in[j]|=1<<i;else in[i]|=1<<j,key|=1ull<<k;
  }
  for(int S=1;S<=full;S++){
   int b=S;while(b){int i=__builtin_ctz(b);b&=b-1;
    int V=S^(1<<i),sum=V==0,W=V&in[i];
    while(W){int j=__builtin_ctz(W);W&=W-1;sum+=F[V][j];}
    F[S][i]=sum;
   }
  }
  int h=0;for(int i=0;i<n;i++)h+=F[full][i];
  if(h<15746&&want[h])want[h]=false,remain--,cout<<h<<","<<key<<","<<t<<"\n";
 }
 cerr<<"remaining "<<remain<<"\n";
 return !!remain;
}
