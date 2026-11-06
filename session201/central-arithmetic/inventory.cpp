// Exact finite arithmetic inventory, independent of the Python progression inverse.
// Build after: python verify.py --header central_forms.h
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "central_forms.h"
int main(int argc,char**argv){
  int N=argc>1?std::atoi(argv[1]):10000001;
  if(N<3 || N>100000001){std::cerr<<"cap outside documented range\n";return 2;}
  std::vector<int> g(N+1);g[1]=1;
  for(int k=2;k<=N;k++){
    if(k%2==0)g[k]=g[k/2];
    else{int h=(k+1)/2;g[k]=g[h]+h-h%2;}
  }
  std::vector<uint8_t> best(N+1,255);uint64_t tests=0;
  for(int f=0;f<FORMCOUNT;f++){
    for(int u=2;1LL*u*u*lower4[f]<=4LL*N;u++){
      long long source[]={g[u],u%2,u,1},row[4]={};
      for(int i=0;i<4;i++)for(int j=0;j<4;j++)row[j]+=source[i]*forms[f][i][j];
      for(int v=u;1LL*u*v*lower4[f]<=4LL*N;v++){
        long long n=row[0]*g[v]+row[1]*(v%2)+row[2]*v+row[3];tests++;
        if(n<3 || n%2!=1){std::cerr<<"invalid physical output\n";return 3;}
        if(n<=N)best[n]=std::min<int>(best[n],lengths[f]);
      }
    }
  }
  std::cout<<"{\n  \"cap\": "<<N<<",\n  \"matrices\": "<<FORMCOUNT
           <<",\n  \"pair_tests\": "<<tests<<",\n  \"by_maximum_central_length\": [\n";
  for(int length=0;length<=8;length++){
    int count=0,last=0;std::vector<int> gaps;
    for(int n=3;n<=N;n+=2)if(best[n]>length){count++;last=n;if(length>=6)gaps.push_back(n);}
    std::cout<<"    {\"length\": "<<length<<", \"gap_count\": "<<count<<", \"last_gap\": "<<last;
    if(length>=6){std::cout<<", \"gaps\": [";bool first=true;for(int n:gaps){if(!first)std::cout<<", ";first=false;std::cout<<n;}std::cout<<"]";}
    std::cout<<"}"<<(length<8?",":"")<<"\n";
  }
  std::cout<<"  ],\n  \"universal_occupancy\": \"OPEN\"\n}\n";
}
