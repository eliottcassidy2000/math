#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

// Independent forward subset DP, compared against every subcube-summed frame.
int main(int argc,char**argv){
  assert(argc==3);int k=stoi(argv[1]);assert(k>=1&&k<=5);
  int m=k*(k-1)/2,size=1<<(2*m),n=2*k;
  vector<uint16_t> expected(size);
  ifstream in(argv[2],ios::binary);
  in.read((char*)expected.data(),size*sizeof(uint16_t));assert(in.gcount()==size*2);
  long long sum=0,squares=0;int mx=0,countmax=0;
  auto start=chrono::steady_clock::now();
  for(int code=0;code<size;code++){
    int adj[10]={};int z=0;
    auto edge=[&](int u,int v,int yes){if(yes)adj[u]|=1<<v;else adj[v]|=1<<u;};
    for(int i=0;i<k;i++)for(int j=i+1;j<k;j++,z++){
      int a=(code>>z)&1,b=(code>>(m+z))&1;
      edge(2*i,2*j,a);edge(2*j+1,2*i+1,a);
      edge(2*i,2*j+1,b);edge(2*j,2*i+1,b);
    }
    int dp[32][10]={};
    for(int i=0;i<n;i++)dp[1<<(i/2)][i]=1;
    for(int mask=1;mask<(1<<k);mask++){
      int avail=0;
      for(int j=0;j<k;j++)if(!(mask&(1<<j)))avail|=3<<(2*j);
      for(int i=0;i<n;i++)if(dp[mask][i]){
        int next=adj[i]&avail;
        while(next){
          int j=__builtin_ctz(next);next&=next-1;
          dp[mask|(1<<(j/2))][j]+=dp[mask][i];
        }
      }
    }
    int h=0;for(int i=0;i<n;i+=2)h+=dp[(1<<k)-1][i];
    assert(h==expected[code]);assert(h%2==1);
    sum+=h;squares+=(long long)h*h;
    if(h>mx){mx=h;countmax=1;}else if(h==mx)countmax++;
  }
  cout<<"{\"k\":"<<k<<",\"frames_checked\":"<<size
      <<",\"all_equal\":true,\"max_h\":"<<mx<<",\"maximizing_labelled_frames\":"<<countmax
      <<",\"sum_h\":"<<sum<<",\"sum_h_squared\":"<<squares
      <<",\"elapsed_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
