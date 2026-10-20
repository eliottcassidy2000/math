#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;

// Exhaustive sum of path indicators over external-frame bit cubes.
// Bit order: A(i,j), i<j, followed by symmetric b(i,j).
int main(int argc, char **argv) {
  int k = argc > 1 ? stoi(argv[1]) : 5;
  assert(1 <= k && k <= 5);
  int m=k*(k-1)/2, nbits=2*m, size=1<<nbits;
  vector<vector<int>> ix(k,vector<int>(k,-1));
  int at=0;
  for(int i=0;i<k;i++)for(int j=i+1;j<k;j++)ix[i][j]=ix[j][i]=at++;
  vector<uint16_t> h(size,0);
  vector<int> p(k);iota(p.begin(),p.end(),0);
  long long paths=0,updates=0;
  auto start=chrono::steady_clock::now();
  do {
    for(int signs=0;signs<(1<<(k-1));signs++){
      unsigned fixed=0,ones=0;
      for(int z=0;z+1<k;z++){
        int i=p[z],j=p[z+1],si=(signs>>z)&1,sj=(signs>>(z+1))&1;
        int bit=ix[i][j],value;
        if(si!=sj){bit+=m;value=1-si;}
        else value=(i<j) != bool(si);
        assert(!(fixed&(1u<<bit)));
        fixed|=1u<<bit;if(value)ones|=1u<<bit;
      }
      unsigned free=((1u<<nbits)-1)^fixed,s=free;
      do {h[ones|s]++;updates++;if(!s)break;s=(s-1)&free;}while(true);
      paths++;
    }
  } while(next_permutation(p.begin(),p.end()));
  map<int,long long> hist;long long sum=0,squares=0;int mx=0;
  vector<long long> asums(1<<m,0);
  for(int f=0;f<size;f++){
    assert(h[f]%2);hist[h[f]]++;sum+=h[f];squares+=(long long)h[f]*h[f];
    mx=max(mx,(int)h[f]);asums[f&((1<<m)-1)]+=h[f];
  }
  cout << "{\"k\":"<<k<<",\"external_frames\":"<<size<<",\"potential_prefixes\":"<<paths
       <<",\"subcube_updates\":"<<updates<<",\"sum_h\":"<<sum<<",\"sum_h_squared\":"<<squares
       <<",\"max_h\":"<<mx<<",\"max_masks\":[";
  bool first=true;
  for(int f=0;f<size;f++)if(h[f]==mx){if(!first)cout<<",";cout<<f;first=false;}
  cout<<"],\"histogram\":{";first=true;
  for(auto [x,y]:hist){if(!first)cout<<",";cout<<"\""<<x<<"\":"<<y;first=false;}
  cout<<"},\"cross_sums_by_A\":[";first=true;
  for(auto x:asums){if(!first)cout<<",";cout<<x;first=false;}
  cout<<"],\"elapsed_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
  if(argc>2){ofstream out(argv[2],ios::binary);out.write((char*)h.data(),h.size()*sizeof(h[0]));}
}
