// Independent literal rooted profiles for every row of the ancestry census.
#define EAR_LIBRARY
#include "ear_probe.cpp"
#include <sstream>
int main(int argc,char**argv){need(argc==3,"ancestry file and profile output");std::ifstream in(argv[1]);std::ofstream out(argv[2]);std::string line;int id=0;
 while(std::getline(in,line)){std::istringstream s(line);I h;int n;s>>h>>n;G T(n);for(auto&x:T)s>>x;auto e=ends(T),st=ends(converse(T));unsigned full=(1u<<n)-1;std::vector<I>hc(1u<<n);hc[0]=1;
  for(unsigned mask=1;mask<=full;++mask)for(int v=0;v<n;++v)hc[mask]+=e[mask*n+v];need(hc[full]==h,"profile source count");
  for(int v=0;v<n;++v){I F=0;for(unsigned mask=1;mask<=full;++mask)F+=e[mask*n+v]*hc[full^mask];out<<id<<' '<<v<<' '<<h<<' '<<st[full*n+v]<<' '<<hc[full^(1u<<v)]<<' '<<F<<'\n';}
  ++id;
 }std::cout<<"profiled classes "<<id<<'\n';
}
