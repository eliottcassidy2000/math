#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using Row = std::array<int,12>;

int modpow(int a,int e,int p){ long long r=1,b=(a%p+p)%p; while(e){if(e&1)r=r*b%p;b=b*b%p;e>>=1;}return (int)r; }

int rank_mod(const std::vector<Row>& source,const std::vector<int>& cols,int p){
  std::vector<std::vector<int>> a(source.size(),std::vector<int>(cols.size()));
  for(size_t i=0;i<source.size();++i)for(size_t j=0;j<cols.size();++j)a[i][j]=(source[i][cols[j]]%p+p)%p;
  int rr=0;
  for(int c=0;c<(int)cols.size();++c){
    int pivot=-1;for(int i=rr;i<(int)a.size();++i)if(a[i][c]){pivot=i;break;}
    if(pivot<0)continue;
    std::swap(a[rr],a[pivot]); int inv=modpow(a[rr][c],p-2,p);
    for(int j=0;j<(int)cols.size();++j)a[rr][j]=(int)((long long)a[rr][j]*inv%p);
    for(int i=0;i<(int)a.size();++i)if(i!=rr&&a[i][c]){int q=a[i][c];for(int j=0;j<(int)cols.size();++j)a[i][j]=(a[i][j]-(long long)q*a[rr][j]%p+p)%p;}
    if(++rr==(int)a.size())break;
  }
  return rr;
}

std::pair<bool,std::string> closure(const std::vector<Row>& rows,int p){
  int live=63; std::string trace;
  while(live){
    std::vector<int> cols;for(int v=0;v<6;++v)if(live&(1<<v)){cols.push_back(2*v);cols.push_back(2*v+1);}
    int r=rank_mod(rows,cols,p), fired=0;
    for(int v=0;v<6;++v)if(live&(1<<v)){
      std::vector<int> cut;for(int c:cols)if(c!=2*v+1)cut.push_back(c);
      if(rank_mod(rows,cut,p)==r-1)fired|=1<<v;
    }
    if(!fired)break;
    if(!trace.empty())trace.push_back('/');
    for(int v=0;v<6;++v)if(fired&(1<<v))trace.push_back(char('A'+v));
    live&=~fired;
  }
  return {live==0,trace};
}

int main(int argc,char**argv){
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  int p=argc>1?std::stoi(argv[1]):3;
  bool gauge=argc>2&&std::string(argv[2])=="gauge";
  const std::array<std::pair<int,int>,7> edges={{{0,3},{1,4},{2,5},{0,1},{1,2},{3,4},{4,5}}};
  uint64_t tests=0,successes=0;std::map<std::string,uint64_t> traces;std::string first;
  int total=1;for(int i=0;i<(gauge?7:8);++i)total*=p;
  std::vector<std::array<int,3>> positions;
  for(int i=0;i<6;++i)for(int j=i;j<6;++j)for(int k=j;k<6;++k)positions.push_back({i,j,k});
  for(int code=0;code<total;++code){
    int z=code;std::array<int,8> q{};
    if(gauge){q[0]=0;for(int i=1;i<8;++i){q[i]=z%p;z/=p;}}
    else for(int i=0;i<8;++i){q[i]=z%p;z/=p;}
    std::array<int,7> es={q[0],q[0],q[0],q[1],q[2],q[3],q[4]};
    for(auto pos:positions){
      std::vector<Row> rows;
      for(int e=0;e<7;++e){Row row{};auto [u,v]=edges[e];int rho=es[e];row[2*u]=1;row[2*u+1]=rho;row[2*v]=-1;row[2*v+1]=-rho;rows.push_back(row);}
      for(int i=0;i<3;++i){Row row{};int v=pos[i],rho=q[5+i];row[2*v]=1;row[2*v+1]=rho;rows.push_back(row);}
      ++tests;auto [ok,tr]=closure(rows,p);if(ok){++successes;++traces[tr];if(first.empty())first=std::to_string(code)+":"+std::to_string(pos[0])+std::to_string(pos[1])+std::to_string(pos[2])+":"+tr;}
    }
  }
  std::cout<<"p="<<p<<";gauge_alpha_zero="<<(gauge?1:0)<<";tests="<<tests<<";successes="<<successes<<"\n";
  std::cout<<"traces=";for(auto const& [tr,n]:traces)std::cout<<tr<<":"<<n<<",";std::cout<<"\nfirst="<<first<<"\n";
}
