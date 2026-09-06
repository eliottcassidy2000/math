// Complete labelled cycle-space and dual-space enumeration, n<=8.
// No graph-isomorphism library, canonical-label heuristic, or sampling.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using U = uint32_t;
using V = std::vector<int>;
using M = std::vector<V>;
void need(bool p, const char* s) { if (!p) throw std::runtime_error(s); }
struct Linear {
  std::array<std::array<U,128>,3> tab{};
  explicit Linear(const std::vector<U>& cols) {
    for(int b=0;b<3;++b) for(int x=1;x<128;++x) {
      int j=__builtin_ctz(unsigned(x)), k=7*b+j;
      tab[b][x]=tab[b][x^(1<<j)] ^ (k<int(cols.size())?cols[k]:0);
    }
  }
  U operator()(U x) const { return tab[0][x&127]^tab[1][(x>>7)&127]^tab[2][x>>14]; }
};
struct Orbits { V label; std::vector<std::vector<U>> members; };
Orbits quotient(U count,const std::vector<Linear>& maps) {
  Orbits o; o.label.assign(count,-1);
  for(U seed=0;seed<count;++seed) if(o.label[seed]<0) {
    int id=int(o.members.size()); o.members.push_back({seed}); o.label[seed]=id;
    auto& q=o.members.back();
    for(size_t k=0;k<q.size();++k) for(const auto& f:maps) {
      U y=f(q[k]); need(y<count,"linear range");
      if(o.label[y]<0) { o.label[y]=id; q.push_back(y); }
      else need(o.label[y]==id,"complete generator orbit");
    }
  }
  return o;
}
template<class T> void array_json(const std::vector<T>& v) {
  std::cout<<'['; for(size_t i=0;i<v.size();++i) { if(i) std::cout<<','; std::cout<<v[i]; } std::cout<<']';
}
void matrix_json(const M& a) {
  std::cout<<'['; for(size_t i=0;i<a.size();++i) { if(i) std::cout<<','; array_json(a[i]); } std::cout<<']';
}
int main(int argc,char** argv) {
  need(argc==2,"one argument n required"); int n=std::stoi(argv[1]); need(3<=n&&n<=8,"n in3..8");
  int idx[8][8], gi[8][8];
  for(auto& row:idx) std::fill(row,row+8,-1);
  for(auto& row:gi) std::fill(row,row+8,-1);
  std::vector<std::pair<int,int>> edges,gauge;
  for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) { idx[i][j]=idx[j][i]=edges.size(); edges.push_back({i,j}); }
  for(int i=1;i<n;++i) for(int j=i+1;j<n;++j) { gi[i][j]=gi[j][i]=gauge.size(); gauge.push_back({i,j}); }
  int m=gauge.size(); U count=U(1)<<m;
  std::vector<U> basis;
  for(auto [i,j]:gauge) basis.push_back((U(1)<<idx[i][j])^(U(1)<<idx[0][i])^(U(1)<<idx[0][j]));
  std::vector<U> full(count,0);
  for(U x=1;x<count;++x) { int k=__builtin_ctz(x); full[x]=full[x^(U(1)<<k)]^basis[k]; }
  auto compress=[&](U x) { U a=0; for(int j=0;j<m;++j) if((x>>idx[gauge[j].first][gauge[j].second])&1) a|=U(1)<<j; return a; };
  std::vector<Linear> primal_maps,dual_maps;
  for(int v=0;v<n-1;++v) {
    V p(n); std::iota(p.begin(),p.end(),0); std::swap(p[v],p[v+1]);
    auto permute=[&](U x) { U a=0; for(int k=0;k<int(edges.size());++k) if((x>>k)&1) a^=U(1)<<idx[p[edges[k].first]][p[edges[k].second]]; return a; };
    std::vector<U> pc,dc;
    for(U x:basis) pc.push_back(compress(permute(x)));
    for(auto [a,b]:gauge) {
      U x=U(1)<<idx[p[a]][p[b]], h=0;
      for(int k=0;k<m;++k) {
        auto [i,j]=gauge[k];
        if(((x>>idx[i][j])^(x>>idx[0][i])^(x>>idx[0][j]))&1) h|=U(1)<<k;
      }
      dc.push_back(h);
    }
    primal_maps.emplace_back(pc); dual_maps.emplace_back(dc);
  }
  Orbits primal=quotient(count,primal_maps),dual=quotient(count,dual_maps);
  int q=primal.members.size(); need(q==int(dual.members.size()),"primal/dual class count");
  V order(q),relabel(q),sizes(q); std::vector<U> mins(q),reps(q),compressed(q);
  std::iota(order.begin(),order.end(),0);
  for(int i=0;i<q;++i) { mins[i]=UINT32_MAX; for(U x:primal.members[i]) mins[i]=std::min(mins[i],full[x]); }
  std::sort(order.begin(),order.end(),[&](int a,int b){return mins[a]<mins[b];});
  for(int i=0;i<q;++i) { relabel[order[i]]=i; reps[i]=mins[order[i]]; compressed[i]=compress(reps[i]); sizes[i]=primal.members[order[i]].size(); }
  for(int& i:primal.label) i=relabel[i];
  // An independent full-S_n path checks every canonical representative and
  // orbit size, including n8 (not just the generator implementation).
  V perm(n),stabilizers(q,0); std::iota(perm.begin(),perm.end(),0);
  int permutations=0;
  do {
    ++permutations;
    std::array<std::array<U,128>,4> tab{};
    for(int b=0;b<4;++b) for(int x=1;x<128;++x) {
      int j=__builtin_ctz(unsigned(x)),k=7*b+j;
      U bit=0;
      if(k<int(edges.size())) bit=U(1)<<idx[perm[edges[k].first]][perm[edges[k].second]];
      tab[b][x]=tab[b][x^(1<<j)]^bit;
    }
    for(int i=0;i<q;++i) {
      U r=reps[i],y=tab[0][r&127]^tab[1][(r>>7)&127]^tab[2][(r>>14)&127]^tab[3][r>>21];
      need(y>=r,"literal full-permutation canonical representative");
      stabilizers[i]+=(y==r);
    }
  } while(std::next_permutation(perm.begin(),perm.end()));
  for(int i=0;i<q;++i) need(int64_t(sizes[i])*stabilizers[i]==permutations,"literal full-permutation orbit size");
  std::vector<U> triangles;
  for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) for(int k=j+1;k<n;++k) triangles.push_back(compress((U(1)<<idx[i][j])^(U(1)<<idx[i][k])^(U(1)<<idx[j][k])));
  M weighted(q,V(q)),psi(q,V(q));
  for(int i=0;i<q;++i) for(U t:triangles) ++weighted[i][primal.label[compressed[i]^t]];
  for(U x=0;x<count;++x) {
    V row(q); for(U t:triangles) ++row[primal.label[x^t]];
    need(row==weighted[primal.label[x]],"all labelled rows match representatives");
  }
  for(int i=0;i<q;++i) for(int j=0;j<q;++j) need(int64_t(sizes[i])*weighted[i][j]==int64_t(sizes[j])*weighted[j][i],"multiplicity reversibility");
  V complement(q),dual_sizes(q),dual_reps(q),eigenvalues(q);
  for(int j=0;j<q;++j) {
    U seed=dual.members[j][0]; dual_reps[j]=seed; dual_sizes[j]=dual.members[j].size();
    complement[j]=dual.label[seed^(count-1)];
    for(U t:triangles) eigenvalues[j]+=1-2*(__builtin_popcount(seed&t)&1);
    for(int i=0;i<q;++i) for(U h:dual.members[j]) psi[i][j]+=1-2*(__builtin_popcount(h&compressed[i])&1);
  }
  std::cout<<"{\"n\":"<<n<<",\"states\":"<<count<<",\"literal_permutations\":"<<permutations<<",\"reps\":"; array_json(reps);
  std::cout<<",\"sizes\":"; array_json(sizes);
  std::cout<<",\"dual_reps\":"; array_json(dual_reps);
  std::cout<<",\"dual_sizes\":"; array_json(dual_sizes);
  std::cout<<",\"complement\":"; array_json(complement);
  std::cout<<",\"eigenvalues\":"; array_json(eigenvalues);
  std::cout<<",\"M\":"; matrix_json(weighted);
  std::cout<<",\"psi\":"; matrix_json(psi); std::cout<<"}\n";
}
