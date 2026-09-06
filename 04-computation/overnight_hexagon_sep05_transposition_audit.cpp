// Independent exact controls for transposition rigidity, no imported helpers.
// Universe: every labelled unoriented Hamilton cycle on K_n, 5 <= n <= 11,
// every two-star split 0 <= r <= n-2; all switching classes for 6 <= n <= 8.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

void require(bool good, const char* message) {
  if (!good) throw std::runtime_error(message);
}
int64_t fact(int n) { int64_t out=1; for(int k=2;k<=n;k++)out*=k; return out; }
int64_t two_star(int n,int r) {
  int64_t s=n-2-r,q=r*s;
  return 2*q*((n-2)*(n-3)-2*q)*fact(n-5);
}

void hamilton_controls(int n) {
  int index[11][11]={}, dim=0;
  for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)index[u][v]=index[v][u]=dim++;
  std::vector<uint64_t> star(n-1);
  for(int r=0;r<=n-2;r++)for(int x=2;x<r+2;x++)
    star[r]|=(uint64_t(1)<<index[0][x])|(uint64_t(1)<<index[1][x]);
  std::vector<int64_t> count(n-1),adj(n-1),near(n-1),far(n-1);
  std::vector<int> p(n-1);std::iota(p.begin(),p.end(),1);
  int64_t total=0,edge_count=0,anti_count=0,negative_edge_count=0,hostile=0;
  uint64_t c5=0;
  if(n==6)for(int v=1;v<=5;v++)c5|=uint64_t(1)<<index[v][v%5+1];
  do {
    if(p.front()>p.back())continue;
    uint64_t mask=uint64_t(1)<<index[0][p.front()];
    mask|=uint64_t(1)<<index[0][p.back()];
    for(int j=0;j<n-2;j++)mask|=uint64_t(1)<<index[p[j]][p[j+1]];
    int position=int(std::find(p.begin(),p.end(),1)-p.begin())+1;
    ++total;
    bool edge=mask&(uint64_t(1)<<index[0][1]);
    edge_count+=edge;anti_count+=n%2;negative_edge_count+=(n%2)^edge;
    if(n==6)hostile+=__builtin_parityll(mask&c5);
    for(int r=0;r<=n-2;r++)if(__builtin_parityll(mask&star[r])) {
      ++count[r];
      if(position==1||position==n-1)++adj[r];
      else if(position==2||position==n-2)++near[r];
      else ++far[r];
    }
  } while(std::next_permutation(p.begin(),p.end()));
  require(total==fact(n-1)/2,"Hamilton universe count");
  require(edge_count==fact(n-2),"single-edge equality control");
  require(anti_count==(n%2)*total,"antibalanced parity control");
  require(negative_edge_count==((n%2)?total-edge_count:edge_count),"negative single-edge control");
  for(int r=0;r<=n-2;r++) {
    int64_t s=n-2-r,q=r*s;
    require(count[r]==two_star(n,r),"literal parity versus T(r)");
    require(adj[r]==2*q*fact(n-4),"adjacent-position count");
    require(near[r]==2*q*fact(n-4),"distance-two count");
    require(far[r]==2*q*((r-1)*(r-2)+(s-1)*(s-2))*fact(n-5),"distant-position count");
  }
  // Different universe: delete vertices 0,1 and count crossing gaps in K_(n-2).
  const int m=n-2;
  std::vector<int> rest(m-1);std::iota(rest.begin(),rest.end(),3);
  std::vector<int64_t> sx(n-1),sx2(n-1),insertions(n-1);
  int64_t base_count=0;
  do {
    if(rest.front()>rest.back())continue;
    std::vector<int> cycle={2};cycle.insert(cycle.end(),rest.begin(),rest.end());
    ++base_count;
    for(int r=0;r<=n-2;r++) {
      int x=0;
      for(int j=0;j<m;j++)x+=((cycle[j]<r+2)!=(cycle[(j+1)%m]<r+2));
      sx[r]+=x;sx2[r]+=x*x;insertions[r]+=2*x*(m+1-x);
    }
  } while(std::next_permutation(rest.begin(),rest.end()));
  require(base_count==fact(m-1)/2,"deleted-cycle universe count");
  for(int r=0;r<=n-2;r++) {
    int64_t q=int64_t(r)*(m-r);
    require(sx[r]*(m-1)==base_count*2*q,"crossing-gap first moment");
    require(sx2[r]*(m-1)*(m-2)==base_count*4*q*(q-1),"crossing-gap second moment");
    require(insertions[r]==two_star(n,r),"independent insertion versus T(r)");
  }
  std::cout<<"n="<<n<<" cycles="<<total<<" splits="<<n-1
           <<" parity_position_gap_moment_agree=1 signed_edge_controls=1\n";
  if(n==6) {
    require(hostile==20&&hostile<fact(n-2),"negative C5 plus positive apex hostile");
    std::cout<<"K6_negative_C5_positive_apex=20 single_edge=24 hostile_retained=1\n";
  }
}

void classification_control(int n) {
  int index[8][8]={},dim=0;
  for(int u=1;u<n;u++)for(int v=u+1;v<n;v++)index[u][v]=index[v][u]=dim++;
  const uint32_t end=uint32_t(1)<<dim,full=end-1;
  std::set<uint32_t> expected={0,full},edges;
  for(int j=0;j<dim;j++)edges.insert(uint32_t(1)<<j);
  for(int v=1;v<n;v++) {
    uint32_t star=0;
    for(int u=1;u<n;u++)if(u!=v)star|=uint32_t(1)<<index[u][v];
    edges.insert(star);
  }
  for(auto h:edges){expected.insert(h);expected.insert(full^h);}
  std::set<uint32_t> actual;
  for(uint32_t h=0;h<end;h++) {
    uint32_t row[8]={};
    for(int u=1;u<n;u++)for(int v=u+1;v<n;v++)if(h&(uint32_t(1)<<index[u][v])){
      row[u]|=uint32_t(1)<<v;row[v]|=uint32_t(1)<<u;
    }
    bool rigid=true;
    for(int u=0;u<n&&rigid;u++)for(int v=u+1;v<n;v++) {
      int r=__builtin_popcount((row[u]^row[v])&~((uint32_t(1)<<u)|(uint32_t(1)<<v)));
      if(r!=0&&r!=1&&r!=n-3&&r!=n-2){rigid=false;break;}
    }
    if(rigid)actual.insert(h);
  }
  require(actual==expected,"full row-disagreement rigidity classification");
  uint32_t matching=(uint32_t(1)<<index[1][2])|(uint32_t(1)<<index[3][4]);
  require(!actual.count(matching),"two matching edges hostile to rigidity");
  std::cout<<"classification_n="<<n<<" switching_classes="<<end
           <<" rigid_classes="<<actual.size()<<" equality_exact=1 matching_hostile=1\n";
}

int main() {
  for(int n=5;n<=11;n++)hamilton_controls(n);
  for(int n=6;n<=8;n++)classification_control(n);
  std::cout<<"RESULT=PASS\n";
}
