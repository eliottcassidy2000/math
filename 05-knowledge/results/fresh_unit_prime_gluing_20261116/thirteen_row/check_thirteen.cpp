// Build: c++ -std=c++17 -O3 check_thirteen.cpp -o /tmp/check_thirteen
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;

constexpr int N=9, INF=1000000;
using Rows=array<int,N>;
struct Counts{
  int h[1<<N]{};
  Counts(const Rows&a,int n){
    int d[1<<N][N]{};
    h[0]=1;
    for(int S=1;S<(1<<n);++S){
      for(int v=0;v<n;++v) if(S>>v&1){
        int R=S^(1<<v);
        if(!R) d[S][v]=1;
        else for(int u=0;u<n;++u)
          if((R>>u&1)&&(a[u]>>v&1)) d[S][v]+=d[R][u];
        h[S]+=d[S][v];
      }
    }
  }
};
bool strong(const Rows&a,int n){
  for(int dir=0;dir<2;++dir){
    int R=1, old=0;
    while(R!=old){
      old=R;
      for(int i=0;i<n;++i) if(old>>i&1){
        if(!dir) R|=a[i];
        else for(int j=0;j<n;++j) if(a[j]>>i&1) R|=1<<j;
      }
    }
    if(R!=(1<<n)-1) return false;
  }
  return true;
}
bool prime(int p){if(p<2)return false; for(int j=2;j*j<=p;++j)if(p%j==0)return false;return true;}
void require(bool b){ if(!b)throw runtime_error("certificate mismatch"); }
int main(){
  struct Stat {int minimum=INF, hits6=0, has13=0, hits5=0, pad=INF, hits4=0, pad4=INF; int n=0, code=0, A=0;};
  map<int,Stat> tab;
  for(int n=3;n<=8;++n){
    vector<pair<int,int>> es;
    for(int j=1;j<n;++j)for(int i=0;i<j-1;++i)es.emplace_back(i,j);
    int F=(1<<n)-1;
    for(int code=0;code<(1<<es.size());++code){
      Rows a{};
      // A Hamilton path 0 -> ... -> n-1 exists in every labeling considered.
      for(int i=0;i<n-1;++i)a[i]=1<<(i+1);
      for(int k=0;k<(int)es.size();++k){
        auto [i,j]=es[k]; if(!(code>>k&1))swap(i,j);
        a[i]|=1<<j;
      }
      int tau=n*(n-1)*(n-2)/6;
      for(int i=0;i<n;++i){
        int deg=__builtin_popcount(unsigned(a[i]));
        tau-=deg*(deg-1)/2;
      }
      if(tau>6||!strong(a,n))continue;
      Counts old(a,n);
      int q=old.h[F];
      if(q<11||!prime(q))continue;
      auto&t=tab[q];
      for(int A=1;A<F;++A){
        Rows b=a; b[n]=A;
        for(int i=0;i<n;++i)if(!(A>>i&1))b[i]|=1<<n;
        Counts z(b,n+1);
        int d=(z.h[2*F+1]-q)/2;
        require(2*d==z.h[2*F+1]-q && d>0);
        t.minimum=min(t.minimum,d);
        if(d==6){
          if(!t.hits6){t.n=n; t.code=code; t.A=A;}
          ++t.hits6;
          if(find(z.h,z.h+2*F+2,13)!=z.h+2*F+2)++t.has13;
        }
        if(d==4||d==5){
          int sums[2]{};
          for(int i=0;i<n;++i)sums[A>>i&1]+=old.h[F^(1<<i)];
          int bnd=min({sums[0],sums[1],q});
          if(d==4){++t.hits4; t.pad4=min(t.pad4,bnd);}
          if(d==5){++t.hits5; t.pad=min(t.pad,bnd);}
        }
      }
    }
  }
  vector<pair<int,int>> expected{{11,4},{13,5},{17,5},{19,6},{23,7},{29,8},{31,10},{37,16},{41,14},{53,14}};
  require(tab.size()==expected.size());
  cout<<"q min_half_gap entries_at_6 has_induced13 entries_at_5 min_pad_at_5 entries_at_4 min_pad_at_4 witness(n,code,A)\n";
  for(auto [q,m]:expected){
    require(tab.count(q)&&tab[q].minimum==m);
    auto [mn,h6,bad,h5,pad,h4,pad4,n,code,A]=tab[q];
    if(q==11)require(h6==66&&h4==22&&h5==0&&pad4==5&&bad==0);
    if(q==13)require(h6==52&&h4==0&&h5==26&&pad==5);
    if(q==17)require(h6==68&&h5==51&&bad==0&&pad==9);
    if(q==19)require(h6==38&&bad==0);
    cout<<q<<" "<<mn<<" "<<h6<<" "<<bad<<" "<<h5<<" "<<(pad==INF?-1:pad)<<" "<<h4<<" "<<(pad4==INF?-1:pad4)<<" "<<n<<","<<code<<","<<A<<"\n";
  }
  // Positive q=13 using literal host code 1105 in the standard upper-bit order.
  Rows w{};
  for(int j=1,k=0;j<6;++j)for(int i=0;i<j;++i,++k){
    if((1105>>k)&1)w[i]|=1<<j;
    else w[j]|=1<<i;
  }
  Counts y(w,6);
  require(y.h[63]==25&&y.h[63^(1<<4)]==13&&y.h[63^(1<<5)]==13);
}
