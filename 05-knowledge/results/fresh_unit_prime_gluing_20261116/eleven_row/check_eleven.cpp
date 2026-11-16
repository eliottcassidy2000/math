// Build: c++ -std=c++17 -O3 check_eleven.cpp -o /tmp/check_eleven
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;

constexpr int N=8, INF=1000000;
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
  struct Stat {int minimum=INF, hits=0, bad=0;};
  map<int,Stat> tab;
  for(int n=3;n<=7;++n){
    vector<pair<int,int>> es;
    for(int j=1;j<n;++j)for(int i=0;i<j;++i)es.emplace_back(i,j);
    int F=(1<<n)-1;
    for(int code=0;code<(1<<es.size());++code){
      Rows a{};
      for(int k=0;k<(int)es.size();++k){
        auto [i,j]=es[k]; if(!(code>>k&1))swap(i,j);
        a[i]|=1<<j;
      }
      int tau=n*(n-1)*(n-2)/6, prev=-1;
      for(int i=0;i<n;++i){
        int deg=__builtin_popcount(unsigned(a[i]));
        if(deg<prev){tau=INF; break;}
        prev=deg; tau-=deg*(deg-1)/2;
      }
      if(tau>5||!strong(a,n))continue;
      int q=Counts(a,n).h[F];
      if(q<11||!prime(q))continue;
      auto&t=tab[q];
      // Mixed contacts only. Other orientation bits are unchanged.
      for(int A=1;A<F;++A){
        Rows b=a; b[n]=A;
        for(int i=0;i<n;++i)if(!(A>>i&1))b[i]|=1<<n;
        Counts z(b,n+1);
        int d=(z.h[2*F+1]-q)/2;
        require(2*d==z.h[2*F+1]-q && d>0);
        t.minimum=min(t.minimum,d);
        if(d==5){
          ++t.hits;
          if(find(z.h,z.h+2*F+2,11)!=z.h+2*F+2)++t.bad;
        }
      }
    }
  }
  vector<pair<int,int>> expected{{11,4},{13,5},{17,5},{19,6},{23,9},{29,8},{31,10}};
  require(tab.size()==expected.size());
  cout<<"q min_half_gap pairs_at_5 pairs_at_5_with_induced_11\n";
  for(auto [q,m]:expected){
    require(tab.count(q)&&tab[q].minimum==m);
    auto [mn,hits,bad]=tab[q];
    if(q==13||q==17)require(hits==12&&bad==0);
    cout<<q<<" "<<mn<<" "<<hits<<" "<<bad<<"\n";
  }
}
