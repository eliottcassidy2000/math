#include <array>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>
using namespace std;
int main() {
  map<tuple<int,int,int,int>,int> out;
  bool failed=false;
  for(int n=3;n<=7;++n){
    vector<pair<int,int>> free;
    for(int i=0;i<n;++i) for(int j=i+2;j<n;++j) free.emplace_back(i,j);
    int N=1<<n, F=N-1;
    for(int t=0;t<(1<<free.size());++t){
      array<int,7> a{};
      for(int i=0;i+1<n;++i) a[i]|=1<<(i+1);
      for(unsigned j=0;j<free.size();++j){
        auto [x,y]=free[j]; if(t>>j&1)swap(x,y); a[x]|=1<<y;
      }
      // The prescribed forward Hamilton path already makes every vertex
      // reachable from vertex zero. Test that all can also reach zero.
      int back=1;
      for(int k=0;k<n;++k) for(int i=0;i<n;++i) if(a[i]&back)back|=1<<i;
      if(back!=F)continue;
      int S[128][7]{},E[128][7]{},H[128]{}; H[0]=1;
      for(int m=1;m<N;++m)for(int i=0;i<n;++i)if(m>>i&1){
        int b=m^(1<<i);
        if(!b)S[m][i]=E[m][i]=1;
        else for(int j=0;j<n;++j)if(b>>j&1){
          if(a[i]>>j&1)S[m][i]+=S[b][j];
          else E[m][i]+=E[b][j];
        }
        H[m]+=S[m][i];
      }
      int g=H[F], r=257-2*g;
      if(g<13 || g>32 || r%3)continue;
      int target=r/3;
      int marks[7]{};
      for(int m=1;m<N;++m)for(int i=0;i<n;++i)marks[i]+=S[m][i]*H[F^m];
      for(int q=1;q<N;++q){
        if(__builtin_popcount((unsigned)q)<2)continue;
        int W=0,e=0,U=0;
        for(int i=0;i<n;++i)if(q>>i&1){W+=marks[i];e+=E[F][i];}
        if(W!=target)continue;
        for(int m=1;m<F;++m){
          int cs=0,ce=0;
          for(int i=0;i<n;++i)if(q>>i&1){cs+=S[m][i];ce+=E[F^m][i];}
          U+=cs*ce;
        }
        ++out[{n,g,W,U+e-g}]; if(U+e==g)failed=true;
      }
    }
  }
  for(auto [a,count]:out){auto [n,g,w,d]=a;cout<<n<<" "<<g<<" "<<w<<" "<<d<<" "<<count<<"\n";}
  return failed?1:0;
}
