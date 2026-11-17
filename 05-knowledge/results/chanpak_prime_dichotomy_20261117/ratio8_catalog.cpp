#include <array>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>
using namespace std;
int main() {
  map<tuple<int,int,int>,int> out;
  for(int n=1;n<=7;++n){
    vector<pair<int,int>> free;
    for(int i=0;i<n;++i) for(int j=i+2;j<n;++j) free.emplace_back(i,j);
    int N=1<<n, F=N-1;
    for(int t=0;t<(1<<free.size());++t){
      array<int,7> a{};
      for(int i=0;i+1<n;++i) a[i]|=1<<(i+1);
      for(unsigned j=0;j<free.size();++j){
        auto [x,y]=free[j]; if(t>>j&1)swap(x,y); a[x]|=1<<y;
      }
      int D[128][7]{},H[128]{}; H[0]=1;
      for(int m=1;m<N;++m)for(int i=0;i<n;++i)if(m>>i&1){
        int b=m^(1<<i);
        if(!b)D[m][i]=1;
        else for(int j=0;j<n;++j)if((b>>j&1) && (a[j]>>i&1))D[m][i]+=D[b][j];
        H[m]+=D[m][i];
      }
      int g=H[F],p=0;
      for(int m=0;m<N;++m)p+=H[m]*H[F^m];
      if(p!=8*g)continue;
      int q=0;
      for(int m=0;m<N;++m){
        int b=F^m;
        for(int s=b;;s=(s-1)&b){q+=H[m]*H[s]*H[b^s]; if(s==0)break;}
      }
      ++out[{g,p,q}];
    }
  }
  for(auto [r,count]:out){auto [h,p,q]=r;cout<<h<<" "<<p<<" "<<q<<" "<<count<<"\n";}
}
