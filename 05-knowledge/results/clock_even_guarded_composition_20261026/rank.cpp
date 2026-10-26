#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;
using I=uint64_t;
const I MOD=1000000007;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
I power(I x,I e){I r=1;for(;e;e>>=1,x=x*x%MOD)if(e&1)r=r*x%MOD;return r;}
struct Perm { vector<int> p, destinations; int reverse=0,path=0; };
int main(){
 for(I d=2;d*d<=MOD;++d)need(MOD%d!=0,"rank modulus must be prime");
 cout<<"{\"modulus\":"<<MOD<<",\"atlases\":[";
 for(int n=2;n<=7;++n){
  int ix[8][8]={},m=0;vector<pair<int,int>> edges;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){ix[i][j]=ix[j][i]=m++;edges.push_back({i,j});}
  int N=1<<m,full=N-1;vector<Perm> ps;vector<int> p(n);iota(p.begin(),p.end(),0);
  do{Perm a;a.p=p;for(auto [i,j]:edges){int k=ix[p[i]][p[j]];a.destinations.push_back(k);if(p[i]>p[j])a.reverse|=1<<k;}
   for(int i=1;i<n;++i)a.path|=1<<ix[p[i-1]][p[i]];ps.push_back(a);
  }while(next_permutation(p.begin(),p.end()));
  vector<int> cl(N,-1),reps,orb;
  for(int t=0;t<N;++t)if(cl[t]<0){
   int id=reps.size(),sz=0;reps.push_back(t);
   for(auto&a:ps){int u=a.reverse;for(int b=t;b;b&=b-1)u^=1<<a.destinations[__builtin_ctz((unsigned)b)];
    need(cl[u]<0||cl[u]==id,"orbit overlap");if(cl[u]<0){cl[u]=id;++sz;}}
   orb.push_back(sz);
  }
  int d=reps.size();vector<vector<I>> L(d,vector<I>(d));vector<int> conv(d);
  for(int i=0;i<d;++i){int t=reps[i];conv[i]=cl[t^full];
   for(auto&a:ps){bool ok=true;for(int k=1;k<n;++k){int u=a.p[k-1],v=a.p[k];
     if(((t>>ix[u][v])&1)!=(u>v)){ok=false;break;}}
    if(ok)++L[i][cl[t^(full^a.path)]];
   }
  }
  need(accumulate(orb.begin(),orb.end(),0)==N,"labelled universe");
  vector<int> cr;for(int i=0;i<d;++i)if(i<=conv[i])cr.push_back(i);
  for(int i=0;i<d;++i)for(int j=0;j<d;++j){
   need(I(orb[i])*L[i][j]==I(orb[j])*L[j][i],"detailed balance");
   need(L[conv[i]][conv[j]]==L[i][j],"converse equivariance");
  }
  int q=cr.size(),rank=0;vector<vector<I>> basis(q);vector<I> v(d,1);vector<int> pivots;
  for(int k=0;k<q;++k){
   vector<I> w;for(int i:cr)w.push_back(v[i]);
   for(int j=0;j<q;++j)if(w[j]){
    if(basis[j].empty()){I inv=power(w[j],MOD-2);for(auto&x:w)x=x*inv%MOD;basis[j]=w;++rank;pivots.push_back(k);break;}
    I c=w[j];for(int a=j;a<q;++a)w[a]=(w[a]+MOD-c*basis[j][a]%MOD)%MOD;
   }
   vector<I> next(d);for(int i=0;i<d;++i)for(int j=0;j<d;++j)next[i]=(next[i]+L[i][j]*v[j])%MOD;v=move(next);
  }
  need(rank==q,"constant is not cyclic on converse quotient");
  if(n>2)cout<<",";
  cout<<"{\"n\":"<<n<<",\"labelled_tournaments\":"<<N<<",\"types\":"<<d
      <<",\"converse_orbits\":"<<q<<",\"rank_mod_prime\":"<<rank
      <<",\"last_pivot_power\":"<<pivots.back()<<"}";
 }
 cout<<"]}\n";
}
