// Hostile independent audit of the candidate all-order H >= u+|v| proof.
// Full adjacency determinants by modular Gaussian elimination; for n<=10,
// Hadamard reconstruction is exact because |det|<=100000<1000003/2.
// Hamilton counts by subset DP, insertion terms by literal directed-path DFS.
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
using I=long long;using Matrix=std::vector<std::vector<I>>;
constexpr I prime=1000003;
I power(I a,I e){I b=1;for(;e;e>>=1,a=a*a%prime)if(e&1)b=b*a%prime;return b;}
I determinant(Matrix a){int n=a.size();I d=1;for(int j=0;j<n;++j){int r=j;while(r<n&&!a[r][j])++r;if(r==n)return 0;if(r!=j){std::swap(a[r],a[j]);d=prime-d;}I p=a[j][j];d=d*p%prime;I inv=power(p,prime-2);for(int i=j+1;i<n;++i){I f=a[i][j]*inv%prime;for(int k=j;k<n;++k)a[i][k]=(a[i][k]-f*a[j][k]%prime+prime)%prime;}}return d>prime/2?d-prime:d;}
void require(bool c,const char*m){if(!c)throw std::runtime_error(m);}
struct Stats{I tournaments=0,roots=0,paths=0,subsets=0,empty_complements=0;std::uint64_t digest=1469598103934665603ULL;} stats;
void mix(I x){stats.digest^=static_cast<std::uint64_t>(x);stats.digest*=1099511628211ULL;}
void audit(const Matrix&a,bool all_roots){
  int n=a.size(),N=1<<n,full=N-1;std::vector<I> u(N),v(N),h(N);u[0]=v[0]=h[0]=1;
  std::vector<std::vector<I>> dp(N,std::vector<I>(n));for(int i=0;i<n;++i)dp[1<<i][i]=1;
  for(int s=1;s<N;++s){
    for(int j=0;j<n;++j){h[s]+=dp[s][j];if(dp[s][j])for(int k=0;k<n;++k)if(!(s>>k&1)&&a[j][k])dp[s|1<<k][k]+=dp[s][j];}
    std::vector<int> vs;for(int i=0;i<n;++i)if(s>>i&1)vs.push_back(i);
    int m=vs.size();Matrix b(m,std::vector<I>(m));for(int i=0;i<m;++i)for(int j=0;j<m;++j)b[i][j]=a[vs[i]][vs[j]];
    v[s]=determinant(b)*(m%2?-1:1);for(int i=0;i<m;++i)b[i][i]=1;u[s]=determinant(b);
    require(u[s]+v[s]>0&&u[s]-v[s]>0,"nonempty algebraic positivity");
    require(h[s]>=u[s]+std::abs(v[s]),"global P inequality");++stats.subsets;
  }
  for(int x=all_roots?0:n-1;x<n;++x){
    int parent=full^(1<<x);I eu=0,ou=0,sv=0,eh=0,plus=0,minus=0;
    std::function<void(int,int)> visit=[&](int last,int used){
      int length=__builtin_popcount(static_cast<unsigned>(used));int rest=parent^used;
      if(length>=2&&a[last][x]){
        ++stats.paths;if(!rest)++stats.empty_complements;
        sv+=v[rest];
        if(length%2==0){eu+=u[rest];eh+=h[rest];plus+=2*h[rest]-u[rest]+v[rest];minus+=2*h[rest]-u[rest]-v[rest];}
        else{ou+=u[rest];plus+=u[rest]+v[rest];minus+=u[rest]-v[rest];}
      }
      for(int j=0;j<n;++j)if(j!=x&&!(used>>j&1)&&a[last][j])visit(j,used|1<<j);
    };
    for(int j=0;j<n;++j)if(a[x][j])visit(j,1<<j);
    require(u[full]==u[parent]+eu-ou,"u insertion sign");
    require(v[full]==-sv,"v insertion sign");
    require(h[full]==h[parent]+2*eh,"H insertion OCF");
    require(h[full]-u[full]-v[full]==h[parent]-u[parent]+plus,"plus gap identity");
    require(h[full]-u[full]+v[full]==h[parent]-u[parent]+minus,"minus gap identity");
    ++stats.roots;
  }
  bool triangle=false;for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k)triangle|=a[i][j]&&a[j][k]&&a[k][i];
  if(!triangle)require(h[full]==1&&u[full]==1&&v[full]==0,"transitive equality");
  else if(n==3)require(h[full]==3&&u[full]==2&&v[full]==-1,"C3 equality");
  else require(h[full]>u[full]+std::abs(v[full]),"nontransitive strict P");
  require((h[full]==u[full]+v[full])==!triangle,"disc equality locus");
  mix(n);mix(h[full]);mix(u[full]);mix(v[full]);++stats.tournaments;
}
Matrix coded(int n,std::uint64_t code){Matrix a(n,std::vector<I>(n));int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k){a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j];}return a;}
int main(){
  for(int n=1;n<=6;++n)for(std::uint64_t c=0;c<(1ULL<<(n*(n-1)/2));++c)audit(coded(n,c),n<=5);
  std::uint64_t seed=0x1371950C0FFEEULL;
  for(int n=7;n<=10;++n)for(int j=0;j<128;++j){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;audit(coded(n,seed),true);}
  std::vector<int> rows={30,8,122,80,2,91,19};Matrix hostile(7,std::vector<I>(7));for(int i=0;i<7;++i)for(int j=0;j<7;++j)hostile[i][j]=rows[i]>>j&1;audit(hostile,true);
  std::cout<<"{\n  \"scope\": \"all labelled tournaments n1..6; all roots n<=5; 128 deterministic tournaments at each n7..10 with all roots; exact s<1 hostile\",\n"
           <<"  \"tournaments\": "<<stats.tournaments<<",\n  \"root_insertions\": "<<stats.roots<<",\n  \"literal_old_paths\": "<<stats.paths<<",\n  \"induced_determinant_checks\": "<<stats.subsets<<",\n  \"empty_complement_terms\": "<<stats.empty_complements<<",\n  \"semantic_fnv64\": \""<<stats.digest<<"\",\n  \"verdict\": \"PASS all three insertion identities, both signed gaps, algebraic positivity, and equality loci\"\n}\n";
  return 0;
}
