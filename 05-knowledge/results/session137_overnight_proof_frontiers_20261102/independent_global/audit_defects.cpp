// Independent direct all-cycle enumeration for the stronger integer defects.
#define main insertion_audit_main
#include "global_gap.cpp"
#undef main
I checked=0,pack2=0,total_cycles=0,total_proper=0;
void defects(const Matrix&a){
  int n=a.size(),N=1<<n;Matrix b=a;for(int i=0;i<n;++i)b[i][i]=1;
  I u=determinant(b),v=determinant(a)*(n%2?-1:1);
  std::vector<std::vector<I>> dp(N,std::vector<I>(n));for(int i=0;i<n;++i)dp[1<<i][i]=1;
  for(int s=1;s<N;++s)for(int j=0;j<n;++j)if(dp[s][j])for(int k=0;k<n;++k)if(!(s>>k&1)&&a[j][k])dp[s|1<<k][k]+=dp[s][j];
  I H=0;for(I x:dp[N-1])H+=x;
  I cycles=0,proper=0;std::vector<int> triangles;
  for(int first=0;first<n;++first){
    std::function<void(int,int)> visit=[&](int last,int used){
      int length=__builtin_popcount(static_cast<unsigned>(used));
      if(length>=3&&a[last][first]){++cycles;proper+=used!=N-1;if(length==3)triangles.push_back(used);}
      for(int j=first+1;j<n;++j)if(!(used>>j&1)&&a[last][j])visit(j,used|1<<j);
    };visit(first,1<<first);
  }
  bool pair=false;for(int s:triangles)for(int t:triangles)pair|=!(s&t);
  require(H-u>=cycles,"all-cycle defect");require(H-u-std::abs(v)>=proper,"proper-cycle defect");
  require((H-u==cycles)==!pair,"all-cycle defect equality");
  require((H-u-std::abs(v)==proper)==!pair,"proper-cycle defect equality");
  if(pair){require(H-u>=cycles+2,"two-packing all-cycle excess");require(H-u-std::abs(v)>=proper+2,"two-packing proper-cycle excess");++pack2;}
  ++checked;total_cycles+=cycles;total_proper+=proper;
}
int main(){
  for(int n=1;n<=6;++n)for(std::uint64_t c=0;c<(1ULL<<(n*(n-1)/2));++c)defects(coded(n,c));
  std::uint64_t seed=0x1371950C0FFEEULL;
  for(int n=7;n<=10;++n)for(int j=0;j<128;++j){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;defects(coded(n,seed));}
  std::cout<<"{\n  \"scope\": \"all labelled n1..6 and128 fixed-seed tournaments each n7..10\",\n  \"tournaments\": "<<checked<<",\n  \"two_disjoint_triangle_cases\": "<<pack2<<",\n  \"literal_directed_cycles\": "<<total_cycles<<",\n  \"literal_proper_cycles\": "<<total_proper<<",\n  \"verdict\": \"PASS both cycle-defect lower bounds, equality iff packing<=1, and extra>=2 for packing>=2\"\n}\n";
}
