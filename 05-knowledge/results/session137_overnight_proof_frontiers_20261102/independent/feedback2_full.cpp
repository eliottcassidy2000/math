// Independent full-matrix determinant audit, never Schur/row-difference.
// For n<=10, Hadamard gives |det(A)|, |det(I+A)| <= 10^5.
// Thus modular Gaussian elimination at prime 1000003, reconstructed into
// (-p/2,p/2), is an exact integer determinant, not a probabilistic check.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using Matrix=std::vector<std::vector<long long>>;
constexpr long long prime=1000003;
long long power(long long a,long long e){long long b=1;for(;e;e>>=1,a=a*a%prime)if(e&1)b=b*a%prime;return b;}
long long determinant(Matrix a){
  const int n=a.size();long long value=1;
  for(int j=0;j<n;++j){
    int p=j;while(p<n&&a[p][j]==0)++p;
    if(p==n)return 0;
    if(p!=j){std::swap(a[p],a[j]);value=prime-value;}
    auto pivot=a[j][j];value=value*pivot%prime;auto inverse=power(pivot,prime-2);
    for(int i=j+1;i<n;++i){
      auto f=a[i][j]*inverse%prime;
      for(int k=j;k<n;++k)a[i][k]=(a[i][k]-f*a[j][k]%prime+prime)%prime;
    }
  }
  return value>prime/2?value-prime:value;
}
bool strong(const Matrix&a){
  const int n=a.size();
  for(int s=0;s<n;++s){
    unsigned reached=1u<<s;bool change=true;
    while(change){change=false;for(int i=0;i<n;++i)if(reached>>i&1u)
      for(int j=0;j<n;++j)if(a[i][j]&&!(reached>>j&1u)){reached|=1u<<j;change=true;}}
    if(reached!=(1u<<n)-1)return false;
  }
  return true;
}
void require(bool p,const char*m){if(!p)throw std::runtime_error(m);}
int main(){
  std::array<int,8> expected_counts={1,6,29,128,540,2224,9040,36480};
  std::array<int,8> expected_max={3,3,8,10,16,20,26,30};
  std::cout<<"{\n  \"method\": \"exact full-matrix modular Gaussian determinants, Hadamard reconstruction\",\n  \"prime\": 1000003,\n  \"rows\": [\n";
  for(int m=1;m<=8;++m){
    int n=m+2,count=0;long long maximum=0;unsigned witness=0;
    for(unsigned word=0;word<(1u<<(2*m));++word){
      Matrix a(n,std::vector<long long>(n));
      for(int i=0;i<m;++i)for(int j=i+1;j<m;++j)a[i][j]=1;
      a[m][m+1]=1;
      for(int i=0;i<m;++i)for(int j=0;j<2;++j){a[i][m+j]=(word>>(2*i+j))&1u;a[m+j][i]=1-a[i][m+j];}
      Matrix b=a;for(int i=0;i<n;++i)b[i][i]=1;
      auto u=determinant(b),v=determinant(a)*(n%2?-1:1),P=u+std::abs(v);
      require(u>0,"u positive");require(u+v>0,"disc positive");
      require(u<=(1+m/2)*(1+m/2)+(m+1)*(m+1),"u bound");
      require(std::abs(v)<=4*(m+1),"v bound");
      if(!strong(a))continue;
      ++count;if(P>maximum){maximum=P;witness=word;}
    }
    require(count==expected_counts[m-1],"strong universe");require(maximum==expected_max[m-1],"strong max P");
    std::cout<<"    {\"n\": "<<n<<", \"presentations\": "<<(1u<<(2*m))<<", \"strong\": "<<count<<", \"max_P\": "<<maximum<<", \"max_word\": "<<witness<<"}"<<(m==8?"\n":",\n");
  }
  std::cout<<"  ]\n}\n";
}
