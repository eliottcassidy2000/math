#include <vector>
#include <cstdint>
#include <iostream>
using namespace std;
using U=uint64_t;
using T=vector<unsigned>;
T join(const T&A,const T&B){
 T C; for(auto v:A)C.push_back(v|(((1u<<B.size())-1)<<A.size()));
 for(auto v:B)C.push_back(v<<A.size());
 return C;
}
T tt(int n){
 T a;for(int i=0;i<n;i++)a.push_back(((1u<<n)-1)^((1u<<(i+1))-1));
 return a;
}
U count(const T&A,bool cycles=false){
 int n=A.size(); unsigned d=(1u<<n)-1;
 vector<U> dp((d+1ull)*n);
 for(int i=0;i<n;i++)if(!cycles||i==n-1)dp[(1ull<<i)*n+i]=1;
 U total=0;
 for(unsigned m=1;m<=d;++m)for(int i=0;i<n;i++)if(m>>i&1){
  U z=dp[1ull*m*n+i]; if(!z)continue;
  for(unsigned k=A[i]&~m;k;k&=k-1){
   int j=__builtin_ctz(k);dp[(m|(1u<<j))*1ull*n+j]+=z;
  }
  if(m==d&&(!cycles||(A[i]>>(n-1)&1)))total+=z;
 }
 return total;
}
int width(int x){return 32-__builtin_clz(x);}
U build(int x,int y,const T&R,unsigned u){
 int p=width(x),q=width(y);
 T A=join(join(tt(p),R),tt(q));
 unsigned M=u<<p;
 for(int i=0;i<p;i++)M|=((x>>(p-1-i))&1)<<i;
 for(int i=0;i<q;i++)M|=(1-((y>>i)&1))<<(i+p+R.size());
 for(unsigned i=0;i<A.size();i++)if(!(M>>i&1))A[i]|=1u<<A.size();
 A.push_back(M);if(count(A,true)!=9)return 0;return count(A);
}
int main(){
 vector<T> R={{254,220,248,240,224,194,128,0},{254,252,184,240,224,192,132,0},{238,252,248,240,225,192,128,0}};
 vector<unsigned> U={182,202,169};
 for(int i=0;i<3;i++)if(build(18,30,R[i],U[i])!=374351+2*i)return 2;
 for(int i=0;i<2;i++)if(build(3,27,R[i],U[i])!=61817+2*i)return 3;
 const T c={10,28,25,16,1};
 if(build(78,78,c,12)!=517727||build(78,78,c,13)!=517729)return 4;
 if(build(2,2,c,15)!=575 || build(2,2,c,23)!=577 || build(2,2,c,14)!=579)return 5;
 cout<<"Ten literal Hamiltonian path counts and ten nine-cycle counts verified, including prime/7p pair in each rail.\n";
}
