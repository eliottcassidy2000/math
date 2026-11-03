// Cheap local test of the incoming nonuniform triangle-replacement interface.
// Keep boundary vertex0, delete another root, toggle at most one old incidence
// at each of three cyclic replacement ports. All retained old arcs are fixed.
#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using I=long long;using V=std::array<I,3>;
std::pair<I,std::vector<I>> dp(const std::vector<unsigned>&r){
 int n=r.size();unsigned full=(1u<<n)-1;std::vector<unsigned>in(n);
 for(int u=0;u<n;++u)for(int v=0;v<n;++v)if(r[u]>>v&1)in[v]|=1u<<u;
 std::vector<I>a((full+1)*n);for(int u=0;u<n;++u)a[(1u<<u)*n+u]=1;
 for(unsigned m=1;m<=full;++m){unsigned t=m;while(t){int u=__builtin_ctz(t);t&=t-1;unsigned prev=m^(1u<<u),x=prev&in[u];
  while(x){int v=__builtin_ctz(x);x&=x-1;a[m*n+u]+=a[prev*n+v];}}}
 I H=0;for(int u=0;u<n;++u)H+=a[full*n+u];return {H,a};
}
V state(const std::vector<unsigned>&r,I H,const std::vector<I>&a){
 int n=r.size();unsigned full=(1u<<n)-1;std::vector<unsigned>rev(n);
 for(int u=0;u<n;++u)for(int v=0;v<n;++v)if(r[u]>>v&1)rev[v]|=1u<<u;
 auto backward=dp(rev).second;std::vector<I>counts(full+1);counts[0]=1;
 for(unsigned m=1;m<=full;++m)for(int u=0;u<n;++u)counts[m]+=a[m*n+u];
 I E=a[full*n],C=0;for(unsigned m=1;m<=full;m+=2)C+=backward[m*n]*counts[full^m];
 return {H,(C-E)/2,(C+E)/2};
}
bool strong(const std::vector<unsigned>&r){unsigned full=(1u<<r.size())-1;
 for(unsigned u=0;u<r.size();++u){unsigned s=1u<<u,old=0;while(s!=old){old=s;for(unsigned v=0;v<r.size();++v)if(s>>v&1)s|=r[v];}if(s!=full)return false;}return true;}
int main(int argc,char**argv){
 if(argc!=3)return 2;std::ifstream in(argv[1]);std::ofstream out(argv[2]);int n;I h,x,y;
 while(in>>n>>h>>x>>y){std::vector<unsigned>T(n);for(auto&v:T)in>>v;
  int tested=0,hits[2]={0,0},residue[2]={0,0},exact[2]={0,0};
  for(int root=1;root<n;++root){std::vector<int>old;for(int u=0;u<n;++u)if(u!=root)old.push_back(u);int m=n-1,N=m+3;unsigned base=0;
   for(int j=0;j<m;++j)if(T[root]>>old[j]&1)base|=1u<<j;
   std::vector<unsigned>options{base};for(int j=0;j<m;++j)options.push_back(base^(1u<<j));
   for(auto a:options)for(auto b:options)for(auto c:options){std::array<unsigned,3>ports{a,b,c};
    if(ports>std::array<unsigned,3>{b,c,a}||ports>std::array<unsigned,3>{c,a,b})continue;
    ++tested;std::vector<unsigned>R(N);
    for(int i=0;i<m;++i)for(int j=0;j<m;++j)if(T[old[i]]>>old[j]&1)R[i]|=1u<<j;
    for(int k=0;k<3;++k){R[m+k]|=1u<<(m+(k+1)%3);for(int j=0;j<m;++j)if(ports[k]>>j&1)R[m+k]|=1u<<j;else R[j]|=1u<<(m+k);}
    auto answer=dp(R);int digit=answer.first==3*h+2?1:(answer.first==3*h+4?2:0);if(!digit)continue;
    if(!strong(R))continue;++hits[digit-1];V v=state(R,answer.first,answer.second),ideal{3*h+2*digit,3*x-digit,3*y};
    bool mod=true;for(int j=0;j<3;++j)mod&=(v[j]-ideal[j])%7==0;
    if(mod)++residue[digit-1];if(v==ideal)++exact[digit-1];
    if(mod){out<<h<<' '<<root<<' '<<a<<' '<<b<<' '<<c;for(I z:v)out<<' '<<z;out<<' '<<N;for(auto z:R)out<<' '<<z;out<<'\n';}
   }
  }
  std::cout<<"H "<<h<<" tested "<<tested<<" Hhits "<<hits[0]<<' '<<hits[1]<<" residue "<<residue[0]<<' '<<residue[1]<<" exact "<<exact[0]<<' '<<exact[1]<<'\n';
 }
}
