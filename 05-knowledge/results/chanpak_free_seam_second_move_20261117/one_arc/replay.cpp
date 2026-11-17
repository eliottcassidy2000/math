#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
const unsigned D[]={54,32,2,119,6,20,55};
vector<int> bits(unsigned x){
 vector<int> out;
 for(int j=bit_width(x)-1;j>=0;--j)out.push_back((x>>j)&1);
 return out;
}
vector<unsigned> make(unsigned u,unsigned w){
 auto l=bits(u),r=bits(w);
 unsigned n=l.size()+7+r.size(), left=l.size(), right=left+7, mark=0;
 vector<unsigned>a(n+1,0);
 for(unsigned j=0;j<left;j++)mark|=l[j]<<j;
 mark|=26u<<left;
 for(unsigned j=0;j<r.size();j++)mark|=(1-r[r.size()-1-j])<<(right+j);
 for(unsigned i=0;i<n;i++){
  for(unsigned j=i+1;j<n;j++){
   bool forward=i<left||j>=right||(D[i-left]>>(j-left)&1u);
   if(forward)a[i]|=1u<<j;else a[j]|=1u<<i;
  }
  if(!(mark>>i&1u))a[i]|=1u<<n;
 }
 a[n]=mark;
 return a;
}
uint64_t count(const vector<unsigned>&a){
 const unsigned n=a.size(),full=(1u<<n)-1;
 for(unsigned i=0;i<n;i++){
  unsigned s=1u<<i;
  for(unsigned t=0;t<n;t++)for(unsigned j=0;j<n;j++)if(s>>j&1u)s|=a[j];
  if(s!=full)throw runtime_error("not strong");
 }
 vector<uint64_t>d(size_t(n)*(full+1),0);
 for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;
 for(unsigned m=1;m<=full;m++){
  for(unsigned s=m;s;s&=s-1){
   unsigned i=countr_zero(s); auto z=d[size_t(m)*n+i];
   if(z)for(unsigned next=a[i]&~m;next;next&=next-1){
    unsigned j=countr_zero(next);d[size_t(m|(1u<<j))*n+j]+=z;
   }
  }
 }
 uint64_t v=0;
 for(unsigned i=0;i<n;i++)v+=d[size_t(full)*n+i];
 return v;
}
void check_single(const vector<unsigned>&a,const vector<unsigned>&b){
 unsigned delta=0;
 for(unsigned i=0;i<a.size();i++)delta+=popcount(a[i]^b[i]);
 if(delta!=2)throw runtime_error("not a single arc");
}
int main(){
 auto a=make(1,650),b=make(650,1),c=make(1,648),d=make(648,1);
 check_single(a,c);check_single(b,d);
 auto ar=array{count(a),count(b),count(c),count(d)};
 if(ar!=array<uint64_t,4>{261163,261165,261167,261169})throw runtime_error("count");
 for(auto x:ar)cout<<x<<" ";
 cout<<endl;
}
