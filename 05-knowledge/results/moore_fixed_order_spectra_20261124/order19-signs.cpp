#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;
constexpr int M=19,Mask=(1<<M)-1;
using Poly=array<int,18>;
int rot(int m,int k){return ((m<<k)|(m>>(M-k)))&Mask;}
int flip(int m){int v=m&1;for(int k=1;k<M;k++)if(m>>k&1)v|=1<<(M-k);return v;}
int norm(int v){v%=3;return v<0?v+3:v;}
int sg(int m,int k){return (m>>k&1)?1:-1;}
int sum(int m){return 2*popcount((unsigned)m)-M;}
int corr(int m,int n,int k){return M-2*popcount((unsigned)(m^rot(n,k)));}
int signature(int m){
 int s=0;for(int k=1;k<=9;k++)s=3*s+norm(corr(m,m,k));
 return s;
}
int canon(int m){int u=m;for(int k=0;k<M;k++){int t=rot(m,k);u=min({u,t,t^Mask});}return u;}
Poly reduce(int m){
 Poly p{};int e=sg(m,18);for(int z=0;z<18;z++)p[z]=norm(sg(m,z)-e);return p;
}
int code(const Poly& a){int s=0;for(int k=17;k>=0;k--)s=3*s+a[k];return s;}
Poly mul(const Poly& a,const Poly& b){
 array<int,19> c{};
 for(int i=0;i<18;i++)if(a[i])for(int j=0;j<18;j++)c[(i+j)%19]+=a[i]*b[j];
 Poly s{};for(int i=0;i<18;i++)s[i]=norm(c[i]-c[18]);return s;
}
Poly inv(const Poly& a){
 Poly R{};R[0]=1;Poly x=a;
 for(int n=387420489-2;n;n/=2){if(n&1)R=mul(R,x);x=mul(x,x);}
 return R;
}
int margin(int g,int f,int h,int e){
 int E[3][3]={{g,f,h},{flip(f),g,e},{flip(h),flip(e),g}};
 int dlow=100,dhigh=-1;
 for(int i=0;i<3;i++)for(int j=i;j<3;j++)for(int k=(i==j);k<M;k++){
    int z=0;for(int q=0;q<3;q++)z+=corr(E[i][q],E[j][q],k);
    const int a=z+3*sg(E[i][j],k);assert(a%6==0);
    int d=(56-a)/2;
    dlow=min(dlow,d);dhigh=max(dhigh,d);

 }
 return min(dlow,56-dhigh);
}
int main(){
 vector<vector<int>> rows(19683);
 unordered_map<int,int> decode;
 for(int v=1;v<Mask;v++){
  assert(decode.emplace(code(reduce(v)),v).second);
  if(canon(v)==v)rows[signature(v)].push_back(v);
 }
 long long examined=0,signs=0,formula=0,examinedG=0;
 array<int,29> hist{};
 for(int c=0;c<512;c++){
  int g=0;for(int i=0;i<9;i++)if(c>>i&1)g|=(1<<(i+1))|(1<<(18-i));
  if(!g)continue;
  examinedG++;
  const int G=norm(sum(g));
  Poly ig=inv(reduce(g));
  for(int f:rows[signature(g)]){
   auto A=mul(reduce(flip(f)),ig);
   for(int h:rows[signature(g)])if(h>=f){
    examined++;
    auto ei=decode.find(code(mul(A,reduce(h))));
    if(ei==decode.end())continue;
    const int e=ei->second;
    if((!G&&norm(sum(e))!=0)||norm(G*sum(e)-sum(f)*sum(h))!=0)continue;
    formula++;
    assert(canon(f)==canon(g)||canon(h)==canon(g)||canon(e)==canon(g));
    const int t=margin(g,f,h,e);
    assert(t>=0&&t<=28);hist[t]++;if(t>=10)signs++;
   }
  }
 }
 for(int i=0;i<29;i++)if(hist[i])cout<<"minimum complementary swap count "<<i<<": "<<hist[i]<<"\n";
 cout<<"rows "<<examinedG<<" pairs "<<examined<<" affine "<<formula<<" valid "<<signs<<endl;
}