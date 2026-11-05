#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using I=long long;using V=std::array<I,5>;
void need(bool p){if(!p)throw std::runtime_error("independent full-gadget inverse failed");}
std::array<I,3> physical_coeff(I L){
 I a=L-1;int m=0;for(I x=a;x;x>>=1)++m;
 I d[3][3]={{1,0,0},{0,0,0},{0,0,0}};
 for(int i=0;i<m;++i){I z[3][3]={};bool toRoot=(a>>i)&1;
  for(int l=0;l<3;++l)for(int r=0;r<3;++r){z[toRoot?1:2][r]+=d[l][r];z[l][r?r:(toRoot?2:1)]+=d[l][r];}
  for(int l=0;l<3;++l)for(int r=0;r<3;++r)d[l][r]=z[l][r];
 }
 I g=d[0][0]+d[0][1]+d[1][0]+d[1][1],Q=1;
 for(int i=0;i<m;++i)if(!(a>>i&1))Q+=1LL<<(m-1-i);
 return{g,Q,1LL<<m};
}
I value(V v,I L){auto c=physical_coeff(L);return v[0]*c[0]+v[1]*L+v[2]*c[1]+v[3]*c[2]-v[4]*(L&1);}
I rev(I x,int k){I z=0;for(int i=0;i<k;++i)if(!(x>>i&1))z|=1LL<<(k-1-i);return z;}
std::vector<I> inverse(V v,I T){
 std::vector<I>out;if(T<1)return out;I K=T/v[1];
 if(value(v,1)==T)out.push_back(1);
 if(T<value(v,2))return out;
 for(int k=1;(1LL<<(k-1))+1<=K;++k){I Z=1LL<<k,lo=Z/2+1,hi=std::min(Z,K);
  if(k<=4){for(I L=lo;L<=hi;++L)if(value(v,L)==T)out.push_back(L);continue;}
  // Choose high width by the inequality, rather than the source's q formula.
  I H=1;while(8*(2*H)*(2*H)<=Z)H*=2;I low=Z/H;
  for(I b=0;b<low;++b){I xlo=std::max<I>(0,(lo-b-1+low-1)/low),xhi=(hi-b-1)/low;
   if(xlo>xhi)continue;
   while(xlo<xhi){I mid=(xlo+xhi)/2;if(value(v,low*mid+b+1)>=T)xhi=mid;else xlo=mid+1;}
   I L=low*xlo+b+1;if(value(v,L)==T)out.push_back(L);
  }
 }
 std::sort(out.begin(),out.end());return out;
}
std::vector<I> gadget(I L){I a=L-1;int m=0;for(I x=a;x;x>>=1)++m;std::vector<I>T(m+1);
 for(int i=0;i<m;++i){int v=i+1;if(a>>i&1)T[v]|=1;else T[0]|=1LL<<v;for(int j=i+1;j<m;++j)T[v]|=1LL<<(j+1);}
 return T;
}
std::array<I,2> literal(const std::vector<I>&T){int n=T.size();std::vector<int>p(n);std::iota(p.begin(),p.end(),0);I h=0,R=0;
 do{int bad=0,r=-1;for(int i=0;i<n;++i){if(p[i]==0)r=i;if(i&&!(T[p[i-1]]>>p[i]&1))++bad;}
  h+=bad==0;int gapbad=r>0&&!(T[p[r-1]]&1);R+=bad==gapbad;
 }while(std::next_permutation(p.begin(),p.end()));return{h,R};
}
int main(){I physical=0,identities=0,monotone=0,inverses=0,shells=0;
 for(I L=1;L<=128;++L){auto a=literal(gadget(L));auto c=physical_coeff(L);need(a[0]==c[0]&&a[1]==c[1]);++physical;}
 for(I L=1;L<=32768;++L){int k=0;for(I a=L-1;a;a>>=1)++k;auto c=physical_coeff(L);need(c[1]==1+rev(L-1,k)&&c[2]==(1LL<<k));++identities;}
 for(int k=5;k<=13;++k){int q=(k+4)/2,p=k-q;I low=1LL<<q,H=1LL<<p;
  for(I b=0;b<low;++b){I u=rev(b,q);
   for(I x=H/2;x<H;++x){I L=low*x+b+1;auto c=physical_coeff(L);
    need(c[1]==1+H*u+rev(x,p));need(c[0]==physical_coeff(x+1)[0]+u*x+(u&1)*(x&1)+physical_coeff(b+1)[0]-1);++shells;}
   for(I A:{1,3})for(I B:{A,A+1})for(I C:{0LL,4*A})for(I E:{0LL,A}){
    V v={A,B,C,2,E};for(I x=H/2;x+1<H;++x){need(value(v,low*x+b+1)<value(v,low*(x+1)+b+1));++monotone;}
   }
  }
 }
 for(I A=1;A<=4;++A)for(I B:{A,A+1,2*A})for(I C:{0LL,A,4*A})for(I D:{0,2})for(I E:{0LL,A/2,A}){
  V v={A,B,C,D,E};std::vector<std::vector<I>>hits(513);
  for(I L=1;L<=512/B;++L){I T=value(v,L);if(T<=512)hits[T].push_back(L);}
  for(I T=0;T<=512;++T){need(inverse(v,T)==hits[T]);++inverses;}
 }
 V big[]={V{9,9,30,30,3},V{39,54,156,216,6},V{100,300,400,1234,100}};
 I Ls[]={1234567,9876543,100000007},Ts[]={119881395,5703260151,235297253052};
 for(int i=0;i<3;++i){need(value(big[i],Ls[i])==Ts[i]);need(inverse(big[i],Ts[i])==std::vector<I>{Ls[i]});}
 I boundary=0;for(V v:{V{1,1,4,0,1},V{3,4,12,2,0},V{7,14,0,5,7}})for(int k=1;k<=18;++k)for(I L:{(1LL<<k)-1,1LL<<k,(1LL<<k)+1}){
  I T=value(v,L);auto out=inverse(v,T);need(std::find(out.begin(),out.end(),L)!=out.end());for(I x:out)need(value(v,x)==T);++boundary;
 }
 std::cout<<"{\"status\":\"PASS\",\"literal_gadget_profiles\":"<<physical<<",\"minimal_width_checks\":"<<identities<<",\"shell_identity_pairs\":"<<shells<<",\"strict_branch_comparisons\":"<<monotone<<",\"exhaustive_inverse_cases\":"<<inverses<<",\"dyadic_boundary_inverses\":"<<boundary<<",\"large_unique_inverses\":3}\n";
}
