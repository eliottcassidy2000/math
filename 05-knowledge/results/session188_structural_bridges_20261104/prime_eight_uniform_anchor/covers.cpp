#include <iostream>
#include <vector>
#include <cstdint>
std::vector<std::vector<unsigned>> solve(const std::vector<unsigned>&Q,unsigned a,unsigned c){
 unsigned n=Q.size(),states=1u<<(2*n),lastN=n+2,all=(1u<<n)-1;std::vector<std::uint8_t>d(std::size_t(states)*4*lastN*3);auto at=[&](unsigned s,unsigned p,unsigned l,unsigned h)->std::uint8_t&{return d[((std::size_t(s)*4+p)*lastN+l)*3+h];};
 for(unsigned i=0;i<n;++i){at(1u<<(2*i),0,i,0)=1;}at(0,1,n,0)=1;at(0,2,n+1,0)=1;
 std::vector<std::vector<unsigned>>ans(3,std::vector<unsigned>(states));
 for(unsigned s=0;s<states;++s){unsigned seen=0,repeated=0;for(unsigned i=0;i<n;++i){unsigned z=(s>>(2*i))&3;if(z)seen|=1u<<i;if(z>=2)++repeated;}if(repeated>=3)continue;
  for(unsigned p=0;p<4;++p)for(unsigned l=0;l<lastN;++l)for(unsigned h=0;h<3;++h){unsigned z=at(s,p,l,h);if(!z)continue;if(seen==all&&p==3)ans[h][s]=(ans[h][s]+z)%8;
   for(unsigned j=0;j<lastN;++j){unsigned ss=s,pp=p;bool edge;
    if(j<n){unsigned digit=(s>>(2*j))&3;if(digit==3||(digit==1&&repeated==2))continue;ss+=1u<<(2*j);edge=j!=l&&(l<n?(Q[l]>>j&1):!((l==n?a:c)>>j&1));}
    else{unsigned bit=1u<<(j-n);if(p&bit)continue;pp|=bit;edge=l<n?((j==n?a:c)>>l&1):(l==n&&j==n+1);}
    if(edge){auto&v=at(ss,pp,j,h);v=(v+z)&7;}
    if(h<2){auto&v=at(ss,pp,j,h+1);v=(v+z)&7;}
   }
  }
 }
 return ans;
}
int main(){std::vector<unsigned>Q={238,4,240,118,195,146,34,74},U(8);for(unsigned i=0;i<8;++i)U[i]=255^Q[i]^(1u<<i);for(auto query:std::vector<std::pair<unsigned,unsigned>>{{1,0},{255,254}}){auto x=solve(Q,query.first,query.second),y=solve(U,query.first,query.second);for(unsigned s=0;s<65536;++s){unsigned seen=0,r=0;for(unsigned i=0;i<8;++i){auto z=(s>>(2*i))&3;if(z)++seen;if(z>1)++r;}if(seen==8&&r<=2){std::cout<<query.first<<' '<<query.second<<' '<<s;for(unsigned h=0;h<3;++h)std::cout<<' '<<x[h][s]<<' '<<y[h][s];std::cout<<'\n';}}}}
