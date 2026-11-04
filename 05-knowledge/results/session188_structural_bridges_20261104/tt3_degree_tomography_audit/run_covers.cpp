#include <iostream>
#include <vector>
#include <cstdint>
using U=std::uint64_t;
// Quotient labels z=0, X=1, Y=2, then individually named exterior vertices.
std::vector<U> run(const std::vector<unsigned>&Q,unsigned side,bool zfirst){
 unsigned n=Q.size(),b=n,rad=b+1,N=1u<<(n-2);auto pos=[&](unsigned M,unsigned x,unsigned y,unsigned ph,unsigned last){return ((((std::size_t(M)*rad+x)*rad+y)*2+ph)*n+last);};std::vector<U>d(std::size_t(N)*rad*rad*2*n),out(rad*rad);
 auto bit=[](unsigned v){return v==0?1u:(v>=3?1u<<(v-2):0u);};
 for(unsigned v=0;v<n;++v){if(!zfirst&&v==0)continue;d[pos(bit(v),v==1,v==2,0,v)]=1;}
 for(unsigned M=0;M<N;++M)for(unsigned x=0;x<=b;++x)for(unsigned y=0;y<=b;++y)for(unsigned ph=0;ph<2;++ph)for(unsigned last=0;last<n;++last){U val=d[pos(M,x,y,ph,last)];if(!val)continue;if(M==N-1&&x&&y&&ph)out[x*rad+y]+=val;
  for(unsigned v=0;v<n;++v){unsigned B=bit(v);if(B&&(M&B))continue;if(v==1&&x==b)continue;if(v==2&&y==b)continue;unsigned pp=ph;bool allowed=false;
   if(!ph){if(zfirst){if(last==0){allowed=v==side;pp=1;}else allowed=last!=v&&(Q[last]>>v&1);}else{if(v==0){allowed=last==side;pp=1;}else allowed=last!=v&&(Q[last]>>v&1);}}
   else allowed=last!=v&&(Q[last]>>v&1);
   if(allowed)d[pos(M|B,x+(v==1),y+(v==2),pp,v)]+=val;
  }
 }
 return out;
}
int main(){unsigned cases;std::cin>>cases;while(cases--){unsigned n;std::cin>>n;std::vector<unsigned>Q(n);for(auto&r:Q)std::cin>>r;auto a=run(Q,1,true),b=run(Q,1,false),c=run(Q,2,true),d=run(Q,2,false);unsigned rad=n+1;for(unsigned x=1;x<=n;++x)for(unsigned y=1;y<=n;++y)std::cout<<x<<' '<<y<<' '<<a[x*rad+y]<<' '<<b[x*rad+y]<<' '<<c[x*rad+y]<<' '<<d[x*rad+y]<<'\n';}}
