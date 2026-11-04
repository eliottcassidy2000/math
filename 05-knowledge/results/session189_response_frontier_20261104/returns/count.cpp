#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
using U=unsigned long long;
int main(int argc,char**argv){
 if(argc!=2){std::cerr<<"usage: count m, 0<=m<=10\n";return 2;}
 int m=std::stoi(argv[1]);if(m<0||m>10){std::cerr<<"certified DP range is 0<=m<=10\n";return 2;}
 int n=m+8, w=4, A=m+5, B=m+6, C=m+7;
 std::vector<unsigned> out(n);
 out[0]=2;out[1]=4;out[2]=9;out[3]=3;
 for(int v=0;v<4;v++)if(v==3)out[v]|=1u<<w;else out[w]|=1u<<v;
 for(int z=5;z<m+5;z++){
   out[w]|=1u<<z;out[z]|=15;
   for(int z2=z+1;z2<m+5;z2++)out[z]|=1u<<z2;
 }
 for(int v=0;v<A;v++)if(v==w)out[v]|=1u<<A;else out[A]|=1u<<v;
 for(int v=0;v<B;v++)if(v==A)out[v]|=1u<<B;else out[B]|=1u<<v;
 for(int v=0;v<C;v++)if(v==B)out[v]|=1u<<C;else out[C]|=1u<<v;
 unsigned N=1u<<n,full=N-1;
 std::vector<U>end((size_t)N*n),start((size_t)N*n),h(N);h[0]=1;
 for(unsigned mask=1;mask<N;mask++){
   for(unsigned vv=mask;vv;vv&=vv-1){int v=__builtin_ctz(vv);unsigned rest=mask^(1u<<v);U e=0,s=0;
     if(!rest)e=s=1;
     else for(unsigned uu=rest;uu;uu&=uu-1){int u=__builtin_ctz(uu);if(out[u]>>v&1)e+=end[(size_t)rest*n+u];else s+=start[(size_t)rest*n+u];}
     end[(size_t)mask*n+v]=e;start[(size_t)mask*n+v]=s;h[mask]+=e;
   }
 }
 U p=0;std::vector<U>e2(n),s2(n),w2(n*n);
 for(unsigned mask=0;mask<N;mask++){unsigned other=full^mask;p+=h[mask]*h[other];
   for(int v=0;v<n;v++){e2[v]+=end[(size_t)mask*n+v]*h[other];s2[v]+=start[(size_t)mask*n+v]*h[other];}
   for(int v: {0,1,2,4})for(int u:{0,1,2,4})w2[v*n+u]+=end[(size_t)mask*n+v]*start[(size_t)other*n+u];
 }
 std::cout<<"{\"m\":"<<m<<",\"rows\":[";for(int v=0;v<n;v++)std::cout<<(v?",":"")<<out[v];std::cout<<"],\"starts\":[";for(int v=0;v<n;v++)std::cout<<(v?",":"")<<start[(size_t)full*n+v];
 std::cout<<"],\"fields\":[";bool first=true;
 for(int v:{0,1,2}){std::vector<U>f={h[full],p,end[(size_t)full*n+v],start[(size_t)full*n+v],h[full^(1u<<v)],e2[v],s2[v],end[(size_t)full*n+w],start[(size_t)full*n+w],h[full^(1u<<w)],e2[w],s2[w],w2[v*n+v],w2[v*n+w],w2[w*n+v],w2[w*n+w]};for(U z:f){std::cout<<(first?"":",")<<z;first=false;}}
 std::cout<<"]}\n";
}
