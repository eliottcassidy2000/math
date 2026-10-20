#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
struct Result {int64_t h;vector<int64_t>d;};
Result evaluate(int k,uint64_t A,uint64_t b,bool grad){int n=2*k,N=1<<k,adj[24]={},z=0;auto edge=[&](int u,int v,int yes){if(yes)adj[u]|=1<<v;else adj[v]|=1<<u;};for(int i=0;i<k;i++)for(int j=i+1;j<k;j++,z++){int a=A>>z&1,c=b>>z&1;edge(2*i,2*j,a);edge(2*j+1,2*i+1,a);edge(2*i,2*j+1,c);edge(2*j,2*i+1,c);}static int64_t f[4096][24],g[4096][24];fill(&f[0][0],&f[N-1][24],0);if(grad)fill(&g[0][0],&g[N-1][24],0);for(int i=0;i<n;i++){f[1<<(i/2)][i]=1;if(grad)g[1<<(i/2)][i]=1-(i&1);}
 for(int mask=1;mask<N;mask++){int av=0;for(int j=0;j<k;j++)if(!(mask>>j&1))av|=3<<(2*j);for(int i=0;i<n;i++){auto v=f[mask][i];if(v){int next=adj[i]&av;while(next){int j=__builtin_ctz(next);next&=next-1;f[mask|(1<<(j/2))][j]+=v;}}}if(grad)for(int i=0;i<n;i++)if(mask>>i/2&1){int sub=mask^(1<<i/2);if(!sub)continue;int next=adj[i]&~av&~(3<<(2*(i/2)));while(next){int j=__builtin_ctz(next);next&=next-1;g[mask][i]+=g[sub][j];}}}
 Result R{0,{}};for(int i=0;i<n;i+=2)R.h+=f[N-1][i];if(!grad)return R;int64_t w[24][24]={};for(int mask=1;mask<N-1;mask++){int sub=(N-1)^mask;for(int i=0;i<n;i++)if(f[mask][i])for(int j=0;j<n;j++)if(g[sub][j])w[i][j]+=f[mask][i]*g[sub][j];}R.d.resize(2*z);z=0;int m=k*(k-1)/2;for(int i=0;i<k;i++)for(int j=i+1;j<k;j++,z++){auto da=w[2*j][2*i]-w[2*i][2*j]+w[2*i+1][2*j+1]-w[2*j+1][2*i+1];auto db=w[2*j+1][2*i]-w[2*i][2*j+1]+w[2*i+1][2*j]-w[2*j][2*i+1];R.d[z]=(A>>z&1)?da:-da;R.d[m+z]=(b>>z&1)?db:-db;}return R;}
int main(int argc,char**argv){int k=argc>1?stoi(argv[1]):8,restarts=argc>2?stoi(argv[2]):1000;if(k<2||k>10)return 2;int m=k*(k-1)/2;uint64_t mask=(1ULL<<m)-1,bestA=0,bestb=0;int64_t best=0;mt19937_64 rng(20261020);auto start=chrono::steady_clock::now();uint64_t evals=0;
 for(int r=0;r<restarts;r++){uint64_t A=rng()&mask,b=rng()&mask;if(r%2&&best){A=bestA;b=bestb;int perturb=2+rng()%(2*k);for(int t=0;t<perturb;t++){int q=rng()%(2*m);if(q<m)A^=1ULL<<q;else b^=1ULL<<(q-m);}}
  for(int step=0;step<100;step++){auto R=evaluate(k,A,b,true);evals++;if(r<2&&step<2)for(int q=0;q<2*m;q++){uint64_t a=A,c=b;if(q<m)a^=1ULL<<q;else c^=1ULL<<(q-m);if(evaluate(k,a,c,false).h!=R.h+R.d[q]){cerr<<"gradient mismatch"<<endl;return 3;}}
   if(R.h>best){best=R.h;bestA=A;bestb=b;cerr<<"restart="<<r<<" step="<<step<<" best="<<best<<" A="<<A<<" b="<<b<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
   auto mx=*max_element(R.d.begin(),R.d.end());if(mx<=0)break;vector<int>qs;for(int q=0;q<2*m;q++)if(R.d[q]==mx)qs.push_back(q);int q=qs[rng()%qs.size()];if(q<m)A^=1ULL<<q;else b^=1ULL<<(q-m);
  }
 }
 cout<<"{\"k\":"<<k<<",\"restarts\":"<<restarts<<",\"evaluations\":"<<evals<<",\"max_h_found\":"<<best<<",\"A\":"<<bestA<<",\"b\":"<<bestb<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<",\"global_optimality_claimed\":false}\n";
}
