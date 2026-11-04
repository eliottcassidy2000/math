#define main inherited_search_main
#include "search.cpp"
#undef main
struct Table{vector<I>e,h;};
Table table(const T&t){int n=t.size(),f=(1<<n)-1;Table d{vector<I>((f+1)*n),vector<I>(f+1)};d.h[0]=1;
 for(int i=0;i<n;i++)d.e[(1<<i)*n+i]=1;
 for(int m=1;m<=f;m++)for(unsigned z=m;z;z&=z-1){int i=countr_zero(z);I a=d.e[m*n+i];d.h[m]+=a;for(unsigned q=t[i]&~m;q;q&=q-1)d.e[(m|(q&-q))*n+countr_zero(q)]+=a;}return d;}
vector<I> births(const T&t){int n=t.size(),f=(1<<n)-1;T rev=t;for(int i=0;i<n;i++)rev[i]=f^t[i]^(1<<i);auto a=table(t),b=table(rev);vector<I>w(n*n),l(n),out(f+1,a.h[f]);
 for(int m=1;m<f;m++)for(unsigned z=m;z;z&=z-1){int i=countr_zero(z);if(!a.e[m*n+i])continue;for(unsigned q=f^m;q;q&=q-1){int j=countr_zero(q);w[i*n+j]+=a.e[m*n+i]*b.e[(f^m)*n+j];}}
 for(int i=0;i<n;i++){l[i]=a.e[f*n+i]-b.e[f*n+i];for(int j=0;j<n;j++)l[i]+=w[i*n+j];}
 for(int m=1;m<=f;m++){int i=countr_zero((unsigned)m),r=m^(1<<i);out[m]=out[r]+l[i];for(unsigned z=r;z;z&=z-1){int j=countr_zero(z);out[m]-=w[i*n+j]+w[j*n+i];}}return out;}
T birth(T t,int mask){int n=t.size();for(int i=0;i<n;i++)if(mask>>i&1)t[i]|=1<<n;t.push_back(((1<<n)-1)^mask);return t;}
#ifndef BIRTH_LIBRARY
int main(int argc,char**argv){int n;I target;int degree=argc>1?stoi(argv[1]):2;while(cin>>n>>target){if(n>17)return 2;T t(n);for(auto&x:t)cin>>x;auto first=births(t);int eligible=0,pruned=0,checked=0,hits=0;cout<<"{\"degree\":"<<degree<<",\"target\":"<<target<<",\"hits\":[";
 for(int m=0;m<(1<<n);m++){int d=popcount((unsigned)m);if(degree>=0&&d!=degree&&n-d!=degree)continue;eligible++;if(first[m]>target){pruned++;continue;}auto u=birth(t,m);auto second=births(u);checked++;
  for(int c=0;c<(1<<(n+1));c++)if(second[c]==target){auto v=birth(u,c);if(hits++)cout<<',';cout<<"{\"first\":"<<m<<",\"second\":"<<c<<",\"intermediate_H\":"<<first[m]<<",\"strong\":"<<(strong(v)?"true":"false")<<",\"rows\":";print(v);cout<<'}'<<flush;}}
 cout<<"],\"eligible\":"<<eligible<<",\"pruned\":"<<pruned<<",\"checked\":"<<checked<<"}\n"<<flush;}return 0;}
#endif
