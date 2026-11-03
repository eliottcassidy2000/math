// Exact subset/endpoint counts of supplied one/two-start parent grammars.
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using I=int64_t;
int main(int argc,char**argv){
 I target=argc>1?stoll(argv[1]):6325;int n,b,index;long long tested=0;
 while(cin>>n>>b>>index){
  vector<int>T(n),pred(n);for(int&i:T)cin>>i;int f=(1<<n)-1;
  for(int i=0;i<n;i++)pred[i]=f^T[i]^(1<<i);
  vector<I>h(f+1),e((f+1)*n);h[0]=1;
  for(int s=1;s<=f;s++)for(int bits=s;bits;bits&=bits-1){
   int v=countr_zero((unsigned)bits),r=s^(1<<v);I z=!r;
   for(int js=r&pred[v];js;js&=js-1)z+=e[r*n+countr_zero((unsigned)js)];
   e[s*n+v]=z;h[s]+=z;
  }
  I F=0;for(int s=1;s<=f;s+=2)F+=e[s*n]*h[f^s];tested++;
  if(h[f]+F-b==target){
   cout<<"{\"parent_order\":"<<n<<",\"budget\":"<<b<<",\"rank\":"<<index<<",\"h\":"<<h[f]<<",\"F\":"<<F<<",\"rows\":[";
   for(int i=0;i<n;i++){if(i)cout<<',';cout<<T[i];}cout<<"]}\n";cout.flush();
  }
 }
 cerr<<"tested "<<tested<<'\n';
}
