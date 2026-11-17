#include <bits/stdc++.h>
using namespace std;
int main(){
 string s;size_t tested=0;
 while(getline(cin,s)){
 istringstream ss(s);unsigned n,target;
 if(!(ss>>n>>target)||n<1||n>21)continue;
 vector<unsigned>a(n);for(auto &v:a)if(!(ss>>v))return 2;
 vector<unsigned>d((1<<n)*n);for(unsigned i=0;i<n;i++)d[(1<<i)*n+i]=1;
 for(unsigned m=1;m<(1u<<n);m++)for(unsigned t=m;t;){
  unsigned b=t&-t; t-=b; unsigned i=__builtin_ctz(b),p=d[m*n+i];
  if(p)for(unsigned z=a[i]&~m;z;){unsigned v=z&-z;z-=v;d[(m|v)*n+__builtin_ctz(v)]+=p;}
 }
 if(accumulate(d.end()-n,d.end(),0u)!=target)return 3;
 tested++;
 }
 cout<<"PASS independent full DP on "<<tested<<" literal targets\n";
 if(tested!=12)return 4;
}
