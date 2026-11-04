#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using I=int64_t;
// Backward recurrence, independent of the forward endpoint discovery code.
int main(){int n;while(cin>>n){if(n<1||n>20)return 2;vector<unsigned>T(n);for(auto&x:T)cin>>x;int f=(1<<n)-1;vector<I>d((f+1)*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;
 for(int m=1;m<=f;m++)for(unsigned z=m;z;z&=z-1){int i=countr_zero(z),r=m^(1<<i);if(!r)continue;for(unsigned q=T[i]&r;q;q&=q-1)d[m*n+i]+=d[r*n+countr_zero(q)];}
 I H=0;for(int i=0;i<n;i++)H+=d[f*n+i];cout<<H<<'\n'<<flush;}return 0;}
