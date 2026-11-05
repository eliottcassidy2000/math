#include <iostream>
#include <vector>
#include <array>
using namespace std;
int main(){int count;cin>>count;int checked=0;for(int z=0;z<count;z++){int n,target;cin>>n>>target;vector<int>r(n);for(int&x:r)cin>>x;long long H=0,C=0;for(int M=1;M<(1<<n);M++){array<long long,14>a{},b{},c{},d{};for(int i=0;i<n;i++)if(M>>i&1)a[i]=1;if(M&1)c[0]=1;long long total=0;for(int step=1;step<n;step++){b.fill(0);d.fill(0);for(int i=0;i<n;i++)if(M>>i&1){int rem=r[i]&M;while(rem){int bit=rem&-rem;rem-=bit;int j=__builtin_ctz(bit);b[j]+=a[i];d[j]+=c[i];}}a=b;c=d;}for(int i=0;i<n;i++)total+=a[i];int sign=((n-__builtin_popcount((unsigned)M))&1)?-1:1;H+=sign*total;if(M&1){long long closed=0;for(int i=0;i<n;i++)if(r[i]&1)closed+=c[i];C+=sign*closed;}}if(H!=target||C!=2){cerr<<"failure "<<target<<" got "<<H<<","<<C<<endl;return 1;}checked++;}cout<<"{\"status\":\"PASS\",\"independent_walk_inclusion_exclusion_graphs\":"<<checked<<",\"every_cycle_count\":2}\n";}
