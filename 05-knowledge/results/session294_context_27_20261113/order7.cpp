#include <cstdio>
#include <algorithm>
#include <array>
#include <vector>
int main(){int n=7,m=n*(n-1)/2; std::vector<std::array<int,7>> triples,p5,p7; for(int k:{3,5,7})for(int b=0;b<128;b++)if(__builtin_popcount((unsigned)b)==k){std::array<int,7>x{};int j=0;for(int i=0;i<n;i++)if(b>>i&1)x[j++]=i;do{(k==3?triples:k==5?p5:p7).push_back(x);}while(std::next_permutation(x.begin()+1,x.begin()+k));}
long triangle=0, Bpass=0, hit=0;auto cy=[](int *a,auto x,int k){for(int j=0;j<k;j++)if(!((a[x[j]]>>x[(j+1)%k])&1))return false;return true;};for(int mask=0;mask<1<<m;mask++){int a[7]={};for(int i=0,k=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(mask>>k&1)a[i]|=1<<j;else a[j]|=1<<i;int t[35],nt=0;for(auto x:triples)if(cy(a,x,3))t[nt++]=1<<x[0]|1<<x[1]|1<<x[2];int A=2*nt;if(A%9!=4)continue;triangle++;int B=0;for(auto x:p5)if(cy(a,x,5))B+=2;for(int i=0;i<nt;i++)for(int j=0;j<i;j++)if(!(t[i]&t[j]))B+=4;if(B%9!=4)continue;Bpass++;int C=0;for(auto x:p7)if(cy(a,x,7))C+=2;if((1+A+B+C)%27==0){hit++;if(hit==1){for(int x:a)printf("%d ",x);printf(" H=%d\n",1+A+B+C);}}}printf("%ld %ld %ld\n",triangle,Bpass,hit);}
