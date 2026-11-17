#include <bits/stdc++.h>
using namespace std;
using V=vector<int>;using ll=long long;
using P=pair<int,int>; // A=2E,g
using L=tuple<int,int,int>; // h,A,g

map<P,set<int>> profiles(int n,const char* file){
 ifstream in(file); if(!in)abort();
 int N=1<<n, read=0; V a(n);
 map<P,set<int>> K;
 while(in>>a[0]){
  for(int j=1;j<n;j++)if(!(in>>a[j]))abort();
  for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(((a[i]>>j&1)+(a[j]>>i&1))!=1)abort();
  read++;
  vector<V>D(N,V(n*n));
  V h(N);h[0]=1;
  for(int i=0;i<n;i++)D[1<<i][i*n+i]=1;
  for(int s=1;s<N;s++)for(int b=0;b<n;b++)for(int j=0;j<n;j++){
   int c=D[s][b*n+j];if(!c)continue;
   for(int d=a[j]&~s;d;d&=d-1){int v=d&-d;D[s|v][b*n+__builtin_ctz(v)]+=c;}
   h[s]+=c;
  }
  for(int b=0;b<n;b++){
   int g=h[(N-1)^(1<<b)],E=0,O=0;
   for(int s=1;s<N;s++)if(s>>b&1){
    auto &v=D[s];int u=accumulate(v.begin()+b*n,v.begin()+(b+1)*n,0);
    (__builtin_popcount(s)%2?O:E)+=u*h[(N-1)^s];
   }
   if(E>0&&O-E==g)K[{2*E,g}].insert(h[N-1]);
  }
 }
 if(read!=(n==7?456:6880))abort();
 return K;
}
set<L> laws(map<P,set<int>> const& K,bool matched){
 set<L>B;
 for(auto const& [p,S]:K)for(int h:S)
  if(!matched || S.count(h+2)||S.count(h-2))B.emplace(h,p.first,p.second);
 return B;
}
int main(int argc,char**argv){
 if(argc!=3)abort();
 int limit=1<<20;
 auto A=profiles(7,argv[1]),B=profiles(8,argv[2]);
 auto a=laws(A,true),b=laws(B,true),c=laws(B,false);
 if(!includes(b.begin(),b.end(),a.begin(),a.end()))abort();
 V f(limit+1,1),prime(limit+1,1);
 for(int z=2;z<=limit;z++){int u=z/2;f[z]=f[u]+(z&1?0:u+(u&1));}
 prime[0]=prime[1]=0;
 for(int z=2;z*z<=limit;z++)if(prime[z])for(int j=z*z;j<=limit;j+=z)prime[j]=0;
 cout<<"cutoff bank laws covered/odds primes/not7 sevenprime/not21 max_zero_run mean singleton"<<endl;
 for(int bank=0;bank<4;bank++){
  set<L>D=bank==0?set<L>{{13,12,3},{15,12,3},{15,18,3},{17,18,3}}:bank==1?a:bank==2?b:c;
  V R(limit+1);
  for(auto [h,A,g]:D)for(int z=0;z<=(limit-h)/A;z++){
   int n=h+A*z+g*(f[z]-1);if(n<=limit)R[n]++;
  }
  int missing=0,last=0;
  for(int n=23;n<=limit;n+=2)if(!R[n]){missing++; last=n;}
  cout<<"prefix "<<bank<<" missing "<<missing<<" last "<<last<<endl;
  for(int k:{12,14,16,18,20}){
   int ctr[3]={},hit[3]={},mx=0,run=0,hi=0;
   ll sum=0; int L=1<<(k-1),U=1<<k;
   for(int n=L+1;n<U;n+=2)sum+=R[n];
   for(int n=L+1;n<U;n+=2){
    bool p=prime[n]&&n!=7, q=n%7==0&&prime[n/7]&&n/7!=3;
    ctr[0]++;ctr[1]+=p;ctr[2]+=q;
    if(R[n]){
     hit[0]++;hit[1]+=p;hit[2]+=q;run=0;
     hi+=R[n]==1;
    } else mx=max(mx,++run);
   }
   cout<<k<<" "<<bank<<" "<<D.size();
   for(int i=0;i<3;i++)cout<<" "<<hit[i]<<"/"<<ctr[i];
   cout<<" "<<mx<<" "<<fixed<<setprecision(6)<<double(sum)/ctr[0]<<" "<<hi<<endl;
  }
 }
}
