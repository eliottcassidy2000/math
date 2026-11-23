#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>
using namespace std;
using Rows=vector<int>;

bool eligible(int p){
 if(p<5||p==7)return false;
 for(int j=2;j*j<=p;j++)if(p%j==0)return false;
 return true;
}
Rows anchor(int n,int mask){
 Rows r(n);int bit=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)
  if(j==i+1||((mask>>bit++)&1))r[i]|=1<<j;
  else r[j]|=1<<i;
 return r;
}
void flip(Rows&r,int u,int v){r[u]^=1<<v;r[v]^=1<<u;}
bool strong(const Rows&r){
 int s=1,p=1,z=0;
 while(s!=z){z=s;for(int i=0;i<(int)r.size();i++)if(s>>i&1)s|=r[i];}
 z=0;
 while(p!=z){z=p;for(int i=0;i<(int)r.size();i++)if(r[i]&p)p|=1<<i;}
 return s==((1<<r.size())-1)&&p==s;
}
struct DP{
 array<array<int,10>,1024> start;
 array<int,1024> h;
 int compute(const Rows&r){
  const int N=1<<r.size();
  h[0]=1;
  for(int s=1;s<N;s++){
   h[s]=0;
   for(int t=s;t;t&=t-1){
    int i=__builtin_ctz(t),m=s^(1<<i),c=!m;
    for(int j=r[i]&m;j;j&=j-1)c+=start[m][__builtin_ctz(j)];
    start[s][i]=c;h[s]+=c;
   }
  }
  return h[N-1];
 }
};
Rows pad(const Rows&c,int k){
 Rows r=c;int n=c.size(),F=(1<<(n+k))-1;
 for(int i=0;i<n;i++)r[i]|=F^((1<<n)-1);
 for(int j=n;j<n+k;j++)r.push_back(F^((1<<(j+1))-1));
 return r;
}
int main(){
 DP dp,aux;
 set<int> adjacent,middle;
 map<int,int> smallPrimeFrames;
 long long testedUpper=0;
 for(int n=4;n<=8;n++){
  int N=1<<n,F=N-1,qty=1<<((n-1)*(n-2)/2),upper=0;
  for(int t=0;t<qty;t++){
   Rows r=anchor(n,t);int h=dp.compute(r);
   if(h%7==0&&eligible(h/7)){
    ++upper;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
     flip(r,i,j);assert(aux.compute(r)!=h/7);flip(r,i,j);
    }
   }
   if(!eligible(h)||!strong(r))continue;
   if(h<45){assert(n<8);++smallPrimeFrames[h];}
   Rows U(n);
   for(int s=1;s<N;s++)for(int t=s;t;t&=t-1){
    int i=__builtin_ctz(t);U[i]+=dp.h[s^F]*dp.start[s][i];
   }
   Rows dual;
   if(any_of(U.begin(),U.end(),[&](int x){return 6*h<x&&x<=7*h;})){
    dual=r;for(int i=0;i<n;i++)dual[i]^=F^(1<<i);
    aux.compute(dual);
   }
   for(int i=0;i<n;i++){
    int E=dual.empty()?0:aux.start[F][i];
    for(int k:{1,2})if(k==1?(U[i]>6*h&&U[i]-E==6*h):U[i]==3*h){
      assert(n==8);
      auto& field=k==1?adjacent:middle;
      if(field.insert(h).second){
       Rows z=pad(r,k);assert(aux.compute(z)==h);
       flip(z,i,n+k-1);assert(aux.compute(z)==7*h&&strong(z));
       if(!dual.empty())aux.compute(dual); // preserve the terminal table for other roots
       cout<<n<<" "<<h<<" "<<k<<" "<<i<<" "<<U[i]<<" "<<E<<"\n";
      }
    }
   }
  }
  testedUpper+=1LL*upper*n*(n-1)/2;
  cout<<"order "<<n<<", "<<qty<<" frames, "<<upper<<" eligible upper frames\n";
 }
 assert(adjacent==set<int>{163});
 assert((middle==set<int>{73,179,191,233,307,389}));
 assert(smallPrimeFrames.size()==11);
 cout<<"global upper-frame arc tests "<<testedUpper<<"\n";
}
