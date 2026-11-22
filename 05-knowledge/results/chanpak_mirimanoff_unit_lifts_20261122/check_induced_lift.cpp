#include<bits/stdc++.h>
using namespace std;using I=long long;using A=vector<int>;
I h(const A&a) {
 int n=a.size(),M=1<<n;static I dp[4096][12];
 memset(dp,0,(M)*sizeof dp[0]);
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int S=0;S<M-1;S++)for(int i=0;i<n;i++)if(dp[S][i]){
  for(int e=a[i]&~S;e;e&=e-1) {int j=__builtin_ctz(e);dp[S|1<<j][j]+=dp[S][i];}
 }
 return accumulate(dp[M-1],dp[M-1]+n,0LL);
}
A add(A a,int x){
 int n=a.size();for(int j=0;j<n;j++)if(!(x>>j&1))a[j]|=1<<n;
 a.push_back(x);return a;
}
void run(int p, A a,int i,int j,int extra){
 int n=a.size();I T=1LL<<p,E=T*(2*T-1); 
 cout<<"START "<<p<<" edge "<<i<<" "<<j<<" rows";for(auto x:a)cout<<" "<<x;cout<<endl;
 auto test=[&](A b,vector<int> addr){
  if(h(b)!=E-1)return;
  b[i]^=1<<j;b[j]^=1<<i;
  if(h(b)!=E+1)return;
  cout<<"EXACT "<<p<<" add";for(int x:addr)cout<<' '<<x;cout<<" upper";for(int x:b)cout<<' '<<x;cout<<endl;
 };
 for(int x=0;x<(1<<n);x++){
  A b=add(a,x);
  if(extra==1)test(b,{x});
  else if(h(b)<=E-1) for(int y=0;y<(2<<n);y++){
   A c=add(b,y);
   if(extra==2)test(c,{x,y});
   else if(h(c)<=E-1) for(int z=0;z<(4<<n);z++)test(add(c,z),{x,y,z});
  }
 }
}

int main(int argc,char** argv){
 int ex=atoi(argv[1]),p=atoi(argv[2]);int n=2*p;
 if(p==5)n=8;
 if(p==7)n=10;
 if(argc<4){
  if(p==3)run(3,{30,52,56,50,32,1},1,5,ex);
  if(p==5)run(5,{160,145,67,7,13,30,59,124},0,2,ex);
  if(p==7)run(7,{988,433,74,34,236,645,938,780,60,286},8,9,ex);
  return 0;
 }
 n=atoi(argv[3]);string deck=argc>=5?argv[4]:"session212/cycle-levels";ifstream f(deck+"/classes"+to_string(n)+".rows");if(!f)throw runtime_error("missing class deck");
 I m=1LL<<(p-1),E=m*(2*m-1);A a(n);
 while(f>>a[0]){
  for(int i=1;i<n;i++)if(!(f>>a[i]))throw runtime_error("truncated class deck");
  if(h(a)!=E-1)continue;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   A b=a;b[i]^=1<<j;b[j]^=1<<i;
   if(h(b)==E+1){run(p,a,i,j,ex);if(ex==3)goto nextclass;}
  }
  nextclass: ;
 }
}
