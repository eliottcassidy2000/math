#include <bits/stdc++.h>
using namespace std;
using U=uint64_t;
const int n=14,N=1<<n;
static U D[N][n];
int norm(int z) {int b=z;for(int i=0;i<6;i++){z=(z>>1)|((z&1)<<6);b=min(z,b);}return b;}
U test(int Sa,int Sb,int C,bool reverse=false){
 int A[n]={};
 for(int a=0;a<7;a++)for(int b=0;b<7;b++){
  int d=(b-a+7)%7;
  if(Sa>>d&1) A[a]|=1<<b;
  if(Sb>>d&1) A[7+a]|=1<<(7+b);
  if(C>>d&1) A[a]|=1<<(7+b); else A[7+b]|=1<<a;
 }
 if(reverse){
  int B[n]={};
  for(int v=0;v<n;v++)for(int w=0;w<n;w++)if(A[v]>>w&1)B[w]|=1<<v;
  copy(B,B+n,A);
 }
 memset(D,0,sizeof(D));
 for(int i=0;i<n;i++)D[1<<i][i]=1;
 for(int s=1;s<N;s++)for(int v=0;v<n;v++){U h=D[s][v];if(!h)continue;
  for(int c=A[v]&~s;c;c&=c-1){int j=__builtin_ctz(c);D[s|(1<<j)][j]+=h;}
 }
 return accumulate(D[N-1],D[N-1]+n,U(0));
}
int main(){
 vector<int> types;for(int s=0;s<8;s++){
  int m=0;for(int k=1;k<=3;k++)m|=1<<((s>>(k-1)&1)?k:7-k); types.push_back(m);}
 int c0=0;for(int k:{1,2,3})c0|=1<<k;
 int c1=0;for(int k:{1,2,4})c1|=1<<k;
 U best=UINT64_MAX,next=UINT64_MAX;int chk=0;vector<tuple<int,int,int,U>> results;
 for(int a:{c0,c1}) for(int b:types) for(int c=1;c<127;c++)if(norm(c)==c){
  auto h=test(a,b,c);assert(h%7==0 && h==test(a,b,c,true));chk++;results.push_back({a,b,c,h});if(h<best)best=h;
 }
 int nb=0;for(auto[a,b,c,h]:results){if(h==best){nb++;cout<<a<<" "<<b<<" "<<c<<" "<<h<<" "<<h/7<<"\n";}
 else next=min(next,h);}
 assert(chk==288 && best==4698337);
 cout<<"checked "<<chk<<", ties "<<nb<<", best "<<best<<", next "<<next<<"\n";
}
