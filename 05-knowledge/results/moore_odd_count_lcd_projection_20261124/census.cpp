#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
vector<U> kernel(vector<U>A,int n){
 vector<int> p; int r=0;
 for(int i=0;i<n;i++){
  int j=r; while(j<int(A.size()) && !(A[j]>>i&1))j++;
  if(j==int(A.size()))continue;
  swap(A[r],A[j]);
  for(j=0;j<int(A.size());j++)if(j!=r && (A[j]>>i&1))A[j]^=A[r];
  p.push_back(i);r++;
 }
 vector<U>B;
 for(int i=0;i<n;i++)if(find(p.begin(),p.end(),i)==p.end()){
  U x=1ULL<<i;
  for(int j=0;j<r;j++)if(A[j]>>i&1)x|=1ULL<<p[j];
  B.push_back(x);
 }
 return B;
}
vector<U> from_lcd(const vector<U>&G,int n){
 int m=G.size();vector<U>A(m),I(m);
 for(int i=0;i<m;i++){
  I[i]=1ULL<<i;assert(!(__builtin_popcountll(G[i])&1));
  for(int j=0;j<m;j++)if(__builtin_popcountll(G[i]&G[j])&1)A[i]|=1ULL<<j;
 }
 for(int i=0;i<m;i++){
  int j=i;while(j<m && !(A[j]>>i&1))j++;
  assert(j<m);swap(A[i],A[j]);swap(I[i],I[j]);
  for(j=0;j<m;j++)if(j!=i&&(A[j]>>i&1)){A[j]^=A[i];I[j]^=I[i];}
 }
 vector<U>B(n);
 for(int i=0;i<n;i++){
  U g=0;for(int j=0;j<m;j++)if(G[j]>>i&1)g|=1ULL<<j;
  U p=0;for(int j=0;j<m;j++)if(__builtin_popcountll(I[j]&g)&1)p^=G[j];
  B[i]=p^(1ULL<<i)^((1ULL<<n)-1);
 }
 return B;
}
bool edge(int n,int a,int b){
 if(a==b)return false;
 if(n==5){int t=(a-b+5)%5;return t==1||t==4;}
 if(n==10){
  vector<int>M;
  for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)M.push_back((1<<i)|(1<<j));
  return !(M[a]&M[b]);
 }
 if(a>=25&&b<25)swap(a,b);
 int i=(a%25)/5,j=(b%25)/5,x=a%5,y=b%5,t=(y-x+5)%5;
 if(a/25==b/25)return i==j&&(a<25?(t==1||t==4):(t==2||t==3));
 return (y-x-i*j+25)%5==0;
}
vector<U> census(const vector<U>&B,int n){
 vector<U> H(n+1);
 U x=0;H[0]=1;
 for(U i=1;i<(1ULL<<B.size());i++){x^=B[__builtin_ctzll(i)];H[__builtin_popcountll(x)]++;}
 return H;
}
void report(const char*label,const vector<U>&B,int n){
 auto H=census(B,n);
 cout<<label<<" dim="<<B.size()<<" weights:";
 for(int i=1;i<=n;i++)if(H[i])cout<<" "<<i<<":"<<H[i];
 cout<<endl;
}
void moves(const vector<U>&A,int n,int D){
 int minoutside=n,minadded=n*n,minexcess=n*n*n,neutral=0;
 const U full=(1ULL<<n)-1;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++)for(int m=k+1;m<n;m++){
  U S=(1ULL<<i)|(1ULL<<j)|(1ULL<<k)|(1ULL<<m);
  U W=A[i]^A[j]^A[k]^A[m];int changed=__builtin_popcountll(W&~S);
  int ds=0,ss=0;
  for(int b=0;b<n;b++){
   U R=A[b]^((S>>b&1)?W:0)^((W>>b&1)?S:0);
   int d=__builtin_popcountll(R&full);
   ds+=d-D;ss+=d*d-D*D;
  }
  minoutside=min(minoutside,changed);minadded=min(minadded,ds);minexcess=min(minexcess,ss);neutral+=ss==0;
  assert(changed>=4*D-12&&ds>=16*D-54&&ss>=0);
 }
 cout<<"weight-4 transvections: min outside "<<minoutside<<", degree sum increase "<<minadded
     <<", squared degree excess "<<minexcess<<", neutral count "<<neutral<<endl;
}
int main(){
 for(int n:{5,10,50}){
  vector<U>A(n),X;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(edge(n,i,j))A[i]|=1ULL<<j;
  int d=__builtin_popcountll(A[0]);
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j)
   assert(__builtin_popcountll(A[i]&A[j])==int(!(A[i]>>j&1)));
  X=A;for(int i=0;i<n;i++)X[i]^=1ULL<<i;
  auto U1=kernel(A,n),V=kernel(X,n),W=kernel(U1,n);
  cout<<"order "<<n<<", degree "<<d<<endl;
  if(n!=5){
   assert(from_lcd(U1,n)==A);
   auto S=census(W,n);
   for(int i=1;i<d;i+=2)assert(!S[i]);
   assert(S[d]==U(n)); moves(A,n,d);
  }
  report("ker A",U1,n);report("ker (A+I)",V,n);report("im A",W,n);
 }
 vector<U> star(10,1);
 star[0]=1022;
 assert(kernel(star,10).size()==8);
 assert(from_lcd(kernel(star,10),10)==star);
 int sum=0;for(U row:star){int d=__builtin_popcountll(row);sum+=d*d;}
 assert(sum==90);
 cout<<"star also has squared-degree sum 90, but nullity 8 instead of 4"<<endl;
}
