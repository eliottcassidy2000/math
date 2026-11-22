#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
using namespace std;using U=unsigned __int128;
static U pows[123],goal;static int N;static int vals[12];static uint64_t tested;
void walk(int q,int width,U P,U H,int l){
 if(q==2) {
  U base=H+(P+1)*pows[width]+2;
  if(base<=goal)return;
  U K=base-goal;uint64_t K0=uint64_t(K);
  int m=K0?__builtin_ctzll(K0):64+__builtin_ctzll(uint64_t(K>>64));
  auto e=[&](int a){
   if(a<1||a>=width)return;
   int b=width-a;
   if(l>0&&vals[0]>b)return;
   if(l==0&&a>b)return;
   tested++;
   if((P+1)*pows[a]+pows[b]==K){
    vals[l]=a;vals[l+1]=b;
    cout<<"HIT";for(int i=0;i<=l+1;i++)cout<<' '<<vals[i];cout<<endl;
   }
  };
  if(m<width/2 || (width%2 && m==width/2)) {e(m);e(width-m);}
  if(width%2==0 && m>width/2)e(width/2);
  return;
 }
 for(int a=1;a<=width-q+1;a++){
  vals[l]=a;
  if(l%2==0)walk(q-1,width-a,(P+1)*pows[a]-1,H,l+1);
  else{
   U K=H+P*(pows[a]-1)+1;
   if(K<goal)walk(q-1,width-a,P*pows[a],K,l+1);
  }
 }
}
int main(int argc,char**argv){
 const int p=stoi(argv[1]);N=2*p;int limit=stoi(argv[2]);
 pows[0]=1;for(int i=1;i<123;i++)pows[i]=pows[i-1]<<1;
 goal=pows[2*p-1]-pows[p-1]-1;
 for(int m=2;m<=limit;m+=2){
  tested=0;walk(m,N,0,1,0);cout<<"runs "<<m<<" tested "<<tested<<endl;
 }
}
