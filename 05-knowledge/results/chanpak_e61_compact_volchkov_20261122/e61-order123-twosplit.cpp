#include<bits/stdc++.h>
using namespace std;using I=__int128_t;
string repr(I n){string r;do{r+='0'+n%10;n/=10;}while(n);reverse(r.begin(),r.end());return r;}
I E=(I(1)<<60)*((I(1)<<61)-1), N=E-1;
I c,h,K[2],visited[12],leaves[12],hit;
int R,root_k;
int tz(I a){
 auto l=(uint64_t)a;
 return l?__builtin_ctzll(l):64+__builtin_ctzll((uint64_t)(a>>64));
}
void split_last(int L,int bit,I x,I f){
 if(L<2)return;
 I B=I(1)<<L, U,V,D;
 if(bit==0){
  U=h*x;V=c;
  D=N-(K[1]+c*(B*x-1)+h*(f-x+(x&1)));
 }else{
  U=h*(x+1);V=c+h;
  D=K[0]+c*B*(x+1)+h*(f+B*(x+1)+2)-N;
 }
 if(D<=0)return;
 int a=tz(U),b=tz(V),d=tz(D);
 int q[3]={d-a,L+b-d,(L+b-a)%2==0?(L+b-a)/2:0};
 for(int j=0;j<3;j++){
  int u=q[j];
  if(u<=0||u>=L||(j>=1&&q[0]==u)||(j==2&&q[1]==u))continue;
  ++leaves[R];
  if((U<<u)+(V<<(L-u))!=D)continue;
  I y=bit?B*(x+1)-(I(1)<<(L-u)):B*x+(I(1)<<(L-u))-1;
  ++hit;cout<<"HIT "<<repr(y)<<endl;
 }
}
void go(int left,int bit,int len,I y,I f){
 int par=y&1;
 if(left==0){
  ++leaves[len];
  if(K[par]+c*y+h*f==N){cout<<"HIT "<<repr(y)<<endl;++hit;}
  return;
 }
 if(len>=R)return;
 I min=K[0]+c*(y<<left)+h*f;
 I max=K[1]+c*(((y+1)<<left)-1)+h*(f+(y+1)*((I(1)<<left)-1)+1);
 if(N<min||N>max)return;
 ++visited[len];
 if(len==R-2){
  I B=I(1)<<left;
  go(0,bit^1,len+1,bit?B*y+B-1:B*y,bit?f:f+(B-1)*y+par);
  split_last(left,bit,y,f);
  return;
 }
 for(int s=1;s<=left;s++){
  if(len+1==R&&s<left)continue;
  I B=I(1)<<s;
  go(left-s,bit^1,len+1,bit?B*y+B-1:B*y,bit?f:f+(B-1)*y+par);
 }
}
int main(int argc,char**argv){
 if(argc!=2||(R=atoi(argv[1]))<1||R>9)return 2;
 const array<I,5> charts[]={{81,25,43,68,115},{411,75,367,392,113},{1021,125,915,940,112}};
 for(auto [c0,h0,k0,k1,w]:charts){
  c=c0;h=h0;K[0]=k0;K[1]=k1;root_k=w;hit=0;
  fill(begin(visited),end(visited),0);fill(begin(leaves),end(leaves),0);
  if(gcd((long long)c,(long long)h)>1)continue;
  for(int s=1;s<=w;s++)go(w-s,0,1,(I(1)<<s)-1,0);
  cout<<"chart "<<repr(c)<<" width "<<root_k<<" at most "<<R<<" runs hits "<<repr(hit)<<endl;
  for(int r=1;r<=R;r++)cout<<r<<": "<<repr(visited[r])<<" nodes, "<<repr(leaves[r])<<" leaves\n";
 }
}
