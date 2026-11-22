#include<bits/stdc++.h>
using namespace std;using I=__int128_t;
string repr(I n){string r;do{r+='0'+n%10;n/=10;}while(n);reverse(r.begin(),r.end());return r;}
I E=(I(1)<<60)*((I(1)<<61)-1), N=E-1;
I c,h,K[2],visited[12],leaves[12],hit;
int R,root_k;
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
