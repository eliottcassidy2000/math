#include <bits/stdc++.h>
#include <gmpxx.h>
using namespace std;using B=mpz_class;using U=unsigned;using I=long;
const int n=10;const U full=1023,port=(1<<1)|(1<<9),ordinary=full^port;
U row[n];I t[4];B E=(B(1)<<60)*((B(1)<<61)-1);
bool cap(int p){
 vector<U> mark;
 for(U s=row[p];s;s&=s-1){
  int v=__builtin_ctz(s);U in=(full^row[p]^(1<<p));
  for(U m=row[v]&in;m;m&=m-1){
   U x=((1<<v)|(m&-m))&ordinary;
   if(!x)return false;
   for(auto y:mark)if(!(x&y))return false;
   mark.push_back(x);
  }
 }
 return true;
}
void walk(U todo,int a,int b,U next){
 if(!todo&&a&&b)t[2*(a-1)+(b-1)]++;
 U av=next&(todo|(a<2 ? 1<<1:0)|(b<2 ? 1<<9:0));
 while(av){
  U m=av&-av;av^=m;int i=__builtin_ctz(m);
  walk(todo&~m,a+(i==1),b+(i==9),row[i]);
 }
}
vector<pair<int,int>> edges;
I screened=0,admissible=0,hit=0,scalars=0;
void tryshape(){
 screened++;
 if(!cap(1)||!cap(9))return;
 admissible++;
 fill(t,t+4,0);
 walk(ordinary,0,0,full);
 I a=t[3],b=t[2]-2*a,c=t[1]-2*a,d=t[0]-2*t[1]-2*t[2]+4*a;
 for(int i=1;i<=121;i++){
  B X=B(1)<<i,K=a*X+c,N=E-1-d-b*X,Y;
  if(!K){if(N)continue;Y=2;}
  else{if(N%K!=0)continue;Y=N/K;}
  if(Y<2||mpz_popcount(Y.get_mpz_t())!=1)continue;
  scalars++;
  int j=mpz_sizeinbase(Y.get_mpz_t(),2)-1;
  cout<<"lower i,j "<<i<<" "<<j<<" poly "<<a<<" "<<b<<" "<<c<<" "<<d<<" rows";
  for(U z:row)cout<<" "<<z;
  cout<<endl;
  for(auto [u,v]:edges){
    row[u]^=1<<v;row[v]^=1<<u;
    if(cap(1)&&cap(9)){
      fill(t,t+4,0);walk(ordinary,0,0,full);
      if(t[0]+t[1]*(Y-2)+t[2]*(X-2)+t[3]*(Y-2)*(X-2)==E+1){
        cout<<"EXACT FLIP "<<u<<" "<<v<<endl;hit++;
      }
    }
    row[u]^=1<<v;row[v]^=1<<u;
  }
 }
}
void rad(int k,int s){tryshape();if(!k)return;for(int j=s;j<(int)edges.size();j++){
auto[u,v]=edges[j];row[u]^=1<<v;row[v]^=1<<u;rad(k-1,j+1);row[u]^=1<<v;row[v]^=1<<u;}
}
int main(int argc,char**argv){
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(((1<<i)|(1<<j))&port?false:true)edges.emplace_back(i,j);
 int r=argc>1?atoi(argv[1]):3;
 for(array<U,n> a: {array<U,n>{990,1020,1016,976,928,137,688,768,608,32},
 array<U,n>{998,1020,1016,977,929,456,656,768,608,32}}){
  copy(a.begin(),a.end(),row);rad(r,0);
  cerr<<screened<<" screened "<<admissible<<" cap-two "<<scalars<<" lower "<<hit<<" exact arcs"<<endl;
 }
}
