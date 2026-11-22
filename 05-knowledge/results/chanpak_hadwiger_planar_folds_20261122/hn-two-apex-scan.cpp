#include<bits/stdc++.h>
using namespace std;using I=__int128_t;using U=uint64_t;
I order_count(const vector<pair<int,int>>& tags,int turn,bool arc){
 I a[2][3][2]{};
 a[1][0][0]=1;
 for(auto [u,v]:tags){
  if(turn)swap(u,v);
  I b[2][3][2]{};
  for(int i=0;i<2;i++)for(int j=0;j<3;j++)for(int k=0;k<2;k++)if(auto d=a[i][j][k]){
   b[!u][j][k]+=d;
   if(j||u)b[i][v?1:2][k]+=d;
   if(k||v)b[i][j][1]+=d;
  }
  memcpy(a,b,sizeof(a));
 }
 I s=0;for(int k=0;k<2;k++){s+=a[1][2][k];if(arc)s+=a[1][0][k];}
 return s;
}
I H(const vector<pair<int,int>>& tags,bool arc){
 return order_count(tags,0,arc)+order_count(tags,1,!arc);
}
vector<pair<int,int>> word(U a,U b,int L,int R){
 vector<pair<int,int>> z;
 U k=a-L,r=b-R;
 if(k){
  for(int i=63-__builtin_clzll(k);i>=0;i--){int t=k>>i&1;z.emplace_back(t,t);}
 }
 z.emplace_back(L,L);
 z.emplace_back(1,0);z.emplace_back(0,1);
 z.emplace_back(!R,!R);
 if(r){
  for(int i=0;i<=63-__builtin_clzll(r);i++){int t=!(r>>i&1);z.emplace_back(t,t);}
 }
 return z;
}
int main(int argc,char**argv){
 for(int j=1;j<argc;j++){
  int p=atoi(argv[j]);if(p<2||p>30)return 2;I t=I(1)<<(p-1),N=t*(2*t-1)-1;
  int hits=0; U tested=0;
  for(U a=0;;a++){
   U b=a+1;
   if((I)8*a*b+1>N)break;
   ++tested;
   for(int L=0;L<=1&&(U)L<=a;L++)for(int R=0;R<=1&&(U)R<=b;R++){
    auto z=word(a,b,L,R);
    // A→B is the lower orientation by the weighted paired-cut identity for b=a+1.
    I u=H(z,1);
    if(u!=N)continue;
    assert(H(z,0)==N+2);
    cout<<"HIT "<<p<<" a "<<a<<" lr "<<L<<R<<" tags ";
    for(auto [x,y]:z)cout<<x<<y<<" ";
    cout<<endl;++hits;
   }
  }
  cout<<"exp "<<p<<" addresses "<<tested<<" hits "<<hits<<endl;
 }
}
