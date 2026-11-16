#include <bits/stdc++.h>
using namespace std;
struct Row { unsigned long long code; int r,e,h,P,Es,Os,Et,Ot; };
int main(int argc,char**argv){
 if(argc<2)return 2;
 std::ifstream in(argv[1]); string l;getline(in,l);
 map<array<int,5>,map<int,Row>> bank;
 while(getline(in,l)){
  auto c=stoull(l);
  int fs[256][8]={},ft[256][8]={},H[256]={};H[0]=1;
  int A[8],B[8];
  for(int i=0;i<8;i++){A[i]=(c>>(8*i))&255;B[i]=255^A[i]^(1<<i);fs[1<<i][i]=ft[1<<i][i]=1;}
  for(int S=1;S<256;S++)for(int i=0;i<8;i++)if(S>>i&1){
    auto T=S^(1<<i);
    if(T){for(int b=T&A[i];b;b&=b-1)fs[S][i]+=fs[T][__builtin_ctz((unsigned)b)];for(int b=T&B[i];b;b&=b-1)ft[S][i]+=ft[T][__builtin_ctz((unsigned)b)];}
    H[S]+=fs[S][i];
  }
  int P=0;for(int S=0;S<256;S++)P+=H[S]*H[255^S];
  for(int i=0;i<8;i++){
    int s[2]={},t[2]={};
    for(int S=1;S<256;S++)if(S>>i&1){int p=__builtin_parity((unsigned)S);s[p]+=fs[S][i]*H[255^S];t[p]+=ft[S][i]*H[255^S];}
    auto key=array<int,5>{H[255],P,s[0],s[1],t[0]+t[1]};
    bank[key].emplace(t[0],Row{c,i,t[0],H[255],P,s[0],s[1],t[0],t[1]});
  }
 }
 long long pairs=0;bool sample=false,unit=false;
 auto show=[](Row y){cout<<y.code<<" root "<<y.r<<" H,P,Es,Os,Ee,Oe "<<y.h<<","<<y.P<<","<<y.Es<<","<<y.Os<<","<<y.Et<<","<<y.Ot<<"\n";};
 for(auto &[key,m]:bank)for(auto &[e,x]:m)if(m.count(e+1)){
   pairs++;
   if(!sample){show(x);show(m.at(e+1));sample=true;}
   if(!unit&&x.Es-x.Os==-1){cout<<"start sign -1\n";show(x);show(m.at(e+1));unit=true;}
 }
 cout<<"profiles "<<bank.size()<<" with matched signed unit differences "<<pairs<<"\n";
}
