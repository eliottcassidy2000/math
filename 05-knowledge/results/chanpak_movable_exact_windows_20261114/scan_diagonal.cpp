#include <bits/stdc++.h>
using namespace std;
struct Rec{uint64_t rows; unsigned short U; unsigned char n;};
using Small=unordered_map<int,Rec>;
unordered_map<uint64_t,Small> tab;
void prof(int n,int *A){
 int F=(1<<n)-1, f=0,h=0,d[8]={0},M[8][8]={},E[256]={0},W[256]={0};
 int ending[256][8]={},starting[256][8]={},H[256]={};
 H[0]=1;
 for(int i=0;i<n;i++)ending[1<<i][i]=starting[1<<i][i]=1;
 for(int s=1;s<=F;s++){
   for(int i=0;i<n;i++){
     H[s]+=ending[s][i]; int vs=(~s)&F; int w=vs&A[i];
     for(int k=w;k;k&=k-1) ending[s|(k&-k)][__builtin_ctz(k)]+=ending[s][i];
     w=vs&(~A[i]);
     for(int k=w;k;k&=k-1) starting[s|(k&-k)][__builtin_ctz(k)]+=starting[s][i];
   }
 }
 h=H[F];
 for(int s=0;s<=F;s++){
    int r=F^s;f+=H[s]*H[r];
    for(int u=s;u;u&=u-1){int i=__builtin_ctz(u);d[i]+=(starting[s][i]-ending[s][i])*H[r];
      for(int v=r;v;v&=v-1){int j=__builtin_ctz(v);M[i][j]+=ending[s][i]*starting[r][j];}}
 }
 uint64_t rows=0;
 for(int i=0;i<n;i++)rows|=uint64_t(A[i])<<(8*i);
 W[0]=f-h;
 auto put=[&](int U){uint64_t key=(uint64_t(h)<<40)|(uint64_t(f)<<20)|W[U];tab[key].emplace(E[U],Rec{rows,(unsigned short)U,(unsigned char)n});};
 put(0);
 for(int s=1;s<=F;s++){
   int r=s&(s-1),i=__builtin_ctz(s);
   W[s]=W[r]+d[i]; E[s]=E[r];
   for(int u=r;u;u&=u-1){int k=__builtin_ctz(u);E[s]+=M[i][k]+M[k][i];}
   put(s);
 }
}
int main(){
 int n,A[8]={},R[8]={},cases=0;
 while(cin>>n){
   for(int i=0;i<n;i++)cin>>A[i];
   prof(n,A);cases++;
   if(n==7)for(int U=0;U<128;U++){
      for(int i=0;i<n;i++)R[i]=A[i]|((U>>i&1)?0:128);
      R[7]=U;prof(8,R);cases++;
   }
 }
 cout<<"cases\t"<<cases<<"\tgroups\t"<<tab.size()<<"\n";
 int best=0; vector<pair<uint64_t,int>> recs;
 for(auto &[key,s]:tab) for(auto &[z,r]:s)if(!s.count(z-2)){
    int p=z;while(s.count(p))p+=2;
    int c=(p-z)/2;
    if(c>best)recs.clear(),best=c;
    if(c==best)recs.emplace_back(key,z);
 }
 cout<<"best "<<best<<" groups "<<recs.size()<<"\n";
 sort(recs.begin(),recs.end()); int report=0;
 for(auto [key,z]:recs){
   cout<<"h f W: "<<(key>>40)<<' '<<((key>>20)&((1<<20)-1))<<' '<<(key&((1<<20)-1))<<"\n";
   for(int j=0;j<best;j++){
      auto r=tab[key][z+2*j];cout<<"E "<<z+2*j<<" | ";
      for(int i=0;i<r.n;i++)cout<<(int)((r.rows>>(8*i))&255)<<' ';
      cout<<"| U "<<r.U<<"\n";
   }
   if(++report==12)break;
 }
}
