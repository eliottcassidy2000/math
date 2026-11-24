#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 int n=stoi(argv[1]),k=n-1,full=(1<<k)-1;
 string s;uint64_t old=~0ull, rows=0, cores=0;
 array<uint64_t,20000> tally{}, witness{};
 int end[512][9],start[512][9],Q[9][9],cv[512],cn[512][9];
 while(cin>>s){
  assert(s.size()==n*(n-1)/2);
  int out[9]={},apex=0,pos=0,shift=0;uint64_t id=0,mask=0;
  for(int i=0;i<k;i++){
   for(int j=i+1;j<k;j++){
    if(s[pos++]=='1'){out[i]|=1<<j;id|=1ull<<shift;}
    else out[j]|=1<<i,mask|=1ull<<(pos-1);
    shift++;
   }
   if(s[pos++]=='0')apex|=1<<i,mask|=1ull<<(pos-1);
  }
  if(id!=old){
   old=id;cores++;
   for(int p=0;p<k;p++)end[0][p]=start[0][p]=0;
   for(int S=1;S<=full;S++){
    int b=S;while(b){
      int i=__builtin_ctz(b);b&=b-1;
      int R=S^(1<<i),e=R==0,t=e;
      int a=R&out[i];while(a){int j=__builtin_ctz(a);a&=a-1;t+=start[R][j];}
      a=R&~out[i];while(a){int j=__builtin_ctz(a);a&=a-1;e+=end[R][j];}
      start[S][i]=t;end[S][i]=e;
    }
   }
   for(int i=0;i<k;i++)for(int j=0;j<k;j++)Q[i][j]=0;
   for(int S=1;S<full;S++){
    const int R=full^S;int a=S;
    while(a){int i=__builtin_ctz(a);a&=a-1;int b=R;
     while(b){int j=__builtin_ctz(b);b&=b-1;Q[i][j]+=end[S][i]*start[R][j];}}
   }
   cv[0]=0;int linear[9],q[9][9];
   for(int i=0;i<k;i++){
    cv[0]+=end[full][i];linear[i]=start[full][i]-end[full][i];
    for(int j=0;j<k;j++){linear[i]+=Q[j][i];q[i][j]=Q[i][j]+Q[j][i];}
    cn[0][i]=0;
   }
   for(int S=1;S<=full;S++){
    int j=__builtin_ctz(S),R=S&(S-1);
    cv[S]=cv[R]+linear[j]-cn[R][j];
    for(int i=0;i<k;i++)cn[S][i]=cn[R][i]+q[j][i];
   }
  }
  int h=cv[apex];assert((h&1)&&h>0&&h<int(tally.size()));
  if(!tally[h]++)witness[h]=mask;
  rows++;
 }
 cerr<<"order "<<n<<" rows "<<rows<<" cached cores "<<cores<<"\n";
 cout<<"count,classes,witness_backward_mask\n";
 for(int h=1;h<(int)tally.size();h+=2)if(tally[h])cout<<h<<','<<tally[h]<<','<<witness[h]<<'\n';
}
