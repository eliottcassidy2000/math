// Independent census: enumerate every receiver subset in binary-reflected Gray order.
// No producer DP, receipt, or passing-set table is read.
#include <array>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
int plusq(int a,int b){return ((((a>>4)+(b>>4))&127)<<4)|((a^b)&15);}
int minusq(int a,int b){return ((((a>>4)-(b>>4))&127)<<4)|((a^b)&15);}
int main(){
 array<int,20> sid{},lab{};array<int,8> oid{},ol{};
 for(int i=0;i<20;i++)if(!(cin>>sid[i]>>lab[i]))throw runtime_error("cheap input");
 for(int i=0;i<8;i++)if(!(cin>>oid[i]>>ol[i]))throw runtime_error("outside input");
 int total=0;
 for(int a=0;a<17;a++)for(int b=a+1;b<18;b++)for(int c=b+1;c<19;c++)for(int d=c+1;d<20;d++){
  array<array<bool,2048>,17> reach{};array<int,16> L{};int j=0;
  for(int i=0;i<20;i++)if(i!=a&&i!=b&&i!=c&&i!=d)L[j++]=lab[i];
  if(j!=16)throw runtime_error("receiver size");
  int q=0,k=0,gray=0;reach[0][0]=true;
  for(unsigned int z=1;z<65536;z++){
   int bit=__builtin_ctz(z);gray^=1<<bit;
   if(gray&(1<<bit)){q=plusq(q,L[bit]);k++;}else{q=minusq(q,L[bit]);k--;}
   reach[k][q]=true;
  }
  int pass=0;array<int,8> wg{},wq{};array<int,8> nh{};
  for(int p=0;p<8;p++){
   for(int r: {8,9})for(int x=0;x<2048;x++)if(!reach[r][x]&&!reach[r-1][minusq(x,ol[p])]){
    if(nh[p]++==0){wg[p]=r;wq[p]=x;}
   }
   if(!nh[p])pass|=1<<p;
  }
  cout<<"R "<<sid[a]<<' '<<sid[b]<<' '<<sid[c]<<' '<<sid[d]<<' '<<pass;
  for(int p=0;p<8;p++)cout<<' '<<wg[p]<<' '<<wq[p]<<' '<<nh[p];cout<<'\n';
  for(int p=0;p<8;p++)for(int t=p+1;t<8;t++){
   if((pass&(1<<p))||(pass&(1<<t)))continue; // proved containment implies grades 8..10.
   int holes=0,fg=0,fq=0;
   for(int r: {8,9,10})for(int x=0;x<2048;x++){
    bool ok=reach[r][x]||reach[r-1][minusq(x,ol[p])]||reach[r-1][minusq(x,ol[t])]||reach[r-2][minusq(x,plusq(ol[p],ol[t]))];
    if(!ok){if(holes++==0){fg=r;fq=x;}}
   }
   if(holes)cout<<"F "<<sid[a]<<' '<<sid[b]<<' '<<sid[c]<<' '<<sid[d]<<' '<<oid[p]<<' '<<oid[t]<<' '<<fg<<' '<<fq<<' '<<holes<<'\n';
  }
  total++;
 }
 if(total!=4845)throw runtime_error("freeze count");
}
