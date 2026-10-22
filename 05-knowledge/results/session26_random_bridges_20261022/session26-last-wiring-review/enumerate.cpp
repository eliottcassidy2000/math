#include <iostream>
#include <vector>
using namespace std;
struct Tri {int m,p,q,v;};
int main(){for(int n=3;n<=8;n++){int pos[8][8];int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)pos[i][j]=(j==i+1)?-1:k++;
vector<Tri> ts;
for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int l=j+1;l<n;l++){
int bits[3]={pos[i][j],pos[j][l],pos[i][l]};int m=0,p=0,q=0;bool validq=true;
for(int z=0;z<3;z++){if(bits[z]>=0){m|=1<<bits[z];if(z<2)p|=1<<bits[z];else q|=1<<bits[z];}else if(z<2)validq=false;}
ts.push_back({m,p,validq?q:-1,(1<<i)|(1<<j)|(1<<l)});}
for(int c=0;c<(1<<k);c++){int ap=(1<<n)-1;for(auto t:ts){int z=c&t.m;if(z==t.p||z==t.q)ap&=t.v;if(!ap)break;}if(ap)cout<<n<<" "<<c<<" "<<ap<<"\n";}}}
