#include <array>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <algorithm>
using namespace std;using ll=long long;
using V=vector<int>;
using P=array<int,3>;
struct D{
 V a; int n; array<int,2>M;
 int C[256][3][3]={};
 D(V r,array<int,2>m):a(r),n(r.size()),M(m){
  int k=1<<n,t[256][8][8]={};
  for(int i=0;i<n;i++)t[1<<i][i][i]=1;
  for(int S=1;S<k;S++)for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   int f=t[S][i][j];if(!f)continue;
   for(int b=a[j]&~S;b;b&=b-1){int x=__builtin_ctz(b);t[S|(1<<x)][i][x]+=f;}
   for(int p=0;p<3;p++)for(int q=0;q<3;q++)
    if((p==2||!(M[p]>>i&1))&&(q==2||(M[q]>>j&1)))C[S][p][q]+=f;
  }
 }
 vector<P> parts(){
  vector<P> z{{0,0,0}};
  for(int i=0;i<n;i++){vector<P>g;
   for(P v:z)for(int j=0;j<3;j++){v[j]^=1<<i;g.push_back(v);v[j]^=1<<i;}
   z.swap(g);
  }
  return z;
 }
};
void readout(D&A,D&B,int sign,int flip){
 ll W[8][8][8]={};
 auto pa=A.parts(),pb=B.parts();
 for(P a:pa)for(P b:pb)for(int v=0;v<8;v++){
   int f=0,l=0;ll h=1;
   for(int i=0;i<3;i++){
    int p=i==0?2:(i==1?sign:1-sign); // previous port; 2=none
    int q=i==2?2:(i==0?sign:1-sign); // next port
    bool m=v>>i&1;
    if(!a[i]&&!b[i]&&!m){if(i==1&&sign!=flip)h=0;continue;}
    if(a[i])h*=A.C[a[i]][p][b[i]||m?2:q];
    if(b[i])h*=B.C[b[i]][a[i]||m?2:p][q];
    if(m&&!a[i]&&p!=2)f|=1<<i;
    if(m&&!b[i]&&q!=2)l|=1<<i;
   }
   W[v][f][l]+=h;
 }
 for(int v=0;v<8;v++)for(int f=0;f<8;f++)for(int l=0;l<8;l++)if(W[v][f][l])printf("%d %d %d %lld\n",v,f,l,W[v][f][l]);
}
int main(){
 D A({2,4,1},{1,4}),B({2,12,17,37,75,87,143,63},{37,67});
 for(int j=0;j<2;j++)for(int f=0;f<2;f++)for(int s=0;s<2;s++){printf("order=%d sign=%d flip=%d\n",j,s,f);if(j==0)readout(A,B,s,f);else readout(B,A,s,f);}
}
