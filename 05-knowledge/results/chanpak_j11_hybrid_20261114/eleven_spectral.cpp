#include <array>
#include <iostream>
using namespace std;
using I=__int128;
using Mat=array<long long,9>;
static Mat A[]={
{0,0,1,1,0,0,1,1,1},{1,0,1,0,0,0,2,2,2},{1,0,1,0,0,0,4,4,4},
{2,1,3,1,1,0,3,3,4},{1,0,1,0,0,0,8,8,8},{4,3,5,1,1,0,7,7,8},
{3,2,3,0,0,0,6,6,8},{4,3,7,1,1,0,5,5,8},{8,7,9,1,1,0,15,15,16},
{8,7,11,1,1,0,11,11,16},{8,7,15,1,1,0,9,9,16}};
int ds[]={1,0,0,4,0,8,0,8,16,16,16};
int st[10],count[10]={}; long long bad[10]={}; bool stopped=false;
void rec(Mat M,I D,int k,int maxk){
 if(k){
  I T=M[0]+M[4]+M[8];
  I S=I(M[0])*M[4]-I(M[1])*M[3]+I(M[0])*M[8]-I(M[2])*M[6]+I(M[4])*M[8]-I(M[5])*M[7];
  if((1-T+S-D)>=0 || (-1-T-S-D)>=0){
   bad[k]++;
   if(!stopped){stopped=true; for(int i=0;i<k;i++)cout<<st[i]<<' ';cout<<"found T "<<(long long)T<<" S "<<(long long)S<<" D "<<(long long)D<<endl;}
  }
  count[k]++;
 }
 if(k>=maxk)return;
 for(int y=0;y<11;y++){
  Mat N{};auto B=A[y];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++) N[3*i+j]=B[3*i]*M[j]+B[3*i+1]*M[3+j]+B[3*i+2]*M[6+j];
  st[k]=y;rec(N,D*ds[y],k+1,maxk);
 }
}
int main(int n,char**s){int k=atoi(s[1]);rec({1,0,0,0,1,0,0,0,1},1,0,k);for(int i=1;i<=k;i++)cout<<i<<' '<<count[i]<<' '<<bad[i]<<endl;}
