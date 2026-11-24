#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;
int md(int x){return (x%5+5)%5;}
int main(){
 bool A[50][50]={};
 vector<pair<int,int>> ed;
 for(int x=0;x<50;x++)for(int y=x+1;y<50;y++){
  int t=x/25,u=y/25,i=x/5%5,j=y/5%5,r=x%5,s=y%5;
  bool z=t==u?(i==j&&(md(r-s)==t+1||md(s-r)==t+1)):(s==md(i*j+r));
  if(z){A[x][y]=A[y][x]=1;ed.push_back({x,y});}
 }
 assert(ed.size()==175);
 int minD=1000,maxD=0,num=0,no=0,fail=0;
 vector<int> P={2,3,4},res; 
 do for(int q=1;q<=2;q++){
  vector<int> Q;for(int j=1;j<5;j++)if(j!=q)Q.push_back(j);
  do for(int f=0;f<1024;f++){
   int pb[5]={0,1,P[0],P[1],P[2]}, qb[5]={0,q,Q[0],Q[1],Q[2]};
   int p[50]={},z[50],len=0;
   for(int j=0;j<5;j++)for(int t=0;t<2;t++){
    int a=0,m=t?25+5*qb[j]:5*pb[j];
    if(len)while(!A[p[len-1]][m+a])a++;
    assert(a<5);
    int s=(t+1)*((f>>(2*j+t)&1)?-1:1);
    for(int k=0;k<5;k++){p[len]=m+md(a+s*k);z[p[len]]=len;len++;}
   }
   int d=0;
   for(auto [u,v]:ed)d+=2*!A[p[49-z[u]]][p[49-z[v]]];
   if(d<minD){minD=d;res.assign(p,p+50);}
   maxD=max(maxD,d);num++;
   int S[50],freq[50]={},pos[50];
   for(int i=0;i<50;i++)S[i]=49-i;
   for(auto [u,v]:ed){
    int i=z[u],j=z[v];if(i>j)swap(i,j);
    if(j>i+1){S[i]--;S[j]++;}
   }
   for(int i=0;i<50;i++){freq[S[i]]++;pos[S[i]]=i;}
   auto arc=[&](int u,int v){return u<v ? (v==u+1||!A[p[u]][p[v]]) : (u!=v+1&&A[p[u]][p[v]]);};
   bool imp=false;
   for(int i=0;i<50;i++)if(freq[S[i]]==1&&freq[49-S[i]]==1)
    for(int j=i+1;j<50;j++)if(freq[S[j]]==1&&freq[49-S[j]]==1)
     if(arc(i,j)==arc(pos[49-S[i]],pos[49-S[j]]))imp=true;
   no+=imp;fail+=!imp;
  }while(next_permutation(Q.begin(),Q.end()));
 }while(next_permutation(P.begin(),P.end()));
 cout<<num<<" "<<minD<<" "<<maxD<<" "<<no<<" "<<fail<<endl;
 assert(num==73728 && minD==124 && maxD==204 && no==73728 && fail==0);
 for(auto x:res)cout<<x<<" ";cout<<endl;
}
